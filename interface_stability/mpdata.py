# coding: utf-8
# Copyright (c) Yifei Mo Group @ University of Maryland, College Park
# Distributed under the terms of the MIT License.

"""
All access to the Materials Project (MP) database goes through this module.

Entries are fetched with the next-gen MP API (``mp-api``). The API key is read
from the ``MP_API_KEY`` environment variable, or from ``PMG_MAPI_KEY`` in
``~/.config/.pmgrc.yaml`` / ``~/.pmgrc.yaml``.

Fetched entries are cached in memory for the life of the process. If a cache
directory is configured (``IFS_CACHE_DIR`` environment variable, or
``PMG_PD_PRELOAD_PATH`` in the pymatgen settings file), they are also cached on
disk as JSON so later runs can work offline. Cached data will not pick up later
changes to the MP database; delete the cache files to refresh.
"""

import json
import os
import uuid

from monty.json import MontyDecoder, MontyEncoder
from pymatgen.core import SETTINGS, Element
from pymatgen.entries.computed_entries import ComputedEntry

# Thermo type used for all MP queries. "GGA_GGA+U" gives entries with the
# composition-based MaterialsProject2020Compatibility corrections, so energies
# fetched for different chemical systems are on the same scale. This matches
# the method used in the papers cited in the README. "GGA_GGA+U_R2SCAN" is the
# current MP website default, but its energies depend on the chemical system
# that was queried.
DEFAULT_THERMO_TYPE = "GGA_GGA+U"
THERMO_TYPES = ("GGA_GGA+U", "GGA_GGA+U_R2SCAN", "R2SCAN")

# Thermo types whose energies depend on the chemical system that was queried: mp-api puts mixed
# GGA/GGA+U/r2SCAN entries on the energy scale of MP's phase diagram of the whole queried system.
# Entries of a subsystem are not taken from a cached larger system for these. Older mp-api versions
# return entries of different subsystems on different energy scales.
CHEMSYS_DEPENDENT_THERMO_TYPES = ("GGA_GGA+U_R2SCAN",)
MIN_MP_API_VERSION_FOR_MIXED = "0.46.5"

_thermo_type = DEFAULT_THERMO_TYPE
_memory_cache = {}


class MPDataError(RuntimeError):
    """
    The Materials Project could not be reached, or it rejected a query.
    """


def set_thermo_type(thermo_type):
    """
    Set the MP thermo type used by all later queries in this process.
    """
    global _thermo_type
    if thermo_type not in THERMO_TYPES:
        raise ValueError("Unknown thermo type {}; choose from {}".format(thermo_type, ", ".join(THERMO_TYPES)))
    if thermo_type in CHEMSYS_DEPENDENT_THERMO_TYPES:
        _check_mp_api_version(MIN_MP_API_VERSION_FOR_MIXED, thermo_type)
    _thermo_type = thermo_type


def _check_mp_api_version(min_version, thermo_type):
    """
    Raise an error if an mp-api older than min_version is installed. A missing mp-api is
    reported when a query is made.
    """
    from importlib.metadata import PackageNotFoundError, version
    from packaging.version import Version

    try:
        installed = version("mp-api")
    except PackageNotFoundError:
        return
    if Version(installed) < Version(min_version):
        raise RuntimeError("The {} thermo type needs mp-api {} or newer to put all entries on one energy scale, "
                           "but mp-api {} is installed. Upgrade it with: pip install -U mp-api"
                           .format(thermo_type, min_version, installed))


def get_thermo_type():
    return _thermo_type


def get_api_key():
    return os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY")


def get_cache_dir():
    return os.environ.get("IFS_CACHE_DIR") or SETTINGS.get("PMG_PD_PRELOAD_PATH")


def _mp_client():
    """
    The mp-api MPRester class and the error class it raises.
    """
    try:
        from mp_api.client import MPRester
        from mp_api.client.core import MPRestError
    except ImportError:
        raise ImportError("The mp-api package is required to fetch data from the Materials Project. "
                          "Install it with: pip install mp-api")
    return MPRester, MPRestError


def _query(method, *args, **kwargs):
    """
    Call an MPRester method with the configured API key. Errors from the Materials Project
    (e.g. a rejected key or no network) are raised as MPDataError.
    """
    import requests

    MPRester, MPRestError = _mp_client()
    api_key = get_api_key()
    if not api_key:
        raise RuntimeError("No Materials Project API key found. Set the MP_API_KEY environment variable, "
                           "or put PMG_MAPI_KEY in your pymatgen settings file (~/.pmgrc.yaml).")
    try:
        with MPRester(api_key) as m:
            return getattr(m, method)(*args, **kwargs)
    except (MPRestError, requests.exceptions.RequestException) as err:
        raise MPDataError("Materials Project query failed: {}".format(err)) from err


def _normalize_chemsys(chemsys):
    return sorted({Element(el).symbol if not isinstance(el, Element) else el.symbol for el in chemsys})


def _cache_path(elements, thermo_type):
    cache_dir = get_cache_dir()
    if not cache_dir:
        return None
    os.makedirs(cache_dir, exist_ok=True)
    name = "{}_{}_Entries.json".format("_".join(elements), thermo_type.replace("+", "p"))
    return os.path.join(cache_dir, name)


def _json_safe(obj):
    """
    Return a copy of obj with dict keys that JSON cannot encode (e.g. the Element keys
    of MP's "oxidation_states" data) turned into strings.
    """
    if isinstance(obj, dict):
        return {k if isinstance(k, (str, int, float, bool)) or k is None else str(k): _json_safe(v)
                for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


def _write_cache(path, entries):
    dicts = []
    for entry in entries:
        data = entry.data
        entry.data = _json_safe(data)
        try:
            dicts.append(entry.as_dict())
        finally:
            entry.data = data
    # Write to a temporary file first so a failed write never leaves a truncated cache. Its name is
    # unique, so several processes can write the same cache file at once.
    tmp_path = "{}.{}.tmp".format(path, uuid.uuid4().hex)
    try:
        with open(tmp_path, "w") as f:
            json.dump(dicts, f, cls=MontyEncoder)
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _read_cache(path):
    """
    The entries in a cache file, or None if the file cannot be used (e.g. it is damaged,
    or was written by an incompatible pymatgen version).
    """
    try:
        with open(path) as f:
            entries = json.load(f, cls=MontyDecoder)
    except Exception:
        return None
    if not isinstance(entries, list) or not all(isinstance(e, ComputedEntry) for e in entries):
        return None
    return entries


def get_entries_in_chemsys(chemsys, use_cache=True):
    """
    Get all MP entries in a chemical system, including all of its subsystems.

    :param chemsys: an iterable of element symbols or Elements, e.g. ["Li", "P", "S"]
    :param use_cache: whether to read from and write to the memory and disk caches
    :return: a new list of ComputedEntry objects. The entry objects themselves are
        shared with the cache, so do not modify them in place; copy them first.
    """
    elements = _normalize_chemsys(chemsys)
    key = (tuple(elements), _thermo_type)
    if use_cache:
        if key in _memory_cache:
            return list(_memory_cache[key])
        # A cached superset of this chemical system already holds every entry we need,
        # unless the energies depend on the chemical system that was queried.
        if _thermo_type not in CHEMSYS_DEPENDENT_THERMO_TYPES:
            for (cached_els, cached_type), cached in _memory_cache.items():
                if cached_type == _thermo_type and set(elements) <= set(cached_els):
                    return [e for e in cached if {el.symbol for el in e.composition.elements} <= set(elements)]

    path = _cache_path(elements, _thermo_type) if use_cache else None
    entries = None
    if path and os.path.isfile(path):
        entries = _read_cache(path)  # None if unusable; fetch again and overwrite it
    if entries is None:
        entries = _query("get_entries_in_chemsys", elements, additional_criteria={"thermo_types": [_thermo_type]})
        if path:
            _write_cache(path, entries)

    if use_cache:
        _memory_cache[key] = entries
    return list(entries)


def get_lowest_energy_entry(criteria):
    """
    Get the lowest-energy MP entry (per atom) among all polymorphs of a formula.
    The entry is taken from the chemical system query, so its energy is on the same
    scale as the other entries fetched by this module.

    :param criteria: a formula such as "Li3PS4", or an MP id such as "mp-985583"
    """
    from pymatgen.core import Composition

    if "-" in criteria:  # an MP id
        entries = _query("get_entries", criteria, additional_criteria={"thermo_types": [_thermo_type]})
        if not entries:
            raise ValueError("MP doesn't have any entry with id {}".format(criteria))
        return min(entries, key=lambda e: e.energy_per_atom)

    comp = Composition(criteria)
    entries = [e for e in get_entries_in_chemsys(comp.elements)
               if e.composition.reduced_composition == comp.reduced_composition]
    if not entries:
        raise ValueError("MP doesn't have any entry that matches the formula {}".format(criteria))
    return min(entries, key=lambda e: e.energy_per_atom)


def seed_cache(chemsys, entries):
    """
    Put entries in the memory cache for a chemical system, so that no MP query is made
    for it. Mainly useful for tests and offline work.
    """
    elements = _normalize_chemsys(chemsys)
    _memory_cache[(tuple(elements), _thermo_type)] = list(entries)


def clear_memory_cache():
    _memory_cache.clear()
