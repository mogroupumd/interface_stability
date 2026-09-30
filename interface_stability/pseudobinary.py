# coding: utf-8
# Copyright (c) Mogroup  @ University of Maryland, College Park
# Distributed under the terms of the MIT License.

import pandas
from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, GrandPotPDEntry
from pymatgen.analysis.reaction_calculator import ComputedReaction, ReactionError
from interface_stability.singlephase import VirtualEntry, element_symbol


__author__ = "Yizhou Zhu"
__copyright__ = ""
__version__ = "2.2"
__maintainer__ = "Yizhou Zhu"
__email__ = "yizhou.zhu@gmail.com"
__status__ = "Production"
__date__ = "Jun 10, 2018"

# Smallest mixing ratio interval the binary search splits further
MIN_RATIO_INTERVAL = 1e-9


class PseudoBinary(object):
    """
    A class for performing analyses on pseudo-binary stability calculations.

    The algorithm is based on the work in the following paper:

    Yizhou Zhu, Xingfeng He, Yifei Mo*, “First Principles Study on Electrochemical and Chemical Stability of the
    Solid Electrolyte-Electrode Interfaces in All-Solid-State Li-ion Batteries”, Journal of Materials Chemistry A, 4,
    3253-3266 (2016)
    DOI: 10.1039/c5ta08574h
    """

    def __init__(self, entry1, entry2, entries=None, sup_el=None):
        comp1 = entry1.composition
        comp2 = entry2.composition
        norm1 = 1.0 / entry1.composition.num_atoms
        norm2 = 1.0 / entry2.composition.num_atoms

        self.entry1 = VirtualEntry.from_composition(entry1.composition * norm1, energy=entry1.energy * norm1,
                                                    name=comp1.reduced_formula)
        self.entry2 = VirtualEntry.from_composition(entry2.composition * norm2, energy=entry2.energy * norm2,
                                                    name=comp2.reduced_formula)

        if not entries:
            entry_mix = VirtualEntry.from_composition(comp1 + comp2)
            entries = [e for e in entry_mix.get_PD_entries(sup_el=sup_el) if e is not entry_mix]
        entries = list(entries) + [entry1, entry2]
        self.PDEntries = entries
        self.PD = PhaseDiagram(entries)

    def pd_mixing(self):
        """
        This function give the phase equilibria of a pseudo-binary in a closed system (PD).
        It will give a complete evolution profile for mixing ratio x change from 0 to 1.
        x is the ratio (both entry norm. to 1 atom/fu) or each entry
        """
        profile = get_full_evolution_profile(self.PD, self.entry1, self.entry2, 0.0, 1.0)
        cleaned = clean_profile(profile)
        return cleaned

    def get_printable_pd_profile(self):
        return self.get_printed_profile(self.pd_mixing())

    def get_printable_gppd_profile(self, chempots, gppd_entries=None):
        profile = self.gppd_mixing(chempots, gppd_entries=gppd_entries)
        n1, n2 = self.get_non_open_fractions(list(chempots.keys())[0])
        return self.get_printed_profile(profile, n1, n2)

    def get_non_open_fractions(self, open_el):
        """
        The fraction of the atoms of entry1 and of entry2 that are not the open element.
        Energies in a grand potential phase diagram are per atom of these elements.
        """
        open_el = Element(element_symbol(open_el))
        return tuple(1 - entry.composition.get_atomic_fraction(open_el) for entry in (self.entry1, self.entry2))

    def get_printed_profile(self, profile, n1=1.0, n2=1.0):
        """
        A general function to generate printable table strings for pseudo-binary mixing results
        :param n1, n2: for a profile from gppd_mixing, the fraction of the atoms of entry1 and entry2
            that are not the open element (see get_mutual_rxn_energies)
        """
        output = ['\n ===  Pseudo-binary evolution profile  === ']
        df = pandas.DataFrame()
        mutual_rxn_e = get_mutual_rxn_energies(profile, n1, n2)

        x1s, x2s, es, mes, pes = [], [], [], [], []

        for (ratio, (decomp, e)), mutual_e in zip(profile, mutual_rxn_e):
            x1s.append(1-ratio)
            x2s.append(ratio)
            es.append(-e*1000)
            mes.append(mutual_e * 1000)
            pes.append(", ".join([x.name for x in decomp]))
        df["x({})".format(self.entry2.name)] = x1s
        df["x({})".format(self.entry1.name)] = x2s
        df["Rxn. E. (meV/atom)"] = es
        df["Mutual Rxn. E. (meV/atom)"] = mes
        df["Phase Equilibria"] = pes

        comments = ["" for _ in range(len(profile))]
        min_loc = list(df[df.columns[2:4]].idxmin())
        if min_loc[0] == min_loc[1]:
            comments[min_loc[0]] = 'Minimum'
        else:
            comments[min_loc[0]] = 'Rxn. E. Min.'
            comments[min_loc[1]] = 'Mutual Rxn. E. Min.'
        df["Comment"] = comments

        print_df = df.to_string(index=False, float_format='{:,.2f}'.format, justify='center')
        output.append(print_df)
        string = '\n'.join(output)
        return string

    def gppd_mixing(self, chempots, gppd_entries=None):
        """
        This function give the phase equilibria of a pseudo-binary in a open system (GPPD).
        It will give a complete evolution profile for mixing ratio x change from 0 to 1.
        x is the ratio of each entry, both entries norm. to 1 atom/fu (open element included, as in pd_mixing).
        The energies are per atom of the mixture other than the open element.
        """
        if len(chempots) != 1:
            raise ValueError("gppd_mixing supports exactly one open element")
        open_el, mu = list(chempots.items())[0]
        open_el = Element(element_symbol(open_el))
        if 0 in self.get_non_open_fractions(open_el):
            raise ValueError("Both phases must contain an element other than the open element {}".format(open_el))
        if not gppd_entries:
            gppd_entries = self.get_gppd_entries(open_el)
        # Use the element reference from the same entries, so all energies share one scale
        el_ref = PhaseDiagram(gppd_entries).el_refs[open_el]
        chempots = {open_el: mu + el_ref.energy_per_atom}
        gppd_entry1 = GrandPotPDEntry(self.entry1, chempots)
        gppd_entry2 = GrandPotPDEntry(self.entry2, chempots)
        gppd = GrandPotentialPhaseDiagram(gppd_entries, chempots)
        profile = get_full_evolution_profile(gppd, gppd_entry1, gppd_entry2, 0.0, 1.0)
        cleaned = clean_profile(profile)
        return cleaned

    def get_gppd_entries(self, open_el):
        open_el = Element(element_symbol(open_el))
        if open_el in (self.entry1.composition + self.entry2.composition):
            gppd_entries = self.PDEntries
        else:
            comp = self.entry1.composition + self.entry2.composition + Composition(open_el.symbol)
            entry_mix = VirtualEntry.from_composition(comp)
            gppd_entries = [e for e in entry_mix.get_PD_entries() if e is not entry_mix]
            # The two phases themselves (with their stabilized / corrected energies)
            gppd_entries += [e for e in self.PDEntries[-2:]]
        return gppd_entries

    def get_gppd_transition_chempots(self, open_el, gppd_entries=None):
        """
        This is to get all possible transition chemical potentials from PD (rather than GPPD)
        Still use pure element ref.
        # May consider supporting negative miu in the future
        """
        open_el = Element(element_symbol(open_el))
        if not gppd_entries:
            gppd_entries = self.get_gppd_entries(open_el)
        pd = PhaseDiagram(gppd_entries)
        vaspref_mius = pd.get_transition_chempots(open_el)
        el_ref = pd.el_refs[open_el]

        elref_mius = [miu - el_ref.energy_per_atom for miu in vaspref_mius]
        return elref_mius

    def gppd_scanning(self, open_el, mu_hi, mu_lo, gppd_entries=None, verbose=False):
        """
        This function is to do a (slightly smarter) screening of GPPD pseudo-binary in a given miu range
        This is a very tedious function, but mainly because GPPD screening itself is very tedious.

        :param open_el: open element
        :param mu_hi:  chemical potential upper bound
        :param mu_lo:  chemical potential lower bound
        :param gppd_entries: Supply GPPD entries manually. If you supply this, I assume you know what you are doing
        :param verbose: if True, list every chemical potential interval in the PE result table; otherwise merge
            neighboring intervals with the same phase equilibria
        :return: a printable string of screening results
        """
        mu_lo, mu_hi = sorted([mu_lo, mu_hi])
        if not gppd_entries:
            gppd_entries = self.get_gppd_entries(open_el)
        n1, n2 = self.get_non_open_fractions(open_el)

        def min_mutual_point(miu):
            """
            (phase equilibria, mutual reaction energy, reaction energy) at the mixing ratio
            with the lowest mutual reaction energy
            """
            profile = self.gppd_mixing({open_el: miu}, gppd_entries)
            mutual = get_mutual_rxn_energies(profile, n1, n2)
            i_min = min(range(len(profile)), key=lambda i: mutual[i])
            ratio, (decomp, e) = profile[i_min]
            return decomp, mutual[i_min], -e

        # Transition chemical potentials strictly inside the range, from high to low
        miu_E_candidates = [miu for miu in self.get_gppd_transition_chempots(open_el, gppd_entries) if
                            mu_lo < miu < mu_hi]
        miu_E_candidates = [mu_hi] + miu_E_candidates + [mu_lo]
        duplicate_index = []

        for i in range(1, len(miu_E_candidates) - 1):
            miu_left = (miu_E_candidates[i] + miu_E_candidates[i - 1]) / 2.0
            miu_right = (miu_E_candidates[i] + miu_E_candidates[i + 1]) / 2.0
            profile_left = self.gppd_mixing({open_el: miu_left}, gppd_entries)
            profile_right = self.gppd_mixing({open_el: miu_right}, gppd_entries)
            if judge_same_decomp(profile_left, profile_right):
                duplicate_index.append(i)
        miu_E_candidates = [miu_E_candidates[i] for i in range(len(miu_E_candidates)) if i not in duplicate_index]

        interval_mu_hi, PE = [], []
        mu_list, E_mutual_list, E_total_list = [], [], []

        for i in range(1, len(miu_E_candidates)):
            miu = (miu_E_candidates[i] + miu_E_candidates[i - 1]) / 2.0
            decomp = min_mutual_point(miu)[0]
            interval_mu_hi.append(miu_E_candidates[i - 1])
            PE.append(", ".join(sorted([x.name for x in decomp])))
        for miu in miu_E_candidates:
            _, e_mutual, e_total = min_mutual_point(miu)
            mu_list.append(miu)
            E_mutual_list.append(e_mutual)
            E_total_list.append(e_total)

        to_be_hidden = []
        if not verbose:
            for i in range(1, len(PE)):
                if PE[i] == PE[i - 1]:
                    to_be_hidden.append(i)

        mu_hi_display_list = [interval_mu_hi[k] for k in range(len(interval_mu_hi)) if k not in to_be_hidden]
        mu_low_display_list = mu_hi_display_list[1:] + [mu_lo]
        PE_display_list = [PE[k] for k in range(len(PE)) if k not in to_be_hidden]

        df1 = pandas.DataFrame()
        df2 = pandas.DataFrame()

        df1['mu_high'] = mu_hi_display_list
        df1['mu_low'] = mu_low_display_list
        df1['phase equilibria'] = PE_display_list

        df2['mu'] = mu_list
        df2['E_mutual(eV/atom)'] = E_mutual_list
        df2['E_total(eV/atom)'] = E_total_list

        print_df1 = df1.to_string(index=False, float_format='{:,.2f}'.format, justify='center')
        print_df2 = df2.to_string(index=False, float_format='{:,.2f}'.format, justify='center')

        output = [' == Phase Equilibria at min E_mutual == ', print_df1,'\n', ' == Reaction Energy ==',
                  print_df2, 'Note: if E_mutual = 0, E_total is at x = 1 or 0']
        string = "\n".join(output)
        return string


"""
The following functions are auxiliary functions.
Most of them are used to solve or clean the mixing PE profile.
"""


def judge_same_decomp(profile1, profile2):
    """
    Judge whether two profiles have identical decomposition products
    """
    if len(profile1) != len(profile2):
        return False
    for step in range(len(profile1)):
        ratio1, (decomp1, e1) = profile1[step]
        ratio2, (decomp2, e1) = profile2[step]
        if abs(ratio1 - ratio2) > 1e-8:
            return False
        names1 = sorted([x.name for x in decomp1])
        names2 = sorted([x.name for x in decomp2])
        if names1 != names2:
            return False
    return True


def get_full_evolution_profile(pd, entry1, entry2, x1, x2):
    """
    This function is used to solve the transition points along a path on convex hull.
    The essence is to use binary search, which is more accurate and faster than brutal force screening
    This is a recursive function.
    :param pd: PhaseDiagram of GrandPotentialPhaseDiagram
    :param entry1 & entry2: mixing entry1/entry2, PDEntry for pd_mixing, GrandPotEntry for gppd_mixing
    :param x1 & x2: The mixing ratio range for binary search.
    :return: An uncleaned but complete profile with all transition points.
    """
    evolution_profile = {}
    entry_left = get_mix_entry(entry1, entry2, x1)
    entry_right = get_mix_entry(entry1, entry2, x2)
    (decomp1, h1) = pd.get_decomp_and_e_above_hull(entry_left)
    (decomp2, h2) = pd.get_decomp_and_e_above_hull(entry_right)
    decomp1 = set(decomp1.keys())
    decomp2 = set(decomp2.keys())
    evolution_profile[x1] = (decomp1, h1)
    evolution_profile[x2] = (decomp2, h2)

    # If the phases at one end are a subset of those at the other end, the whole range lies in the
    # phase region of the other end, so there is no transition inside it.
    if decomp1 <= decomp2 or decomp2 <= decomp1:
        return evolution_profile

    intersect = decomp1 & decomp2
    if len(intersect) > 0:
        # This is try to catch a single transition point: where the path crosses the boundary made of the
        # phases shared by both ends. It is only accepted if the mixture there really decomposes into
        # shared phases, since the path can also leave through other phase regions.
        x = _get_transition_ratio(entry_left, entry_right, intersect, x1, x2)
        if x is not None:
            (decomp_x, h_x) = pd.get_decomp_and_e_above_hull(get_mix_entry(entry1, entry2, x))
            decomp_x = set(decomp_x.keys())
            if decomp_x <= intersect:
                evolution_profile[x] = (decomp_x, h_x)
                return evolution_profile

    if x2 - x1 < MIN_RATIO_INTERVAL:
        return evolution_profile

    x_mid = (x1 + x2) / 2.0
    entry_mid = get_mix_entry(entry1, entry2, x_mid)
    (decomp_mid, h_mid) = pd.get_decomp_and_e_above_hull(entry_mid)
    decomp_mid = set(decomp_mid.keys())
    evolution_profile[x_mid] = (decomp_mid, h_mid)
    part1 = get_full_evolution_profile(pd, entry1, entry2, x1, x_mid)
    part2 = get_full_evolution_profile(pd, entry1, entry2, x_mid, x2)
    evolution_profile.update(part1)
    evolution_profile.update(part2)
    return evolution_profile


def _get_transition_ratio(entry_left, entry_right, phases, x1, x2):
    """
    The mixing ratio strictly between x1 and x2 (the ratios of entry_left and entry_right) at which the
    mixture can be made of the given phases alone, or None if there is no such ratio.
    """
    try:
        rxn = ComputedReaction([entry_left, entry_right], list(phases))
    except ReactionError:
        return None
    c1 = _entry_coeff(rxn, entry_left)
    c2 = _entry_coeff(rxn, entry_right)
    if c1 is None or c2 is None or c1 * c2 <= 0:
        return None
    x = (c1 * x1 + c2 * x2) / (c1 + c2)
    return x if x1 < x < x2 else None


def _entry_coeff(rxn, entry):
    """
    The coefficient of an entry in a reaction, in units of the entry's own composition
    (not its reduced composition). None if the entry is not in the reaction.
    Entries are matched on their composition, which for GrandPotPDEntry excludes the open element.
    """
    comp = entry.composition
    reduced = comp.reduced_composition
    if reduced not in rxn.all_comp:
        return None
    return rxn.get_coeff(reduced) * reduced.num_atoms / comp.num_atoms


def clean_profile(evolution_profile):
    """
    This function is to clean the calculated profile from binary search. Redundant trial results are pruned out,
    with only the two ends and the transition points left.
    """
    raw_data = sorted(evolution_profile.items(), key=lambda item: item[0])
    clean_set = [raw_data[0]]
    for i in range(1, len(raw_data) - 1):
        x, (decomp, h) = raw_data[i]
        x_cpr, (decomp_cpr, h_cpr) = clean_set[-1]
        if set(decomp_cpr) <= set(decomp):
            continue
        else:
            clean_set.append(raw_data[i])
    # Always keep the far end, even when it is in the same phase region as the points before it
    if len(raw_data) > 1:
        clean_set.append(raw_data[-1])
    return clean_set


def get_mutual_rxn_energies(profile, n1=1.0, n2=1.0):
    """
    The mutual reaction energy at each point of a cleaned mixing profile: the reaction energy of the
    mixture minus that of entry1 and entry2 decomposing on their own, per atom of the mixture.

    The energies in the profile are per atom of the mixture that the phase diagram counts. In a grand
    potential phase diagram these are the atoms other than the open element, so the two ends have to be
    weighted by the counted atoms they bring to the mixture, not just by the mixing ratio x.
    :param profile: a profile from pd_mixing or gppd_mixing, which starts at x = 0 and ends at x = 1
    :param n1, n2: the counted atoms per atom of entry1 and entry2. 1 for pd_mixing; for gppd_mixing,
        the fraction of their atoms that are not the open element.
    :return: a list of mutual reaction energies (eV/atom), one for each point of the profile
    """
    x_first, (_, h2) = profile[0]
    x_last, (_, h1) = profile[-1]
    if x_first != 0 or x_last != 1:
        raise ValueError("The profile must start at x = 0 and end at x = 1")
    mutual = []
    for x, (_, h) in profile:
        n = x * n1 + (1 - x) * n2
        mutual.append(-h + (x * n1 * h1 + (1 - x) * n2 * h2) / n)
    return mutual


def get_mix_entry(entry1, entry2, x):
    """
    Mixing PDEntry or GrandPotEntry for the binary search algorithm.
    :return: the mixture x * entry1 + (1 - x) * entry2
    """
    x1, x2 = x, 1 - x
    if isinstance(entry1, GrandPotPDEntry):
        mid_ori_entry = VirtualEntry.from_mixing({entry1.original_entry: x1, entry2.original_entry: x2})
        return GrandPotPDEntry(mid_ori_entry, entry1.chempots)
    else:
        return VirtualEntry.from_mixing({entry1: x1, entry2: x2})
