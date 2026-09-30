"""
Offline tests on the made-up data set in synthetic.py. Expected values are worked out by hand
from the energies listed there. No Materials Project access is needed.
"""
import os
import unittest
from unittest import mock

import matplotlib
import numpy

matplotlib.use("Agg")

from pymatgen.core import Composition, Element  # noqa: E402
from pymatgen.analysis.phase_diagram import PhaseDiagram  # noqa: E402
from pymatgen.entries.computed_entries import ComputedEntry  # noqa: E402

from interface_stability import mpdata  # noqa: E402
from interface_stability.pseudobinary import PseudoBinary, get_mix_entry, get_mutual_rxn_energies  # noqa: E402
from interface_stability.singlephase import VirtualEntry  # noqa: E402
from interface_stability.tests import synthetic  # noqa: E402

ALL_ELEMENTS = ["Li", "P", "S", "Co", "O"]


class OfflineTestCase(unittest.TestCase):
    def setUp(self):
        mpdata.clear_memory_cache()
        mpdata.seed_cache(ALL_ELEMENTS, synthetic.get_entries())

    def tearDown(self):
        mpdata.clear_memory_cache()

    def stabilized(self, formula, sup_el=None):
        entry = VirtualEntry.from_composition(formula)
        entries = entry.get_PD_entries(sup_el=sup_el)
        entry.stabilize(entries=entries)
        return entry, entries


class TestMPData(OfflineTestCase):
    def test_subsystem_is_served_from_seeded_superset(self):
        entries = mpdata.get_entries_in_chemsys(["Li", "S"])
        self.assertEqual(sorted(e.name for e in entries), ["Li", "Li2S", "S"])

    def test_lowest_energy_entry(self):
        self.assertEqual(mpdata.get_lowest_energy_entry("Li6P2S8").entry_id, "syn-Li3PS4")
        with self.assertRaisesRegex(ValueError, "LiS4"):
            mpdata.get_lowest_energy_entry("LiS4")

    def test_returned_list_is_a_copy(self):
        mpdata.get_entries_in_chemsys(["Li", "S"]).append("junk")
        self.assertNotIn("junk", mpdata.get_entries_in_chemsys(["Li", "S"]))

    def test_disk_cache_round_trip(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            os.environ["IFS_CACHE_DIR"] = tmp
            try:
                path = mpdata._cache_path(["Li", "S"], mpdata.get_thermo_type())
                import json
                from monty.json import MontyEncoder

                with open(path, "w") as f:
                    json.dump(mpdata.get_entries_in_chemsys(["Li", "S"]), f, cls=MontyEncoder)
                mpdata.clear_memory_cache()
                entries = mpdata.get_entries_in_chemsys(["Li", "S"])
                self.assertEqual(sorted(e.name for e in entries), ["Li", "Li2S", "S"])
            finally:
                del os.environ["IFS_CACHE_DIR"]

    def test_disk_cache_write_with_element_keys(self):
        # MP entries carry data["oxidation_states"] keyed by Element, which plain JSON cannot encode.
        import tempfile

        entries = mpdata.get_entries_in_chemsys(["Li", "S"])
        for e in entries:
            e.data = {"oxidation_states": {el: 0.0 for el in e.composition.elements}}
        with tempfile.TemporaryDirectory() as tmp:
            os.environ["IFS_CACHE_DIR"] = tmp
            try:
                path = mpdata._cache_path(["Li", "S"], mpdata.get_thermo_type())
                mpdata._write_cache(path, entries)
                self.assertIsInstance(next(iter(entries[0].data["oxidation_states"])), Element)
                mpdata.clear_memory_cache()
                cached = mpdata.get_entries_in_chemsys(["Li", "S"])
                self.assertEqual(sorted(e.name for e in cached), ["Li", "Li2S", "S"])
                self.assertEqual(cached[0].data["oxidation_states"], {cached[0].composition.elements[0].symbol: 0.0})
            finally:
                del os.environ["IFS_CACHE_DIR"]

    def test_bad_thermo_type(self):
        with self.assertRaises(ValueError):
            mpdata.set_thermo_type("PBE")

    def test_unusable_disk_cache_is_fetched_again(self):
        import tempfile

        fetched = synthetic.get_entries()
        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch.object(mpdata, "get_cache_dir", return_value=tmp), \
                mock.patch.object(mpdata, "_query", return_value=fetched) as query:
            path = mpdata._cache_path(["Li", "S"], mpdata.get_thermo_type())
            # A damaged file, and one that decodes to something other than entries
            for content in ("[{", '[{"@module": "no.such.module", "@class": "Entry"}]'):
                with open(path, "w") as f:
                    f.write(content)
                mpdata.clear_memory_cache()
                self.assertEqual(mpdata.get_entries_in_chemsys(["Li", "S"]), fetched)
                # The file is rewritten, and no temporary file is left behind
                self.assertEqual(os.listdir(tmp), [os.path.basename(path)])
                self.assertEqual(len(mpdata._read_cache(path)), len(fetched))
            self.assertEqual(query.call_count, 2)

    def test_mixed_thermo_type_does_not_use_superset(self):
        # Mixed GGA/GGA+U/r2SCAN energies depend on the chemical system queried
        subset = [e for e in synthetic.get_entries() if e.composition.chemical_system in ("Li", "S", "Li-S")]
        mpdata.set_thermo_type("GGA_GGA+U_R2SCAN")
        try:
            mpdata.seed_cache(ALL_ELEMENTS, synthetic.get_entries())
            with mock.patch.object(mpdata, "get_cache_dir", return_value=None), \
                    mock.patch.object(mpdata, "_query", return_value=subset) as query:
                self.assertEqual(mpdata.get_entries_in_chemsys(["Li", "S"]), subset)
            query.assert_called_once()
        finally:
            mpdata.set_thermo_type(mpdata.DEFAULT_THERMO_TYPE)

    def test_mixed_thermo_type_needs_recent_mp_api(self):
        with mock.patch("importlib.metadata.version", return_value="0.45.15"):
            with self.assertRaisesRegex(RuntimeError, mpdata.MIN_MP_API_VERSION_FOR_MIXED):
                mpdata.set_thermo_type("GGA_GGA+U_R2SCAN")
        self.assertEqual(mpdata.get_thermo_type(), mpdata.DEFAULT_THERMO_TYPE)

    def test_mp_errors_are_reported_as_mpdataerror(self):
        from mp_api.client.core import MPRestError

        class FailingRester:
            def __init__(self, api_key):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def get_entries_in_chemsys(self, *args, **kwargs):
                raise MPRestError("REST query returned with error status code 401")

        mpdata.clear_memory_cache()
        with mock.patch.object(mpdata, "_mp_client", return_value=(FailingRester, MPRestError)), \
                mock.patch.object(mpdata, "get_api_key", return_value="key"), \
                mock.patch.object(mpdata, "get_cache_dir", return_value=None):
            with self.assertRaisesRegex(mpdata.MPDataError, "401"):
                mpdata.get_entries_in_chemsys(["Li", "S"])


class TestSinglePhase(OfflineTestCase):
    def test_phase_equilibria(self):
        entry = VirtualEntry.from_composition("Li7P3S11")
        decomp, _ = entry.get_decomp_entries_and_e_above_hull()
        self.assertEqual(sorted(e.name for e in decomp), ["Li3PS4", "P2S5"])
        text = entry.get_printable_PE_data_in_pd()
        self.assertIn("Li7P3S11 -> 2.333 Li3PS4 + 0.3333 P2S5", text)

    def test_stabilize_puts_entry_on_hull(self):
        entry, entries = self.stabilized("Li3PS4")
        pd = PhaseDiagram([e for e in entries if e is not entry])
        _, e_hull = pd.get_decomp_and_e_above_hull(entry, allow_negative=True)
        self.assertAlmostEqual(e_hull, 0, places=7)
        self.assertLess(e_hull, 0)
        # stabilize does not touch the shared MP entries
        mp_li3ps4 = [e for e in mpdata.get_entries_in_chemsys(["Li", "P", "S"]) if e.name == "Li3PS4"][0]
        self.assertEqual(mp_li3ps4.correction, 0)

    def test_energy_correction_is_per_atom(self):
        entry = VirtualEntry.from_composition("Li3PS4", energy=-10)
        entry.energy_correction(0.1)
        self.assertAlmostEqual(entry.energy, -10 + 0.8)

    def test_decomposition_in_gppd(self):
        # Li3PS4 -> 3 Li + 0.5 P2S5 + 1.5 S.
        # dE = 3.5 * -0.3 - 8 * -0.95 = 6.55 eV at mu_Li = 0, so 6.55 + 3 * (-5) = -8.45 eV at mu_Li = -5
        entry = VirtualEntry.from_composition("Li3PS4")
        chempot = {"Li": -5}
        entries = entry.get_gppd_entries(chempot)
        entry.stabilize(entries=entries)
        decomp, rxn = entry.get_decomposition_in_gppd(chempot, entries=entries)
        self.assertEqual(sorted(e.name for e in decomp), ["P2S5", "S"])
        self.assertAlmostEqual(rxn.calculated_reaction_energy, -8.45, places=5)
        text = entry.get_printable_PE_and_decomposition_in_gppd(chempot, entries=entries)
        self.assertIn("Reaction energy: -8.45 eV per Li3PS4", text)

    def test_gppd_does_not_modify_entries(self):
        entry = VirtualEntry.from_composition("Li3PS4")
        entries = entry.get_gppd_entries({"Li": -5})
        entry.stabilize(entries=entries)
        energies = [e.energy for e in entries]
        entry.get_decomposition_in_gppd({"Li": -5}, entries=entries)
        entry.get_phase_evolution_profile("Li", entries=entries, allowpmu=True)
        self.assertEqual(energies, [e.energy for e in entries])

    def test_stability_window(self):
        # Reduction: Li3PS4 + 5 Li -> P + 4 Li2S, dE = -16.8 + 7.6 = -9.2 eV, so mu = -9.2 / 5 = -1.84
        # Oxidation: mu = -6.55 / 3 = -2.1833
        entry, entries = self.stabilized("Li3PS4", sup_el=["Li"])
        hi, lo = entry.get_stability_window("Li", entries=entries)
        self.assertAlmostEqual(hi, -1.84, places=5)
        self.assertAlmostEqual(lo, -6.55 / 3, places=5)

    def test_stability_window_posmu(self):
        # Li3P is the most Li-rich phase: stable up to Li metal (mu = 0), or without an upper bound once
        # positive mu is allowed. 3 Li + P -> Li3P, dE = 4 * -0.9 = -3.6 eV, so the lower bound is mu = -1.2
        entry, entries = self.stabilized("Li3P", sup_el=["Li"])
        hi, lo = entry.get_stability_window("Li", entries=entries)
        self.assertAlmostEqual(hi, 0, places=6)
        self.assertAlmostEqual(lo, -1.2, places=6)
        hi, lo = entry.get_stability_window("Li", entries=entries, allowpmu=True)
        self.assertIsNone(hi)
        self.assertAlmostEqual(lo, -1.2, places=6)

    def test_stability_window_without_lower_bound(self):
        # 2 Li + S -> Li2S, dE = 3 * -1.4 = -4.2 eV, so S is stable below mu_Li = -2.1
        entry, entries = self.stabilized("S", sup_el=["Li"])
        hi, lo = entry.get_stability_window("Li", entries=entries)
        self.assertAlmostEqual(hi, -2.1, places=5)
        self.assertIsNone(lo)

    def test_evolution_profile(self):
        entry, entries = self.stabilized("Li3PS4", sup_el=["Li"])
        profile = entry.get_phase_evolution_profile("Li", entries=entries)
        self.assertEqual([round(s["evolution"], 6) for s in profile], [8, 5, 0, -3])
        self.assertEqual([sorted(e.name for e in s["entries"]) for s in profile],
                         [["Li2S", "Li3P"], ["Li2S", "P"], ["Li3PS4"], ["P2S5", "S"]])
        # 3 Li + P -> Li3P, dE = 4 * -0.9 = -3.6 eV, so the Li3P / P transition is at mu = -1.2
        ref = profile[0]["element_reference"].energy_per_atom
        self.assertAlmostEqual(ref, synthetic.ELEMENT_ENERGIES["Li"])
        self.assertAlmostEqual(profile[1]["chempot"] - ref, -1.2, places=6)

    def test_printable_evolution_profile(self):
        entry, entries = self.stabilized("Li3PS4", sup_el=["Li"])
        text = entry.get_printable_evolution_profile("Li", entries=entries, plot_rxn_e=False)
        self.assertIn("Li3PS4 + 8 Li -> Li3P + 4 Li2S", text)
        # Reaction energy at mu = 0: (-3.6 - 16.8 + 7.6) eV / 8 atoms = -1.60 eV/atom
        self.assertIn("-1.60", text)

    def test_printable_evolution_profile_posmu(self):
        # Used to fail with UnboundLocalError
        entry, entries = self.stabilized("Li3PS4", sup_el=["Li"])
        text = entry.get_printable_evolution_profile("Li", entries=entries, plot_rxn_e=False, allowpmu=True)
        self.assertIn("inf", text)

    def test_evolution_plot_is_saved(self):
        import tempfile

        entry, entries = self.stabilized("Li3PS4", sup_el=["Li"])
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "rxn.png")
            entry.get_printable_evolution_profile("Li", entries=entries, save_path=path)
            self.assertTrue(os.path.getsize(path) > 0)

    def test_voltage_profile(self):
        entry, entries = self.stabilized("Li3PS4", sup_el=["Li"])
        oe_list, v_list = entry.get_vc_plot_data("Li", entries=entries)
        self.assertEqual([round(x, 6) for x in oe_list], [-3, 0, 5, 8])
        self.assertEqual([round(v, 4) for v in v_list[:3]], [2.1833, 1.84, 1.2])
        self.assertIn("Voltage ref. to Li", entry.get_printable_vc_plot_data("Li", oe_list, v_list))
        entry.get_voltage_profile_plot("Li", oe_list, v_list, 1)
        entry.get_voltage_profile_plot("Mg", oe_list, v_list, 2)

    def test_unknown_working_ion_needs_valence(self):
        entry, entries = self.stabilized("Li3PS4", sup_el=["Li"])
        with self.assertRaises(ValueError):
            entry.get_vc_plot_data("S", entries=entries)

    def test_import_leaves_matplotlib_settings(self):
        defaults = matplotlib.rc_params()
        for key in ("font.size", "mathtext.default"):
            self.assertEqual(matplotlib.rcParams[key], defaults[key], key)


class TestPseudoBinary(OfflineTestCase):
    def make_pb(self, f1, f2, e1_correction=0.0):
        e1 = VirtualEntry.from_composition(f1)
        e2 = VirtualEntry.from_composition(f2)
        mix = VirtualEntry.from_composition(Composition(f1) + Composition(f2))
        entries = [e for e in mix.get_PD_entries() if e is not mix]
        e1.stabilize(entries=entries)
        e2.stabilize(entries=entries)
        e1.energy_correction(e1_correction)
        return PseudoBinary(e1, e2, entries=entries), entries

    def assert_profile_matches_sampling(self, pd, entry1, entry2, profile, n=2001):
        """
        Each point of the profile has the phases that the phase diagram gives there, and each phase
        region met when sampling the mixing line is bounded by a point of the profile.
        """
        def phases(x):
            return frozenset(e.name for e in pd.get_decomposition(get_mix_entry(entry1, entry2, x).composition))

        points = [frozenset(e.name for e in decomp) for _, (decomp, _) in profile]
        for (x, _), names in zip(profile, points):
            self.assertEqual(names, phases(x), x)
        for x in numpy.linspace(0, 1, n):
            region = phases(x)
            self.assertTrue(any(p <= region for p in points), (x, sorted(region)))

    def test_does_not_modify_given_entries(self):
        pb, entries = self.make_pb("Li2S", "P2S5")
        self.assertEqual(len(pb.PDEntries), len(entries) + 2)
        self.assertIn(pb, {pb})  # hashable

    def test_mix_entry(self):
        a = VirtualEntry.from_composition("Li", energy=-1)
        b = VirtualEntry.from_composition("S", energy=-3)
        mix = get_mix_entry(a, b, 0.25)
        self.assertAlmostEqual(mix.composition["Li"], 0.25)
        self.assertAlmostEqual(mix.energy, -2.5)

    def test_pd_mixing(self):
        # Li2S + P2S5 -> Li3PS4 at x(Li2S) = 4.5 / 8 atoms = 0.5625.
        # E = -0.95 - (4.5 * -1.4 + 3.5 * -0.3) / 8 = -0.03125 eV/atom
        pb, _ = self.make_pb("Li2S", "P2S5")
        profile = pb.pd_mixing()
        self.assertEqual(len(profile), 3)
        x, (decomp, e) = profile[1]
        self.assertAlmostEqual(x, 0.5625, places=6)
        self.assertEqual([d.name for d in decomp], ["Li3PS4"])
        self.assertAlmostEqual(e, 0.03125, places=6)
        text = pb.get_printable_pd_profile()
        self.assertIn("-31.25", text)
        self.assertIn("Minimum", text)

    def test_pd_mixing_no_reaction(self):
        pb, _ = self.make_pb("Li3PS4", "Co3O4")
        profile = pb.pd_mixing()
        self.assertEqual(len(profile), 2)

    def test_pd_mixing_many_transitions(self):
        # Li + P2S5 crosses several phase fields; every transition point must be on the hull
        # of the mixing line, and the recursion must place each x at its real composition.
        pb, _ = self.make_pb("Li", "P2S5")
        profile = pb.pd_mixing()
        self.assertGreater(len(profile), 3)
        for x, (decomp, e) in profile:
            mix = get_mix_entry(pb.entry1, pb.entry2, x)
            self.assertAlmostEqual(pb.PD.get_decomp_and_e_above_hull(mix)[1], e, places=6)
        self.assert_profile_matches_sampling(pb.PD, pb.entry1, pb.entry2, profile)

    def test_pd_mixing_five_elements(self):
        # A made-up 5-element system (found by comparing with sampling) where the search used to skip the
        # phase regions with Co + CoSO2 near x = 0.82, and gave the wrong phases at the point it reported there.
        compounds = {"LiCo4P4O6": -1.874, "Co3P4S6": -1.784, "Co3S3O6": -1.645, "Li6P4S5O4": -1.916,
                     "Li2P6S2O3": -1.851}
        entries = [ComputedEntry(el, 0.0) for el in ("Li", "P", "S", "O", "Co")]
        entries += [ComputedEntry(f, e * Composition(f).num_atoms) for f, e in compounds.items()]
        e1 = VirtualEntry.from_composition("LiCo3P4SO3")
        e2 = VirtualEntry.from_composition("CoS4O2")
        e1.stabilize(entries=entries)
        e2.stabilize(entries=entries)
        pb = PseudoBinary(e1, e2, entries=entries)
        profile = pb.pd_mixing()
        self.assert_profile_matches_sampling(pb.PD, pb.entry1, pb.entry2, profile)
        self.assertEqual([round(x, 3) for x, _ in profile], [0, 0.515, 0.818, 0.823, 0.929, 1])
        self.assertEqual(sorted(d.name for d in profile[2][1][0]),
                         ["Co3(P2S3)2", "CoSO2", "Li6P4S5O4", "LiCo4(P2O3)2"])

    def test_pd_mixing_keeps_end_above_hull(self):
        # Li7P3S11 + Li2S -> 3 Li3PS4 at x(Li7P3S11) = 21 / 24. The Li7P3S11 composition lies on the Li3PS4-P2S5
        # tie line at (7/3 * 8 * -0.95 + 1/3 * 7 * -0.3) / 21 = -0.877778 eV/atom, so the reaction energy there
        # is 0.875 * -0.877778 + 0.125 * -1.4 + 0.95 = -1/144 eV/atom. Putting Li7P3S11 0.02 eV/atom above the
        # hull lowers the reaction energy by 0.875 * 0.02 but leaves the mutual reaction energy unchanged,
        # and the Li7P3S11 end (which then decomposes into Li3PS4 + P2S5) must stay in the profile.
        for correction in (0.0, 0.02):
            pb, _ = self.make_pb("Li7P3S11", "Li2S", e1_correction=correction)
            profile = pb.pd_mixing()
            self.assertEqual([round(x, 6) for x, _ in profile], [0, 0.875, 1])
            x, (decomp, e) = profile[1]
            self.assertEqual([d.name for d in decomp], ["Li3PS4"])
            self.assertAlmostEqual(-e, -1 / 144 - 0.875 * correction, places=6)
            self.assertAlmostEqual(get_mutual_rxn_energies(profile)[1], -1 / 144, places=6)
            self.assertIn("-6.94", pb.get_printable_pd_profile())

    def test_gppd_mixing(self):
        # At mu_Li = -1.97 Li3PS4 is stable. The Li3PS4 point is at the same x as in the closed system,
        # and its energy per non-Li atom is -0.03125 / 0.625 = -0.05 eV.
        pb, _ = self.make_pb("Li2S", "P2S5")
        chempots = {"Li": -1.97}
        profile = pb.gppd_mixing(chempots)
        self.assertEqual(chempots, {"Li": -1.97})  # not modified
        names = [sorted(d.name for d in decomp) for _, (decomp, _) in profile]
        self.assertIn(["Li3PS4"], names)
        x, (decomp, e) = profile[names.index(["Li3PS4"])]
        self.assertAlmostEqual(x, 0.5625, places=6)
        self.assertAlmostEqual(e, 0.05, places=6)
        self.assertIn("Minimum", pb.get_printable_gppd_profile({"Li": -1.97}))

    def test_gppd_mutual_rxn_energy(self):
        # At mu_Li = -2.15 Li2S is oxidized on its own: Li2S -> S + 2 Li, 4.2 + 2 * -2.15 = -0.1 eV per S atom.
        # Li3PS4 is still stable, so at x(Li2S) = 0.5625 the mixture forms Li3PS4 without exchanging Li:
        # -0.03125 eV per atom, i.e. -0.05 eV per non-Li atom (0.625 per atom). Li2S on its own would give
        # 0.5625 / 3 * -0.1 = -0.01875 eV per atom of the mixture, so the mutual reaction energy is
        # (-0.03125 + 0.01875) / 0.625 = -0.02 eV per non-Li atom. It used to come out as +6.25 meV.
        pb, _ = self.make_pb("Li2S", "P2S5")
        n1, n2 = pb.get_non_open_fractions("Li")
        self.assertAlmostEqual(n1, 1 / 3)
        self.assertAlmostEqual(n2, 1)
        profile = pb.gppd_mixing({"Li": -2.15})
        self.assertEqual([sorted(d.name for d in decomp) for _, (decomp, _) in profile],
                         [["P2S5"], ["Li3PS4"], ["S"]])
        mutual = get_mutual_rxn_energies(profile, n1, n2)
        self.assertAlmostEqual(profile[1][1][1], 0.05, places=6)
        self.assertAlmostEqual(profile[-1][1][1], 0.1, places=6)
        self.assertAlmostEqual(mutual[1], -0.02, places=6)
        self.assertIn("-20.00", pb.get_printable_gppd_profile({"Li": -2.15}))

    def test_gppd_keeps_both_ends(self):
        # At mu_Li = -3, Li3PS4 -> P2S5 + S and Li2S -> S on their own, and a mixture of them does no more.
        # Both ends stay in the profile whichever phase comes first.
        for f1, f2 in (("Li3PS4", "Li2S"), ("Li2S", "Li3PS4")):
            pb, _ = self.make_pb(f1, f2)
            profile = pb.gppd_mixing({"Li": -3})
            self.assertEqual([x for x, _ in profile], [0, 1])
            ends = {f2: profile[0][1], f1: profile[1][1]}
            self.assertEqual(sorted(d.name for d in ends["Li3PS4"][0]), ["P2S5", "S"])
            self.assertEqual(sorted(d.name for d in ends["Li2S"][0]), ["S"])
            self.assertEqual(get_mutual_rxn_energies(profile, *pb.get_non_open_fractions("Li")), [0, 0])

    def test_gppd_mixing_matches_sampling(self):
        from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram, GrandPotPDEntry

        pb, _ = self.make_pb("Li2S", "P2S5")
        li_ref = PhaseDiagram(pb.PDEntries).el_refs[Element("Li")].energy_per_atom
        for mu in (-0.5, -1.5, -2.0, -2.15, -3.0):
            chempots = {Element("Li"): mu + li_ref}
            gppd = GrandPotentialPhaseDiagram(pb.PDEntries, chempots)
            g1, g2 = GrandPotPDEntry(pb.entry1, chempots), GrandPotPDEntry(pb.entry2, chempots)
            self.assert_profile_matches_sampling(gppd, g1, g2, pb.gppd_mixing({"Li": mu}), n=501)

    def test_gppd_open_element_only(self):
        pb, _ = self.make_pb("Li", "P2S5")
        with self.assertRaisesRegex(ValueError, "other than the open element"):
            pb.gppd_mixing({"Li": -1})

    def test_gppd_open_element_not_in_phases(self):
        # Used to fail: get_GPPD_entries did not exist
        pb, _ = self.make_pb("P2S5", "S")
        entries = pb.get_gppd_entries("Li")
        self.assertIn(Element("Li"), {el for e in entries for el in e.composition.elements})
        self.assertIn(pb.PDEntries[-1], entries)
        profile = pb.gppd_mixing({"Li": -3})
        self.assertGreaterEqual(len(profile), 2)

    def test_transition_chempots(self):
        pb, _ = self.make_pb("Li2S", "P2S5")
        mus = pb.get_gppd_transition_chempots("Li")
        self.assertEqual(list(mus), sorted(mus, reverse=True))
        for expected in (-1.2, -1.84, -2.1, -6.55 / 3):
            self.assertTrue(any(abs(m - expected) < 1e-6 for m in mus), expected)

    def test_gppd_scanning(self):
        pb, _ = self.make_pb("Li2S", "P2S5")
        text = pb.gppd_scanning("Li", 0, -4)
        lines = text.splitlines()
        self.assertIn("mu_high", lines[1])
        self.assertLess(lines[1].index("mu_high"), lines[1].index("mu_low"))
        li3ps4_row = [ln for ln in lines if ln.strip().endswith("Li3PS4")][0]
        self.assertEqual(li3ps4_row.split()[:2], ["-1.84", "-2.18"])


class TestScripts(OfflineTestCase):
    def run_script(self, main, argv):
        import contextlib
        import io
        import sys

        old = sys.argv
        sys.argv = ["prog"] + argv
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                main()
        finally:
            sys.argv = old
            mpdata.set_thermo_type(mpdata.DEFAULT_THERMO_TYPE)
        return out.getvalue()

    def test_phase_stability_commands(self):
        from interface_stability.scripts.phase_stability import main

        self.assertIn("2.333 Li3PS4 + 0.3333 P2S5", self.run_script(main, ["stability", "Li7P3S11"]))
        self.assertIn("-8.45", self.run_script(main, ["mu", "Li3PS4", "Li", "-5"]))
        self.assertIn("Li3PS4 + 5 Li -> P + 4 Li2S",
                      self.run_script(main, ["evolution", "Li3PS4", "Li", "--noplot"]))
        self.assertIn("inf", self.run_script(main, ["evolution", "-posmu", "Li3PS4", "Li", "--noplot"]))
        self.assertIn("1.84", self.run_script(main, ["plotvc", "Li3PS4", "Li", "--noplot"]))

    def test_bad_input_is_reported_without_traceback(self):
        from interface_stability.scripts import phase_stability, pseudo_binary

        with self.assertRaisesRegex(SystemExit, "Error: .*Xx"):
            self.run_script(phase_stability.main, ["stability", "Li3Xx"])
        with self.assertRaisesRegex(SystemExit, "Error: .*other than the open element"):
            self.run_script(pseudo_binary.main, ["gppd", "Li", "Li3PS4", "Li", "-1"])

    def test_pseudo_binary_commands(self):
        from interface_stability.scripts.pseudo_binary import main

        self.assertIn("-31.25", self.run_script(main, ["pd", "Li2S", "P2S5"]))
        self.assertIn("Li3PS4", self.run_script(main, ["gppd", "Li2S", "P2S5", "Li", "-1.97"]))
        self.assertIn("Reaction Energy", self.run_script(main, ["gppd_screen", "Li2S", "P2S5", "Li", "-4", "0"]))
        self.assertIn("Li2S", self.run_script(main, ["gppd", "P2S5", "S", "Li", "-1"]))


if __name__ == "__main__":
    unittest.main()
