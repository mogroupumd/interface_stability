"""
Offline tests on the made-up data set in synthetic.py. Expected values are worked out by hand
from the energies listed there. No Materials Project access is needed.
"""
import os
import unittest

import matplotlib

matplotlib.use("Agg")

from pymatgen.core import Composition, Element  # noqa: E402
from pymatgen.analysis.phase_diagram import PhaseDiagram  # noqa: E402

from interface_stability import mpdata  # noqa: E402
from interface_stability.pseudobinary import PseudoBinary, get_mix_entry  # noqa: E402
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

    def test_bad_thermo_type(self):
        with self.assertRaises(ValueError):
            mpdata.set_thermo_type("PBE")


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


class TestPseudoBinary(OfflineTestCase):
    def make_pb(self, f1, f2):
        e1 = VirtualEntry.from_composition(f1)
        e2 = VirtualEntry.from_composition(f2)
        mix = VirtualEntry.from_composition(Composition(f1) + Composition(f2))
        entries = [e for e in mix.get_PD_entries() if e is not mix]
        e1.stabilize(entries=entries)
        e2.stabilize(entries=entries)
        return PseudoBinary(e1, e2, entries=entries), entries

    def test_does_not_modify_given_entries(self):
        pb, entries = self.make_pb("Li2S", "P2S5")
        self.assertEqual(len(pb.PDEntries), len(entries) + 2)

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

    def test_pseudo_binary_commands(self):
        from interface_stability.scripts.pseudo_binary import main

        self.assertIn("-31.25", self.run_script(main, ["pd", "Li2S", "P2S5"]))
        self.assertIn("Li3PS4", self.run_script(main, ["gppd", "Li2S", "P2S5", "Li", "-1.97"]))
        self.assertIn("Reaction Energy", self.run_script(main, ["gppd_screen", "Li2S", "P2S5", "Li", "-4", "0"]))
        self.assertIn("Li2S", self.run_script(main, ["gppd", "P2S5", "S", "Li", "-1"]))


if __name__ == "__main__":
    unittest.main()
