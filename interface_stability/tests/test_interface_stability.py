"""
Live tests against the Materials Project. They are skipped unless an API key is set
(MP_API_KEY environment variable, or PMG_MAPI_KEY in ~/.pmgrc.yaml).

They check the README examples qualitatively: the MP data has changed since the README
was written, so exact numbers are not compared.
"""
import unittest

import matplotlib

matplotlib.use("Agg")

from interface_stability import mpdata  # noqa: E402
from interface_stability.pseudobinary import PseudoBinary  # noqa: E402
from interface_stability.singlephase import VirtualEntry  # noqa: E402


@unittest.skipIf(not mpdata.get_api_key(), "No Materials Project API key set")
class InterfaceStabilityLiveTest(unittest.TestCase):
    def test_get_entries(self):
        entries = mpdata.get_entries_in_chemsys(["Li", "S"])
        self.assertIn("Li2S", {e.composition.reduced_formula for e in entries})

    def test_lgps_phase_equilibria(self):
        entry = VirtualEntry.from_composition("Li10GeP2S12")
        decomp, _ = entry.get_decomp_entries_and_e_above_hull()
        self.assertEqual(sorted(e.name for e in decomp), ["Li3PS4", "Li4GeS4"])

    def test_li3ps4_stability_window(self):
        entry = VirtualEntry.from_composition("Li3PS4")
        entries = entry.get_PD_entries(sup_el=["Li"])
        entry.stabilize(entries=entries)
        hi, lo = entry.get_stability_window("Li", entries=entries)
        # README: stable from about 1.7 V to 2.3 V vs. Li
        self.assertTrue(-2.0 < hi < -1.4, hi)
        self.assertTrue(-2.7 < lo < -2.0, lo)

    def test_licoo2_li3ps4_reacts(self):
        e1 = VirtualEntry.from_composition("LiCoO2")
        e2 = VirtualEntry.from_composition("Li3PS4")
        mix = VirtualEntry.from_composition("LiCoO2Li3PS4")
        entries = [e for e in mix.get_PD_entries() if e is not mix]
        e1.stabilize(entries=entries)
        e2.stabilize(entries=entries)
        pb = PseudoBinary(e1, e2, entries=entries)
        profile = pb.pd_mixing()
        min_e = max(e for _, (_, e) in profile)
        # README: about -400 meV/atom
        self.assertGreater(min_e, 0.2)


if __name__ == "__main__":
    unittest.main()
