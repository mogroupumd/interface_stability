import unittest
from pymatgen import MPRester, SETTINGS

#@unittest.skipIf(not SETTINGS.get("PMG_MAPI_KEY"), "PMG_MAPI_KEY environment variable not set.")

class InterfaceStabilityTest(unittest.TestCase):
    def test_get_api_key(self):
        self.assertTrue(SETTINGS.get("PMG_MAPI_KEY"))


    def test_get_structure_by_material_id(self):
        with MPRester() as m:
            s1 = m.get_structure_by_material_id("mp-1")
        self.assertEqual(s1.formula, "Cs1")

if __name__ == "__main__":
    unittest.main()