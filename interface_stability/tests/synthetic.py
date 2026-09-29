"""
A small made-up Li-P-S (+ Co, O) data set for offline tests. The energies are not real,
but chosen so the expected answers can be worked out by hand.

Element energies (eV/atom): Li -1.9, P -5.4, S -4.1, Co -7.0, O -4.9
Formation energies (eV/atom): Li2S -1.4, Li3P -0.9, P2S5 -0.3, Li3PS4 -0.95, Co3O4 -1.3

Li3PS4 lies below the Li2S-P2S5 tie line, whose energy at that composition is
(4.5 * -1.4 + 3.5 * -0.3) / 8 = -0.91875 eV/atom.
"""

from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry

ELEMENT_ENERGIES = {"Li": -1.9, "P": -5.4, "S": -4.1, "Co": -7.0, "O": -4.9}

FORMATION_ENERGIES = {
    "Li2S": -1.4,
    "Li3P": -0.9,
    "P2S5": -0.3,
    "Li3PS4": -0.95,
    "Co3O4": -1.3,
}


def make_entry(formula, formation_energy_per_atom=0.0, entry_id=None):
    comp = Composition(formula)
    energy = sum(amt * ELEMENT_ENERGIES[el.symbol] for el, amt in comp.items())
    energy += formation_energy_per_atom * comp.num_atoms
    return ComputedEntry(comp, energy, entry_id=entry_id)


def get_entries():
    entries = [make_entry(el, 0.0, entry_id="syn-el-" + el) for el in ELEMENT_ENERGIES]
    entries += [make_entry(f, ef, entry_id="syn-" + f) for f, ef in FORMATION_ENERGIES.items()]
    return entries
