# coding: utf-8
# Copyright (c) Yifei Mo Group @ University of Maryland, College Park
# Distributed under the terms of the MIT License.


import re

import pandas

import matplotlib.pyplot as plt
from pymatgen.core import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.entries.computed_entries import ComputedEntry, ConstantEnergyAdjustment

from interface_stability import mpdata

__author__ = "Yizhou Zhu"
__copyright__ = ""
__version__ = "3.0"
__maintainer__ = "Yizhou Zhu"
__email__ = "yizhou.zhu@gmail.com"
__status__ = "Production"
__date__ = "Jun 10, 2018"

plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['font.size'] = 15

# Working ions and their charges, used to convert chemical potential to voltage
COMMON_WORKING_IONS = {'Li': 1, 'Na': 1, 'K': 1, 'Mg': 2, 'Ca': 2, 'Zn': 2, 'Al': 3}


def add_energy_adjustment(entry, value, name="Constant energy adjustment"):
    """
    Add a constant energy (eV per formula unit of the entry) to an entry, in place.
    """
    entry.energy_adjustments.append(ConstantEnergyAdjustment(value, name=name))


def shifted_copy(entry, shift_per_atom, name="Chemical potential shift"):
    """
    Return a copy of an entry with its energy shifted by shift_per_atom (eV/atom).
    The original entry is not changed.
    """
    new = ComputedEntry(entry.composition, entry.uncorrected_energy,
                        energy_adjustments=list(entry.energy_adjustments),
                        parameters=entry.parameters, data=entry.data, entry_id=entry.entry_id)
    new.name = entry.name
    add_energy_adjustment(new, shift_per_atom * entry.composition.num_atoms, name=name)
    return new


def element_symbol(el):
    """
    Accept an Element or a symbol string and return the symbol string.
    """
    return el.symbol if isinstance(el, Element) else Element(el).symbol


class VirtualEntry(ComputedEntry):
    def __init__(self, composition, energy, name=None):
        super(VirtualEntry, self).__init__(Composition(composition), energy)
        if name:
            self.name = name

    @classmethod
    def from_composition(cls, comp, energy=0, name=None):
        return cls(Composition(comp), energy, name=name)

    @classmethod
    def from_mixing(cls, mixing_dict):
        comp = Composition()
        energy = 0
        for i in mixing_dict.keys():
            comp += Composition({el: i.composition[el] * mixing_dict[i] for el in i.composition.keys()})
            energy += mixing_dict[i] * i.energy
        return cls(Composition(comp), energy)

    @classmethod
    def from_mp(cls, criteria):
        entry = cls.get_mp_entry(criteria)
        return cls(entry.composition, energy=entry.energy, name=entry.name)

    @staticmethod
    def get_mp_entry(criteria):
        """
        Here always return the lowest energy among all polymorphs.
        Criteria can be a formula or an mp-id
        """
        return mpdata.get_lowest_energy_entry(criteria)

    @property
    def chemsys(self):
        return [_.symbol for _ in self.composition.elements]

    def get_PD_entries(self, sup_el=None, exclusions=None, trypreload=False):
        """
        :param sup_el: a list for extra element dimension, using str format
        :param exclusions: a list of manually exclusion entries, can use entry name or mp_id
        :param trypreload: Kept for backward compatibility and ignored. Entries are always cached in memory,
            and also on disk if IFS_CACHE_DIR or PMG_PD_PRELOAD_PATH is set. See interface_stability.mpdata.
        :return: all related entries to construct phase diagram, with this entry appended.
        """
        chemsys = self.chemsys + [element_symbol(el) for el in sup_el] if sup_el else self.chemsys
        chemsys = list(set(chemsys))

        entries = self.get_PD_entries_from_MP(chemsys)
        entries.append(self)
        if exclusions:
            entries = [e for e in entries if e.name not in exclusions]
            entries = [e for e in entries if e.entry_id not in exclusions]
        return entries

    @staticmethod
    def get_PD_entries_from_MP(chemsys):
        return mpdata.get_entries_in_chemsys(chemsys)

    @staticmethod
    def get_PD_entries_from_preload_file(chemsys):
        """
        Kept for backward compatibility. Caching is now handled by interface_stability.mpdata.
        """
        return mpdata.get_entries_in_chemsys(chemsys)

    def get_decomp_entries_and_e_above_hull(self, entries=None, exclusions=None, trypreload=None):
        if not entries:
            entries = self.get_PD_entries(exclusions=exclusions, trypreload=trypreload)
        pd = PhaseDiagram(entries)
        decomp_entries, hull_energy = pd.get_decomp_and_e_above_hull(self, allow_negative=True)
        return decomp_entries, hull_energy

    def stabilize(self, entries=None):
        """
        Stabilize an entry by putting it on the convex hull
        (1e-8 eV below it, so the phase diagram counts it as stable).
        """
        entries = [e for e in entries if e is not self] if entries else None
        if not entries:
            entries = self.get_PD_entries_from_MP(self.chemsys)
        decomp_entries, hull_energy = self.get_decomp_entries_and_e_above_hull(entries=entries)
        add_energy_adjustment(self, -(hull_energy * self.composition.num_atoms + 1e-8), name="Stabilization")
        return None

    def energy_correction(self, e):
        """
        Correction term is applied by per atom.
        """
        add_energy_adjustment(self, e * self.composition.num_atoms, name="Manual energy correction")
        return None

    def get_printable_PE_data_in_pd(self, entries=None):
        decomp, hull_e = self.get_decomp_entries_and_e_above_hull(entries=entries)
        output = ['-' * 60]
        PE = list(decomp.keys())
        output.append("Reduced formula of the given composition: " + self.composition.reduced_formula)
        output.append("Calculated phase equilibria: " + "\t".join(i.name for i in PE))
        rxn = ComputedReaction([self], PE)
        rxn.normalize_to(self.composition.reduced_composition)
        output.append(str(rxn))
        output.append('-' * 60)
        string = '\n'.join(output)
        return string

    def GPComp(self, chempot):
        """
        Non-open element composition, which excluded the open element part.
        """
        open_els = {element_symbol(el) for el in chempot}
        GPComp = Composition({el: amt for el, amt in self.composition.items() if el.symbol not in open_els})
        return GPComp

    def get_gppd_entries(self, chempot, exclusions=None, trypreload=False):
        return self.get_PD_entries(sup_el=list(chempot.keys()), exclusions=exclusions, trypreload=trypreload)

    def get_decomposition_in_gppd(self, chempot, entries=None, exclusions=None, trypreload=False):
        """
        :param chempot: {open element: chemical potential in eV, referenced to the pure element}
        :return: (decomposition entries, ComputedReaction). The open element entries in the reaction
            have their energy set to the given chemical potential, so the reaction energy is the
            grand potential change.
        """
        chempot = {element_symbol(el): mu for el, mu in chempot.items()}
        gppd_entries = entries if entries \
            else self.get_gppd_entries(chempot, exclusions=exclusions, trypreload=trypreload)
        pd = PhaseDiagram(gppd_entries)
        stable_entries = list(pd.stable_entries)
        el_ref = {el: pd.el_refs[Element(el)].energy_per_atom for el in chempot}
        chempot_vaspref = {el: chempot[el] + el_ref[el] for el in chempot}

        # Copies of the open element references, with energy equal to the chemical potential
        open_el_entries = [shifted_copy(pd.el_refs[Element(el)], chempot[el]) for el in chempot]

        GPPD = GrandPotentialPhaseDiagram(stable_entries, {Element(el): mu for el, mu in chempot_vaspref.items()})
        GPComp = self.GPComp(chempot)
        decomp_GP_entries = GPPD.get_decomposition(GPComp)
        decomp_entries = [gpe.original_entry for gpe in decomp_GP_entries]
        rxn = ComputedReaction([self] + open_el_entries, decomp_entries)
        rxn.normalize_to(self.composition)
        return decomp_entries, rxn

    def get_printable_PE_and_decomposition_in_gppd(self, chempot, entries=None, exclusions=None, trypreload=False):
        oes = list(chempot.keys())
        output = ['-' * 60, "Reduced formula of the given composition: " + self.composition.reduced_formula]
        for oe in oes:
            output.append("Open element : " + element_symbol(oe))
            output.append("Chemical potential: {:.5g} eV referenced to pure phase".format(chempot[oe]))
        output.append('-' * 60)
        decomp_entries, rxn = self.get_decomposition_in_gppd(chempot, entries=entries, exclusions=exclusions,
                                                             trypreload=trypreload)
        formula = self.composition.reduced_composition
        rxn.normalize_to(formula)
        rxn_e = round(rxn.calculated_reaction_energy, 5)
        output.append("Reaction:" + str(rxn))
        output.append("Reaction energy: {:.5g} eV per {}".format(rxn_e, formula.reduced_formula))
        output.append('-' * 60)
        string = '\n'.join(output)
        return string

    def get_phase_evolution_profile(self, oe, allowpmu=False, entries=None, exclusions=None):
        """
        The phase evolution of this composition when open to element oe, from high to low chemical potential.
        Chemical potentials in the profile are absolute (same energy scale as the entries).
        The 'element_reference' of each stage is the pure element entry (not shifted).

        :param allowpmu: also include chemical potentials above the pure element (positive mu)
        """
        oe = Element(oe)
        pd_entries = entries if entries else self.get_PD_entries(sup_el=[oe], exclusions=exclusions)
        offset = 30 if allowpmu else 0
        originals = {}
        shifted_entries = []
        for e in pd_entries:
            if offset and e.composition.is_element and oe in e.composition:
                # Raise the pure element energy so that positive mu (ref. to the element) becomes reachable
                new = shifted_copy(e, offset)
                originals[id(new)] = e
                shifted_entries.append(new)
            else:
                shifted_entries.append(e)
        pd = PhaseDiagram(shifted_entries)
        evolution_profile = pd.get_element_profile(oe, self.composition.reduced_composition)
        el_ref = evolution_profile[0]['element_reference']
        el_ref = originals.get(id(el_ref), el_ref)
        for stage in evolution_profile:
            stage['element_reference'] = el_ref
        evolution_profile[0]['chempot'] -= offset
        return evolution_profile

    def get_stability_window(self, oe, allowpmu=False, entries=None):
        """
        The chemical potential range (ref. to the pure element) where this phase is stable
        against gain or loss of element oe.
        :return: (mu_high, mu_low). mu_low is None if there is no lower bound, and
            (None, None) if the phase is not stable at any chemical potential.
        """
        profile = self.get_phase_evolution_profile(oe=oe, allowpmu=allowpmu, entries=entries)
        chempots = [_['chempot'] for _ in profile]
        evolutions = [_['evolution'] for _ in profile]
        index = min(range(len(evolutions)), key=lambda i: abs(evolutions[i]))

        if abs(evolutions[index]) < 1e-8:
            ref = profile[0]['element_reference'].energy_per_atom
            if index < len(profile) - 1:
                return (chempots[index] - ref, chempots[index + 1] - ref)
            else:
                return (chempots[index] - ref, None)
        else:
            return (None, None)

    def get_evolution_phases_table_string(self, open_el, pure_el_ref, PE_list, oe_amt_list, mu_trans_list, allowpmu):
        mu_h_list = ['inf' if allowpmu else 0] + mu_trans_list
        mu_l_list = mu_h_list[1:] + ['-inf']
        df = pandas.DataFrame()
        df['mu_high (eV)'] = mu_h_list
        df['mu_low (eV)'] = mu_l_list
        df['d(n_{})'.format(element_symbol(open_el))] = oe_amt_list
        PE_names = []
        rxns = []
        for PE in PE_list:
            rxn = ComputedReaction([self, pure_el_ref], PE)
            rxn.normalize_to(self.composition.reduced_composition)
            PE_names.append(', '.join(sorted([_.name for _ in PE])))
            rxns.append(str(rxn))
        df['Phase equilibria'] = PE_names
        df['Reaction'] = rxns
        print_df = df.to_string(index=False, float_format='{:,.2f}'.format, justify='center')
        return print_df

    def get_rxn_e_table_string(self, pure_el_ref, open_el, PE_list, oe_amt_list, mu_trans_list, plot_rxn_e,
                               save_path=None):
        """
        :param plot_rxn_e: whether to plot the reaction energy
        :param save_path: if given, save the plot to this file instead of showing it
        """
        open_el = element_symbol(open_el)
        neg_flag = (max(mu_trans_list) > 1e-6) if mu_trans_list else False
        rxn_trans_list = [mu for mu in mu_trans_list]
        rxn_e_list = []
        ext = 0.2
        rxn_trans_list = [rxn_trans_list[0] + ext] + rxn_trans_list if neg_flag else [0] + rxn_trans_list
        for data in zip(oe_amt_list, PE_list, rxn_trans_list):
            oe_amt, PE, ext_miu = data
            rxn = ComputedReaction([self, pure_el_ref], PE)
            rxn.normalize_to(self.composition.reduced_composition)
            rxn_e_list.append(rxn.calculated_reaction_energy - oe_amt * ext_miu)
        rxn_trans_list = rxn_trans_list + [rxn_trans_list[-1] - ext]
        rxn_e_list = rxn_e_list + [rxn_e_list[-1] + ext * oe_amt_list[-1]]
        rxn_e_list = [e / self.composition.reduced_composition.num_atoms for e in rxn_e_list]
        df = pandas.DataFrame()
        df["miu_{} (eV)".format(open_el)] = rxn_trans_list
        df["Rxn energy (eV/atom)"] = rxn_e_list

        if plot_rxn_e:
            plt.figure(figsize=(8, 6))
            ax = plt.gca()
            ax.invert_xaxis()
            ax.axvline(0, linestyle='--', color='k', linewidth=0.5, zorder=1)
            ax.plot(rxn_trans_list, rxn_e_list, '-', linewidth=1.5, color='cornflowerblue', zorder=3)
            ax.scatter(rxn_trans_list[1:-1], rxn_e_list[1:-1], edgecolors='cornflowerblue', facecolors='w',
                       linewidth=1.5, s=50, zorder=4)
            ax.set_xlabel('Chemical potential ref. to {}'.format(open_el))
            ax.set_ylabel('Reaction energy (eV/atom)')
            ax.set_xlim([float(rxn_trans_list[0]), float(rxn_trans_list[-1])])
            if save_path:
                plt.savefig(save_path, bbox_inches='tight')
                plt.close()
            else:
                plt.show()
        print_df = df.to_string(index=False, float_format='{:,.2f}'.format, justify='center')

        return print_df

    def get_printable_evolution_profile(self, open_el, entries=None, plot_rxn_e=True, allowpmu=False, save_path=None):
        evolution_profile = self.get_phase_evolution_profile(open_el, entries=entries, allowpmu=allowpmu)

        PE_list = [list(stage['entries']) for stage in evolution_profile]
        oe_amt_list = [stage['evolution'] for stage in evolution_profile]
        pure_el_ref = evolution_profile[0]['element_reference']

        miu_trans_list = [stage['chempot'] for stage in evolution_profile][1:]  # The first chempot is always useless
        miu_trans_list = sorted(miu_trans_list, reverse=True)
        miu_trans_list = [miu - pure_el_ref.energy_per_atom for miu in miu_trans_list]

        table1 = self.get_evolution_phases_table_string(open_el, pure_el_ref, PE_list, oe_amt_list, miu_trans_list,
                                                        allowpmu)
        table2 = self.get_rxn_e_table_string(pure_el_ref, open_el, PE_list, oe_amt_list, miu_trans_list, plot_rxn_e,
                                             save_path=save_path)

        output = ['-' * 60, "Reduced formula of the given composition: " + self.composition.reduced_formula,
                  '\n === Evolution Profile ===', str(table1), '\n === Reaction energy ===', str(table2),
                  'Note:\nChemical potential referenced to element phase.',
                  'Reaction energy is normalized to per atom of the given composition.']
        string = '\n'.join(output)
        return string

    def get_vc_plot_data(self, open_el, valence=None, entries=None, allowpmu=True):
        open_el = element_symbol(open_el)
        if valence:
            ioncharge = valence
        else:
            if open_el not in COMMON_WORKING_IONS:
                raise ValueError('Working ion {} not supported. You can provide charge manually'.format(open_el))
            else:
                ioncharge = COMMON_WORKING_IONS[open_el]

        evolution_profile = self.get_phase_evolution_profile(open_el, entries=entries, allowpmu=allowpmu)
        oe_list = []
        v_list = []
        for i in range(len(evolution_profile)):
            step = evolution_profile[-i - 1]
            oe_content = step['evolution']
            miu_vasp = step['chempot']
            oe_list.append(oe_content)
            v_list.append(miu_vasp)
        v_ref = v_list[-1]
        v_list = [-i + v_ref for i in v_list]
        v_list = [v / ioncharge for v in v_list]

        return oe_list, v_list

    def get_printable_vc_plot_data(self, open_el, oe_list, v_list):
        open_el = element_symbol(open_el)
        df = pandas.DataFrame()
        oes, vs = [], []
        for i in range(len(oe_list) - 1):
            oes.append(oe_list[i])
            oes.append(oe_list[i + 1])
            vs.append(v_list[i])
            vs.append(v_list[i])
        df["d n({})".format(open_el)] = oes
        df["Voltage ref. to {} (V)".format(open_el)] = vs
        print_df = df.to_string(index=False, float_format='{:,.2f}'.format, justify='center')
        return print_df

    def get_voltage_profile_plot(self, open_el, oe_list, v_list, valence):
        open_el = element_symbol(open_el)
        X, Y = [], []
        for i in range(len(oe_list) - 1):
            X += [oe_list[i], oe_list[i + 1]]
            Y += [v_list[i], v_list[i]]
        fig, ax = plt.subplots(1, 1)
        plt.plot(X, Y)
        ylabel = 'Potential ref. to {} / '.format(open_el)
        sup = r'${}^{{{}+}}$'.format(open_el, valence) if valence > 1 else r'${}^{{+}}$'.format(open_el)
        ax.set_ylabel(ylabel + sup)
        s1 = re.sub("([0-9.]+)", "_{\\1}", self.name)
        formula = r'$\mathregular{' + s1 + '}$'
        ax.set_xlabel(r'$\Delta$n({}) per {}'.format(open_el, formula))
        ax.legend([formula])
        if min(ax.get_ylim()) < 0:
            ax.axhline(0, linestyle='--', color='k', linewidth=0.5, zorder=1)
        else:
            ax.set_ylim(bottom=0)
        return plt
