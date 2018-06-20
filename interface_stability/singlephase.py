# coding: utf-8
# Copyright (c) Yifei Mo Group @ University of Maryland, College Park
# Distributed under the terms of the MIT License.


import os
import json
import re
import warnings
import pandas

import matplotlib.pyplot as plt
from matplotlib import rc
from monty.json import MontyDecoder, MontyEncoder
from pymatgen import Composition, SETTINGS, Element, MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.entries.computed_entries import ComputedEntry

__author__ = "Yizhou Zhu"
__copyright__ = ""
__version__ = "2.2"
__maintainer__ = "Yizhou Zhu"
__email__ = "yizhou.zhu@gmail.com"
__status__ = "Production"
__date__ = "Jun 10, 2018"

PD_PRELOAD_PATH = SETTINGS.get("PMG_PD_PRELOAD_PATH")
# if PD_PRELOAD_PATH is None:
#     trypreload = False

plt.rcParams['mathtext.default'] = 'regular'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica'], 'size': 15})


class VirtualEntry(ComputedEntry):
    def __init__(self, composition, energy, attribute=None, name=None):
        super(VirtualEntry, self).__init__(Composition(composition), energy, attribute=attribute)
        if name:
            self.name = name

    @classmethod
    def from_composition(cls, comp, energy=0, name=None):
        attribute = 'An initialized virtual entry'
        return cls(Composition(comp), energy, attribute, name=name)

    @classmethod
    def from_mixing(cls, mixing_dict):
        attribute = 'An virtual entry based on linear combination of parent entries\nParents:\n'
        comp = Composition("")
        energy = 0
        for i in mixing_dict.keys():
            comp += Composition({el: i.composition[el] * mixing_dict[i] for el in i.composition.keys()})
            energy += mixing_dict[i] * i.energy
            attribute += str(i.composition) + "\tEnergy:{}\tAmount:{}\n".format(i.energy, mixing_dict[i])
        return cls(Composition(comp), energy, attribute)

    @classmethod
    def from_mp(cls, criteria):
        attribute = "Entry from MP database"
        entry = cls.get_mp_entry(criteria)
        return cls(entry.composition, energy=entry.energy, name=entry.name, attribute=attribute)

    @staticmethod
    def get_mp_entry(criteria):
        """
        Here always return the lowest energy among all polymorphs.
        Criteria can be a formula or an mp-id
        """
        with MPRester() as m:
            entries = m.get_entries(criteria)
        entries = sorted(entries, key=lambda e: e.energy_per_atom)
        if len(entries) == 0:
            raise ValueError("MP doesn't have any entry that matches the given formula/MP-id!")
        return entries[0]

    @property
    def chemsys(self):
        return [_.symbol for _ in self.composition.elements]

    def get_PD_entries(self, sup_el=None, exclusions=None, trypreload=False):
        """
        :param sup_el: a list for extra element dimension, using str format
        :param exclusions: a list of manually exclusion entries, can use entry name or mp_id
        :param trypreload: If try to reload from cached search results.
        Warning: if set to True, the return result may not be consistent with the updated MP database.
        :return: all related entries to construct phase diagram.
        """

        chemsys = self.chemsys + sup_el if sup_el else self.chemsys
        chemsys = list(set(chemsys))

        if trypreload:
            entries = self.get_PD_entries_from_preload_file(chemsys)
        else:
            entries = self.get_PD_entries_from_MP(chemsys)
        entries.append(self)
        if exclusions:
            entries = [e for e in entries if e.name not in exclusions]
            entries = [e for e in entries if e.entry_id not in exclusions]
        return entries

    @staticmethod
    def get_PD_entries_from_MP(chemsys):
        with MPRester() as m:
            entries = m.get_entries_in_chemsys(chemsys)

        return entries

    @staticmethod
    def get_PD_entries_from_preload_file(chemsys):
        """
        If you use this method, the results may be incompatible with the most updated MP database.
        """
        if PD_PRELOAD_PATH is None:
            warnings.warn("\nYou are trying load locally cached entries. "
                          "\nYou should set up a valid folder for local cache, "
                          "\nand put the path as PMG_PD_PRELOAD_PATH in ~/.pmgrc.yaml")
        if not os.path.isdir(PD_PRELOAD_PATH):
            warnings.warn("\nPMG_PD_PRELOAD_PATH is not a valid folder path."
                          "\nPlease reset PMG_PD_PRELOAD_PATH in ~/.pmgrc.yaml")

        el_list = [x.symbol for x in sorted(set(chemsys))]
        load_path = os.path.join(PD_PRELOAD_PATH, "_".join(el_list) + "_Entries.json")
        try:
            with open(load_path) as f:
                entries = json.load(f, cls=MontyDecoder)
        except (IOError, EOFError):
            entries = VirtualEntry.get_PD_entries_from_MP(chemsys)
            with open(load_path, 'w') as f:
                json.dump(entries, f, cls=MontyEncoder)
        return entries

    def get_decomp_entries_and_e_above_hull(self, entries=None, exclusions=None, trypreload=None):
        if not entries:
            entries = self.get_PD_entries(exclusions=exclusions, trypreload=trypreload)
        pd = PhaseDiagram(entries)
        decomp_entries, hull_energy = pd.get_decomp_and_e_above_hull(self)
        return decomp_entries, hull_energy

    def stabilize(self, entries=None):
        """
        Stabilize an entry by putting it on the convex hull
        """
        self.attribute = "A virtual entry stabilized with \"on the hull\" energy"
        decomp_entries, hull_energy = self.get_decomp_entries_and_e_above_hull(entries=entries)
        self.correction -= (hull_energy * self.composition.num_atoms + 1e-8)
        return None

    def energy_correction(self, e):
        """
        Correction term is applied by per atom.
        """
        self.correction += e * self.composition.num_atoms
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
        GPComp = Composition({el: amt for el, amt in self.composition.items() if el.symbol not in chempot.keys()})
        return GPComp

    def get_gppd_entries(self, chempot, exclusions=None, trypreload=False):
        return self.get_PD_entries(sup_el=list(chempot.keys()), exclusions=exclusions, trypreload=trypreload)

    def get_decomposition_in_gppd(self, chempot, entries=None, exclusions=None, trypreload=False):
        gppd_entries = entries if entries \
            else self.get_gppd_entries(chempot, exclusions=exclusions, trypreload=trypreload)
        pd = PhaseDiagram(gppd_entries)
        gppd_entries = pd.stable_entries
        open_el_entries = [_ for _ in gppd_entries if
                           _.is_element and _.composition.elements[0].symbol in chempot.keys()]
        el_ref = {_.composition.elements[0].symbol: _.energy_per_atom for _ in open_el_entries}
        chempot_vaspref = {_: chempot[_] + el_ref[_] for _ in chempot}
        for open_entry in open_el_entries:
            open_entry.correction += chempot_vaspref[open_entry.composition.elements[0].symbol]

        GPPD = GrandPotentialPhaseDiagram(gppd_entries, chempot_vaspref)
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
            output.append("Open element : " + oe)
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

    def get_phase_evolution_profile(self, oe, allowpmu=False, entries=None,exclusions=None):
        pd_entries = entries if entries else self.get_PD_entries(sup_el=[oe],exclusions=exclusions)
        offset = 30 if allowpmu else 0
        for e in pd_entries:
            if e.composition.is_element and oe in e.composition.keys():
                e.correction += offset * e.composition.num_atoms
        pd = PhaseDiagram(pd_entries)
        evolution_profile = pd.get_element_profile(oe, self.composition.reduced_composition)
        el_ref = evolution_profile[0]['element_reference']
        el_ref.correction -= el_ref.composition.num_atoms * offset
        evolution_profile[0]['chempot'] -= offset
        return evolution_profile


    def get_stability_window(self,oe,allowpmu=False, entries=None):
        profile = self.get_phase_evolution_profile(oe=oe,allowpmu=allowpmu,entries=entries)
        chempots = [_['chempot'] for _ in profile]
        evolutions = [_['evolution'] for _ in profile]
        index = evolutions.index(sorted(evolutions,key=lambda x: abs(x))[0])

        if abs(evolutions[index]) < 1e-8:
            ref = profile[0]['element_reference'].energy_per_atom
            if index < len(profile)-1:
                return (chempots[index]-ref,chempots[index+1]-ref)
            else:
                return (chempots[index]-ref,None)
        else:
            return (None, None)



        #index = chempots.index([evolution])
        # print (evolutions[index],chempots[index])

        #return



    def get_evolution_phases_table_string(self, open_el, pure_el_ref, PE_list, oe_amt_list, mu_trans_list, allowpmu):
        if not allowpmu:
            mu_h_list = [0] + mu_trans_list
        mu_l_list = mu_h_list[1:] + ['-inf']
        df = pandas.DataFrame()
        df['mu_high (eV)'] = mu_h_list
        df['mu_low (eV)'] = mu_l_list
        df['d(n_{})'.format(open_el)] = oe_amt_list
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

    def get_rxn_e_table_string(self, pure_el_ref, open_el, PE_list, oe_amt_list, mu_trans_list, plot_rxn_e):
        neg_flag = (max(mu_trans_list) > 1e-6)
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
        rxn_e_list = [e / self.composition.num_atoms for e in rxn_e_list]
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
            plt.show()
        print_df = df.to_string(index=False, float_format='{:,.2f}'.format, justify='center')

        return print_df

    def get_printable_evolution_profile(self, open_el, entries=None, plot_rxn_e=True, allowpmu=False):
        evolution_profile = self.get_phase_evolution_profile(open_el, entries=entries, allowpmu=allowpmu)

        PE_list = [list(stage['entries']) for stage in evolution_profile]
        oe_amt_list = [stage['evolution'] for stage in evolution_profile]
        pure_el_ref = evolution_profile[0]['element_reference']

        miu_trans_list = [stage['chempot'] for stage in evolution_profile][1:]  # The first chempot is always useless
        miu_trans_list = sorted(miu_trans_list, reverse=True)
        miu_trans_list = [miu - pure_el_ref.energy_per_atom for miu in miu_trans_list]

        table1 = self.get_evolution_phases_table_string(open_el, pure_el_ref, PE_list, oe_amt_list, miu_trans_list,
                                                        allowpmu)
        table2 = self.get_rxn_e_table_string(pure_el_ref, open_el, PE_list, oe_amt_list, miu_trans_list, plot_rxn_e)

        output = ['-' * 60, "Reduced formula of the given composition: " + self.composition.reduced_formula,
                  '\n === Evolution Profile ===', str(table1), '\n === Reaction energy ===', str(table2),
                  'Note:\nChemical potential referenced to element phase.',
                  'Reaction energy is normalized to per atom of the given composition.']
        string = '\n'.join(output)
        return string

    def get_vc_plot_data(self, open_el, valence=None, entries=None, allowpmu=True):
        common_working_ion = {Element('Li'): 1, Element('Na'): 1, Element('K'): 1, Element('Mg'): 2, Element('Ca'): 2,
                              Element('Zn'): 2, Element('Al'): 3}
        if valence:
            ioncharge = valence
        else:
            if open_el not in common_working_ion.keys():
                raise ValueError('Working ion {} not supported. You can provide charge manually'.format(open_el.symbol))
            else:
                ioncharge = common_working_ion[open_el]

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
        X, Y = [], []
        for i in range(len(oe_list) - 1):
            X += [oe_list[i], oe_list[i + 1]]
            Y += [v_list[i], v_list[i]]
        fig, ax = plt.subplots(1, 1)
        plt.plot(X, Y)
        ylabel = 'Potential ref. to {} / '.format(open_el, open_el, valence)
        sup = r'${}^{{{}+}}$'.format(open_el.symbol, valence) if valence > 1 else r'${}^{{+}}$'.format(open_el)
        ax.set_ylabel(ylabel + sup)
        s1 = re.sub("([0-9]+)", "_{\\1}", self.name)
        formula = '$\mathregular{' + s1 + '}$'
        ax.set_xlabel('$\Delta$n({}) per {}'.format(open_el, formula))
        ax.legend([formula])
        if min(ax.get_ylim()) < 0:
            ax.axhline(0, linestyle='--', color='k', linewidth=0.5, zorder=1)
        else:
            ax.set_ylim(bottom=0)
        return plt
