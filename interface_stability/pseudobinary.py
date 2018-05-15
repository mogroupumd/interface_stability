# coding: utf-8
# Copyright (c) Mogroup  @ University of Maryland, College Park
# Distributed under the terms of the MIT License.

import pandas
from pymatgen import Composition, Element
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, GrandPotPDEntry
from pymatgen.analysis.reaction_calculator import ComputedReaction, ReactionError
from interface_stability.singlephase import VirtualEntry


__author__ = "Yizhou Zhu"
__copyright__ = ""
__version__ = "2.2"
__maintainer__ = "Yizhou Zhu"
__email__ = "yizhou.zhu@gmail.com"
__status__ = "Production"
__date__ = "Jun 10, 2018"


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
            entries = entry_mix.get_PD_entries(sup_el=sup_el)
        entries += [entry1, entry2]
        self.PDEntries = entries
        self.PD = PhaseDiagram(entries)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

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
        return self.get_printed_profile(self.gppd_mixing(chempots, gppd_entries=gppd_entries))

    def get_printed_profile(self, profile):
        """
        A general function to generate printable table strings for pseudo-binary mixing results
        """
        output = ['\n ===  Pseudo-binary evolution profile  === ']
        df = pandas.DataFrame()
        rxn_e = []
        mutual_rxn_e = []
        E0 = -profile[0][1][1]
        E1 = -profile[-1][1][1]

        x1s, x2s, es, mes, pes = [], [], [], [], []

        for item in profile:
            ratio, (decomp, e) = item
            x1s.append(1-ratio)
            x2s.append(ratio)
            es.append(-e*1000)
            mes.append((-e - ratio * E1 - (1 - ratio) * E0) * 1000)
            pes.append(", ".join([x.name for x in decomp]))
            rxn_e.append(-e)
            mutual_rxn_e.append(((-e - ratio * E1 - (1 - ratio) * E0) * 1000))
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
        x is the ratio (both entry norm. to 1 atom/fu(w/o open element) ) or each entry
        """
        open_el = list(chempots.keys())[0]
        el_ref = VirtualEntry.get_mp_entry(open_el)
        chempots[open_el] = chempots[open_el] + el_ref.energy_per_atom
        gppd_entry1 = GrandPotPDEntry(self.entry1, {Element[_]: chempots[_] for _ in chempots})
        gppd_entry2 = GrandPotPDEntry(self.entry2, {Element[_]: chempots[_] for _ in chempots})
        if not gppd_entries:
            gppd_entries = self.get_gppd_entries(open_el)
        gppd = GrandPotentialPhaseDiagram(gppd_entries, chempots)
        profile = get_full_evolution_profile(gppd, gppd_entry1, gppd_entry2, 0, 1)
        cleaned = clean_profile(profile)
        return cleaned

    def get_gppd_entries(self, open_el):
        if open_el in (self.entry1.composition + self.entry2.composition).keys():
            gppd_entries = self.PDEntries
        else:
            comp = self.entry1.composition + self.entry2.composition + Composition(open_el.symbol)
            gppd_entries = VirtualEntry.from_composition(comp).get_GPPD_entries(open_el)
        return gppd_entries

    def get_gppd_transition_chempots(self, open_el, gppd_entries=None):
        """
        This is to get all possible transition chemical potentials from PD (rather than GPPD)
        Still use pure element ref.
        # May consider supporting negative miu in the future
        """
        if not gppd_entries:
            gppd_entries = self.get_gppd_entries(open_el)
        pd = PhaseDiagram(gppd_entries)
        vaspref_mius = pd.get_transition_chempots(Element(open_el))
        el_ref = VirtualEntry.get_mp_entry(open_el)

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
        :param verbose: whether to prune the PE result table
        :return: a printable string of screening results
        """
        mu_lo, mu_hi = sorted([mu_lo, mu_hi])
        miu_E_candidates = [miu for miu in self.get_gppd_transition_chempots(open_el) if
                            (miu - mu_lo) * (miu - mu_hi) <= 0]
        miu_E_candidates = [mu_hi] + miu_E_candidates + [mu_lo]
        duplicate_index = []
        if not gppd_entries:
            gppd_entries = self.get_gppd_entries(open_el)

        for i in range(1, len(miu_E_candidates) - 1):
            miu_left = (miu_E_candidates[i] + miu_E_candidates[i - 1]) / 2.0
            miu_right = (miu_E_candidates[i] + miu_E_candidates[i + 1]) / 2.0
            profile_left = self.gppd_mixing({open_el: miu_left}, gppd_entries)
            profile_right = self.gppd_mixing({open_el: miu_right}, gppd_entries)
            if judge_same_decomp(profile_left, profile_right):
                duplicate_index.append(i)
        miu_E_candidates = [miu_E_candidates[i] for i in range(len(miu_E_candidates)) if i not in duplicate_index]

        mu_hi, miu_low, PE = [], [], []
        mu_list, E_mutual_list, E_total_list = [], [], []

        for i in range(1, len(miu_E_candidates)):
            miu = (miu_E_candidates[i] + miu_E_candidates[i - 1]) / 2.0
            profile = self.gppd_mixing({open_el: miu}, gppd_entries)
            E0 = -profile[0][1][1]
            E1 = -profile[-1][1][1]
            min_mutual = min(profile, key=lambda step: (-step[1][1] - step[0] * E1 - (1 - step[0]) * E0))
            mu_hi.append(miu_E_candidates[i - 1])
            miu_low.append(miu_E_candidates[i])
            PE.append(", ".join(sorted([x.name for x in min_mutual[1][0]])))
        for i in range(len(miu_E_candidates)):
            profile_transition = self.gppd_mixing({open_el: miu_E_candidates[i]}, gppd_entries)
            E0 = -profile_transition[0][1][1]
            E1 = -profile_transition[-1][1][1]
            min_mutual_transition = min(profile_transition,
                                        key=lambda step: (-step[1][1] - step[0] * E1 - (1 - step[0]) * E0))
            mu_list.append(miu_E_candidates[i])

            E_mutual_list.append(
                (-min_mutual_transition[1][1] - min_mutual_transition[0] * E1 - (1 - min_mutual_transition[0]) * E0))
            E_total_list.append(-min_mutual_transition[1][1])

        to_be_hidden = []
        if not verbose:
            for i in range(1, len(PE)):
                if PE[i] == PE[i - 1]:
                    to_be_hidden.append(i)

        mu_hi_display_list = [mu_hi[k] for k in range(len(mu_hi)) if k not in to_be_hidden]
        mu_low_display_list = mu_hi_display_list[1:] + [mu_lo]
        PE_display_list = [PE[k] for k in range(len(PE)) if k not in to_be_hidden]

        df1 = pandas.DataFrame()
        df2 = pandas.DataFrame()

        df1['mu_low'] = mu_hi_display_list
        df1['mu_high'] = mu_low_display_list
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
    entry_left = get_mix_entry({entry1: x1, entry2: 1 - x1})
    entry_right = get_mix_entry({entry1: x2, entry2: 1 - x2})
    (decomp1, h1) = pd.get_decomp_and_e_above_hull(entry_left)
    (decomp2, h2) = pd.get_decomp_and_e_above_hull(entry_right)
    decomp1 = set(decomp1.keys())
    decomp2 = set(decomp2.keys())
    evolution_profile[x1] = (decomp1, h1)
    evolution_profile[x2] = (decomp2, h2)

    if decomp1 == decomp2:
        return evolution_profile

    intersect = decomp1 & decomp2
    if len(intersect) > 0:
        # This is try to catch a single transition point
        try:
            rxn = ComputedReaction([entry_left, entry_right], list(intersect))
            if not {entry_left, entry_right} < set(rxn.all_entries):
                return evolution_profile

            c1 = rxn.coeffs[rxn.all_entries.index(entry_left)]
            c2 = rxn.coeffs[rxn.all_entries.index(
                entry_right)]  # I know this is tedious but this is the only way I found that works..
            x = (c1 * x1 + c2 * x2) / (c1 + c2)
            if c1 * c2 == 0:
                return evolution_profile
            entry_mid = VirtualEntry.from_mixing({entry_left: c1 / (c1 + c2), entry_right: c2 / (c1 + c2)})
            h_mid = pd.get_decomp_and_e_above_hull(entry_mid)[1]
            evolution_profile[x] = (intersect, h_mid)
            return evolution_profile
        except ReactionError:
            pass

    x_mid = (x1 + x2) / 2.0
    entry_mid = get_mix_entry({entry1: 0.5, entry2: 0.5})
    (decomp_mid, h_mid) = pd.get_decomp_and_e_above_hull(entry_mid)
    decomp_mid = set(decomp_mid.keys())
    evolution_profile[x_mid] = (decomp_mid, h_mid)
    part1 = get_full_evolution_profile(pd, entry1, entry2, x1, x_mid)
    part2 = get_full_evolution_profile(pd, entry1, entry2, x_mid, x2)
    evolution_profile.update(part1)
    evolution_profile.update(part2)
    return evolution_profile


def clean_profile(evolution_profile):
    """
    This function is to clean the calculated profile from binary search. Redundant trial results are pruned out,
    with only transition points left.
    """
    raw_data = list(evolution_profile.items())
    raw_data.sort()
    clean_set = [raw_data[0]]
    for i in range(1, len(raw_data)):
        x, (decomp, h) = raw_data[i]
        x_cpr, (decomp_cpr, h_cpr) = clean_set[-1]
        if set(decomp_cpr) <= set(decomp):
            continue
        else:
            clean_set.append(raw_data[i])
    return clean_set


def get_mix_entry(mix_dict):
    """
    Mixing PDEntry or GrandPotEntry for the binary search algorithm.
    """
    entry1, entry2 = mix_dict.keys()
    x1, x2 = mix_dict[entry1], mix_dict[entry2]
    if type(entry1) == GrandPotPDEntry:
        mid_ori_entry = VirtualEntry.from_mixing({entry1.original_entry: x1, entry2.original_entry: x2})
        return GrandPotPDEntry(mid_ori_entry, entry1.chempots)
    else:
        return VirtualEntry.from_mixing({entry1: x1, entry2: x2})
