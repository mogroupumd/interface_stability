#!/usr/bin/env python3

import argparse
from pymatgen import Composition
from interface_stability.pseudobinary import PseudoBinary
from interface_stability.singlephase import VirtualEntry


def input_handling(args):
    comp1 = Composition(args.composition_1)
    comp2 = Composition(args.composition_2)
    entry1 = VirtualEntry.from_composition(comp1)
    entry2 = VirtualEntry.from_composition(comp2)
    entry1.stabilize()
    entry2.stabilize()
    entry1.energy_correction(args.e1)
    entry2.energy_correction(args.e2)
    return entry1, entry2


def chemical_stability(args):
    entry1, entry2 = input_handling(args)
    print("-" * 100, "\nThe starting phases compositions are ", entry1.name, 'and', entry2.name)
    print("All mixing ratio based on all formula already normalized to ONE atom per fu!")
    pb = PseudoBinary(entry1, entry2)
    print(pb.get_printable_pd_profile())
    return 0


def electrochemical_stability(args):
    entry1, entry2 = input_handling(args)
    oe = args.open_element
    mu = args.chemical_potential
    chempots = {oe: mu}
    print("-" * 100, "\nThe starting phases compositions are ", entry1.name, 'and', entry2.name)
    print("All mixing ratio based on all formula already normalized to ONE atom per fu!")

    print("Chemical potential is miu_{} = {}, using elementary phase as reference.".format(oe, mu))
    print('-' * 60)
    pb = PseudoBinary(entry1, entry2)
    print(pb.get_printable_gppd_profile(chempots))
    return 0


def electrochemical_stability_screening(args):
    entry1, entry2 = input_handling(args)
    oe = args.open_element
    miu_low = args.miu_low
    miu_high = args.miu_high
    pb = PseudoBinary(entry1, entry2)
    print(pb.gppd_scanning(oe, miu_high, miu_low))


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
--BRIEF INTRO--        
    This script will calculate the phase equilibria of a pseudo binary (linear combination of two entries)
    Either in a closed system (pd) or a system with an open element (gppd)
    They reflect the chemical/electrochemical stability of the pseudo binary, respectively.
    This script works based on several sub-commands with their own options. 
    To see the options for the sub-commands, use "pseudo_stability.py sub-command -h".
    """, epilog="""
--REMINDER--
    To use this script, you need to set following variable in ~/.pmgrc.yaml:
    PMG_MAPI_KEY : the API key for MP to fetch data from MP website.
    PMG_PD_PRELOAD_PATH : the local path for saved pickle files.
    """)
    subparsers = parser.add_subparsers()
    parent_comp_mp = argparse.ArgumentParser(add_help=False)
    parent_comp_mp.add_argument("composition_1", type=str, help="The first phase composition of the pseudo-binary")
    parent_comp_mp.add_argument("composition_2", type=str, help="The second phase composition of the pseudo-binary")

    parent_comp_mp.add_argument("-e1", type=float, default=0.0, help="The energy correction of entry1 ref. to hull ")
    parent_comp_mp.add_argument("-e2", type=float, default=0.0, help="The energy correction of entry2 ref. to hull ")

    parent_oe = argparse.ArgumentParser(add_help=False)
    parent_oe.add_argument("open_element", type=str, help="The open element")

    parent_miu = argparse.ArgumentParser(add_help=False)
    parent_miu.add_argument("chemical_potential", type=float, help="The chemical potential of open element."
                                                                   "Default referenced to pure phase, "
                                                                   "ref can be changed with -vaspref")

    parser_pd = subparsers.add_parser("pd", parents=[parent_comp_mp],
                                      help="The chemical stability / phase equilibria info of the pseudo-binary, "
                                           "calculated in PD")
    parser_pd.set_defaults(func=chemical_stability)

    parser_gppd = subparsers.add_parser("gppd", parents=[parent_comp_mp, parent_oe, parent_miu],
                                        help="The electrochemical stability / phase equilibria info of the pseudo-"
                                             "binary, calculated in GPPD")
    parser_gppd.set_defaults(func=electrochemical_stability)

    parser_gppd_screen = subparsers.add_parser("gppd_screen", parents=[parent_comp_mp, parent_oe],
                                               help="The electrochemical stability in a given chemical potential range")
    parser_gppd_screen.add_argument("miu_low", type=float, help="lower chemical potential for gppd screening")
    parser_gppd_screen.add_argument("miu_high", type=float, help="upper chemical potential for gppd screening")
    parser_gppd_screen.set_defaults(func=electrochemical_stability_screening)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
