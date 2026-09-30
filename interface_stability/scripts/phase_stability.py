import argparse
from pymatgen.core import Composition
from interface_stability import mpdata
from interface_stability.singlephase import VirtualEntry, COMMON_WORKING_IONS


def get_phase_equilibria_from_composition(args):
    """
    Provides the phase equilibria of a phase with given composition
    """
    comp = Composition(args.composition)
    entry = VirtualEntry.from_composition(comp)
    print(entry.get_printable_PE_data_in_pd())
    return 0


def get_phase_equilibria_and_decomposition_energy_under_mu_from_composition(args):
    """
    Provide the phase equilibria and decomposition energy when open to one element with given miu
    Chemical potential is referenced to pure phase of open element.
    """
    comp = Composition(args.composition)
    chempot = {args.open_element: args.chemical_potential}
    entry = VirtualEntry.from_composition(comp)
    entries = entry.get_gppd_entries(chempot)
    entry.stabilize(entries=entries)
    print(entry.get_printable_PE_and_decomposition_in_gppd(chempot, entries=entries))
    return 0


def get_phase_evolution_profile(args):
    """
    Provides the phase equilibria and decomposition energy evolution process of a phase when open to a specific element
    Chemical potential is referenced to pure phase of open element.
    """
    comp = Composition(args.composition)
    entry = VirtualEntry.from_composition(comp)
    oe = args.open_element
    entries = entry.get_PD_entries(sup_el=[oe])
    entry.stabilize(entries=entries)
    print(entry.get_printable_evolution_profile(oe, entries=entries, allowpmu=args.posmu,
                                                plot_rxn_e=not args.noplot, save_path=args.save))
    return 0


def plot_vc(args):
    """
    Get the plot data of voltage profile.
    """
    comp = Composition(args.composition)
    entry = VirtualEntry.from_composition(comp)
    oe = args.open_element
    entries = entry.get_PD_entries(sup_el=[oe])
    entry.stabilize(entries=entries)
    valence = args.valence if args.valence else COMMON_WORKING_IONS.get(oe)
    if not valence:
        raise SystemExit("Working ion {} not supported. Provide its valence with -v".format(oe))
    oe_list, v_list = entry.get_vc_plot_data(oe, valence=valence, entries=entries, allowpmu=args.posmu)
    print(entry.get_printable_vc_plot_data(oe, oe_list, v_list))
    plot = entry.get_voltage_profile_plot(oe, oe_list, v_list, valence)
    if args.save:
        plot.savefig(args.save, bbox_inches='tight')
    elif not args.noplot:
        plot.show()
    return 0


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""
--BRIEF INTRO--
    This script will analyze the stability of a phase with any given input composition
    Either in a closed system (phase diagram) or in a system with an open element(grand potential phase diagram)
    This script works based on several sub-commands with their own options.
    To see the options for the sub-commands, use "phase_stability sub-command -h".
    """, epilog="""
--REMINDER--
    To use this script, you need a Materials Project API key (https://next-gen.materialsproject.org/api):
    MP_API_KEY : [Mandatory] environment variable with your API key
                 (or PMG_MAPI_KEY in ~/.pmgrc.yaml)
    IFS_CACHE_DIR : [Optional] a local directory to cache downloaded entries
                 (or PMG_PD_PRELOAD_PATH in ~/.pmgrc.yaml)
    """)
    parser.add_argument("--thermo-type", default=mpdata.DEFAULT_THERMO_TYPE, choices=mpdata.THERMO_TYPES,
                        help="Materials Project thermo type to use (default: %(default)s)")

    parent_comp_mp = argparse.ArgumentParser(add_help=False)
    parent_comp_mp.add_argument("composition", type=str, help="The composition for analysis")
    parent_oe = argparse.ArgumentParser(add_help=False)
    parent_oe.add_argument("open_element", type=str, help="The open element")
    parent_mu = argparse.ArgumentParser(add_help=False)
    parent_mu.add_argument("chemical_potential", type=float, help="The chemical potential of open element."
                                                                   "Referenced to pure phase")

    parent_posmu = argparse.ArgumentParser(add_help=False)
    parent_posmu.add_argument("-posmu", action='store_true', default=False,
                             help="Allow mu range to go beyond 0 and become positive")

    parent_plot = argparse.ArgumentParser(add_help=False)
    parent_plot.add_argument("--save", type=str, default=None, metavar="FILE",
                             help="Save the figure to FILE (e.g. plot.png) instead of showing it")
    parent_plot.add_argument("--noplot", action='store_true', default=False, help="Do not make a figure")

    subparsers = parser.add_subparsers()

    parser_stability = subparsers.add_parser("stability", parents=[parent_comp_mp],
                                             help="Obtain the phase equilibria of a phase with given composition")
    parser_stability.set_defaults(func=get_phase_equilibria_from_composition)

    parser_evolution = subparsers.add_parser("evolution", parents=[parent_comp_mp, parent_oe, parent_posmu,
                                                                   parent_plot],
                                             help="Obtain the evolution profile at a given composition when open to an element")

    parser_evolution.set_defaults(func=get_phase_evolution_profile)

    parser_mu = subparsers.add_parser("mu", parents=[parent_comp_mp, parent_oe, parent_mu],
                                       help="Obtain the phase equilibria & decomposition energy of a phase with given composition when open to an element")
    parser_mu.set_defaults(func=get_phase_equilibria_and_decomposition_energy_under_mu_from_composition)

    parser_plot_vc = subparsers.add_parser("plotvc", parents=[parent_comp_mp, parent_oe, parent_posmu, parent_plot],
                                           help="Plot the voltage profile of at a given composition")
    parser_plot_vc.add_argument('-v', '--valence', type=int, default=None, help='Valence of Working ion')

    parser_plot_vc.set_defaults(func=plot_vc)

    args = parser.parse_args()

    if hasattr(args, "func"):
        # Bad input (e.g. an unknown element) raises ValueError, and Materials Project errors RuntimeError
        try:
            mpdata.set_thermo_type(args.thermo_type)
            args.func(args)
        except (RuntimeError, ImportError, ValueError) as err:
            raise SystemExit("Error: {}".format(err))
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
