"""Routines handling the run of ADG."""
from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
from builtins import input
from builtins import range
from adg.tools import reversed_enumerate
import os
import argparse
import shutil
import networkx as nx
import adg.mbpt
import adg.bmbpt
import adg.pbmbpt
import adg.diag


def parse_command_line(cli_args):
    """Return run commands from the Command Line Interface.

    Args:
        cli_args: Command-line arguments submitted with the program.

    Returns:
        Namespace: Appropriate commands to manage the program's run.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="AUTOMATIC DIAGRAM GENERATOR v%s\n\n" % adg.__version__
        + "Generates diagrams at a given order for a "
        + "range of many-body formalisms",
        epilog="If a theory is chosen in batch mode, all flags associated\n"
        + "to other theories will be deactivated to avoid conflicts.\n\n")

    basic_args = parser.add_argument_group(
        title="Basic arguments",
        description="For ADG to run, you need either to run it in\n"
        + "interactive mode or to provide a theory of interest and an order\n"
        + "at which to run the calculation.")
    run_args = parser.add_argument_group(
        title="Run arguments",
        description="Those arguments control the content of the ouput files.\n"
        )
    mbpt_args = parser.add_argument_group(
        title="MBPT-specific arguments",
        description="Arguments available only for MBPT calculations.\n")
    bmbpt_args = parser.add_argument_group(
        title="BMBPT-specific arguments",
        description="Arguments available only for (P)BMBPT calculations.\n")

    basic_args.add_argument(
        "-o", "--order", type=int, nargs='+', choices=list(range(1, 10)),
        help="order of the diagrams (>=1) or N_A, N_B, N_C truncation for BIMSRG")
    basic_args.add_argument(
        "-t", "--theory", type=str, choices=['MBPT', 'BMBPT', 'PBMBPT', 'BIMSRG'],
        help="theory of interest: MBPT, BMBPT, PBMBPT, BIMSRG")
    basic_args.add_argument(
        "-i", "--interactive", action="store_true",
        help="execute ADG in interactive mode")

    bmbpt_args.add_argument(
        "-can", "--canonical", action="store_true",
        help="consider only canonical diagrams")
    bmbpt_args.add_argument(
        "-nobs", "--nbody_observable", type=int, choices=list(range(1, 4)), default=2,
        help="maximal n-body character of the observable [1-3], default = 2")
    bmbpt_args.add_argument(
        "-3NF", "--with_3NF", action="store_true",
        help="use two and three-body forces for (P)BMBPT diagrams")
    bmbpt_args.add_argument(
        "-dt", "--draw_tsds", action="store_true",
        help="draw Time-Structure Diagrams (BMBPT or PBMBPT)")

    run_args.add_argument(
        "-d", "--draw_diags", action="store_true",
        help="draw the diagrams using FeynMF")
    run_args.add_argument(
        "-c", "--compile", action="store_true",
        help="compile the LaTeX output file with PDFLaTeX")

    mbpt_args.add_argument(
        "-cd", "--cd_output", action="store_true",
        help="produce computer-readable output for automated framework (MBPT)")

    args = parser.parse_args(cli_args)

    if (not args.interactive) and ((args.order is None)
                                   or (args.theory is None)):
        print("\nPlease either run the interactive mode, or the batch mode by")
        print("providing the theory and the order for the desired diagrams.\n")
        print("Use 'adg -h' for help.\n")
        exit()

    # Avoid conflicting flags
    if args.theory not in ('BMBPT', 'PBMBPT') and not args.interactive:
        args.canonical = None
        args.with_3NF = None
        args.nbody_observable = 2
        args.draw_tsds = None
    if args.theory != 'MBPT' and not args.interactive:
        args.cd_output = None
    if args.theory != 'BIMSRG' and len(args.order) > 1:
        print("Order of the diagrams has to be a single int. Exiting.")
        exit()

    if args.theory == 'BIMSRG':
        if len(args.order) == 1:
            print("Using the traditional IMSRG truncation: N_A = N_B = N_C.")
            dummy_order = args.order[0]
            args.order = (dummy_order, dummy_order, dummy_order)
        elif len(args.order) == 2:
            print("Truncate commutator and operators differently: N_A = N_B.")
            ope_order, comm_order = args.order[0], args.order[-1]
            args.order = (ope_order, ope_order, comm_order)
        elif len(args.order) == 3:
            print("Using the N_A, N_B, N_C truncation scheme.")
            a_order, b_order, comm_order = args.order[0], args.order[1], args.order[-1]
            args.order = (a_order, b_order, comm_order)
        else:
            print("Truncation order has to be three int at most. Exiting.")
            exit()
        if args.draw_diags and max(args.order) >= 7:
            args.draw_diags = False
            print("Diagram drawing is not supported for BIMSRG above order 6.")
            print("The '-d' flag has been deactivated.")
    else:
        dummy_order = args.order[0]
        args.order = dummy_order

    return args


def interactive_interface(commands):
    """Run the interactive interface mode, return the appropriate commands.

    Args:
        commands (Namespace): Flags for the run management.

    Returns:
        Namespace: Flags initialized through keyboard input.

    """
    theories = ["BMBPT", "MBPT", "PBMBPT", "BIMSRG"]

    commands.theory = input('MBPT, BMBPT, PBMBPT or BIMSRG?\n').upper()
    while commands.theory not in theories:
        print("Invalid theory!")
        commands.theory = input('MBPT, BMBPT, PBMBPT or BIMSRG?\n').upper()

    if commands.theory != 'BIMSRG':
        try:
            commands.order = int(input('Order of the diagrams? [1-9]\n'))
        except ValueError:
            print("Please enter an integer value! Program exiting.")
            exit()
        while commands.order < 1 or commands.order > 9:
            print("Perturbative order too small or too high!")
            commands.order = int(input('Order of the diagrams? [1-9]\n'))
    else:
        commands.order = [0, 0, 0]
        for position, operator in enumerate('A', 'B', 'C'):
            commands.order[position] = get_bimsrg_truncation_order(operator)
        commands.order = tuple(commands.order)

    if commands.theory in ("BMBPT", "PBMBPT"):
        commands.canonical = input(
            "Consider only canonical diagrams? (y/N)").lower() == 'y'
        commands.with_3NF = input(
            "Include three-body forces? (y/N)").lower() == 'y'
        try:
            commands.nbody_observable = int(input(
                "Maximal n-body character of the observable? [1-3]"))
        except ValueError:
            print("Please enter an integer value! Program exiting.")
            exit()
        while commands.nbody_observable < 1 or commands.nbody_observable > 3:
            print("The observable must be 1-body, 2-body or 3-body!")
            commands.nbody_observable = int(input(
                "Maximal n-body character of the observable? [1-3]"))
        commands.draw_tsds = input(
            "Draw time-structure diagrams? (y/N)").lower() == 'y'

    commands.draw_diags = input(
        "Generate diagrams FeynMF instructions in TeX file? (y/N) "
        ).lower() == 'y'

    if commands.theory == "BIMSRG" and commands.draw_diags \
            and max(commands.order) >= 7:
        commands.draw_diags = False
        print("Diagram drawing is not supported for BIMSRG above order 6.")
        print("The '-d' flag has been deactivated.")

    if commands.theory == "MBPT":
        commands.cd_output = input(
            "Produce a CD output file? (y/N) ").lower() == 'y'
    commands.compile = input("Compile pdf? (y/N) ").lower() == 'y'

    return commands


def get_bimsrg_truncation_order(operator):
    """Return the truncation order of a given operator from the user input.

    Args:
        operator (str): The letter corresponding to the operator name.py

    Returns:
        int: The truncation rank of the operator.
    """
    try:
        order = int(input('Order of the %s operator? [1-9]\n' % operator))
    except ValueError:
        print("Please enter an integer value! Program exiting.")
        exit()
    while order < 1 or order > 9:
        print("Order too small or too high!")
        order = int(input('Order of the %s operator? [1-9]\n' % operator))
    return order


def attribute_directory(commands):
    """Create missing directories and return the working directory.

    Args:
        commands (Namespace): Flags for the run management.

    Returns:
        str: Path to the result folder.

    >>> com = argparse.Namespace()
    >>>
    >>> com.theory, com.order = 'BMBPT', 4
    >>> com.with_3NF, com.nbody_observable, com.canonical = False, 2, False
    >>>
    >>> attribute_directory(com)
    'BMBPT/Order-4_2body_observable'
    >>>
    >>> com.theory, com.order = 'BMBPT', 5
    >>> com.with_3NF, com.nbody_observable, com.canonical = True, 3, False
    >>>
    >>> attribute_directory(com)
    'BMBPT/Order-5_3body_observable_with3N'
    >>>
    >>> com.theory, com.order = 'MBPT', 3
    >>> com.with_3NF, com.nbody_observable, com.canonical = False, 2, False
    >>>
    >>> attribute_directory(com)
    'MBPT/Order-3'
    >>>
    >>> com.theory, com.order = 'BIMSRG', (1,2,3)
    >>>
    >>> attribute_directory(com)
    'BIMSRG/Order_1_2_3'

    """
    directory = '%s/' % commands.theory
    if commands.theory != "BIMSRG":
        directory += 'Order-%i' % commands.order
        if commands.canonical:
            directory += '_canonical'
        if commands.theory in ('BMBPT', 'PBMBPT'):
            directory += '_%ibody_observable' % commands.nbody_observable
        if commands.with_3NF:
            directory += '_with3N'
    else:
        directory += 'Order_%i_%i_%i' % commands.order
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(directory+"/Diagrams"):
        os.makedirs(directory+"/Diagrams")
    return directory


def generate_diagrams(commands, id_generator):
    """Return a list with diagrams of the appropriate type.

    Args:
        commands (Namespace): Flags for the run management.
        id_generator (UniqueID): A unique ID number generator.

    Returns:
        list: All the diagrams of the appropriate Class and order.

    """
    if commands.theory == "MBPT":
        diagrams = adg.mbpt.diagrams_generation(commands.order)
    elif commands.theory in ("BMBPT", "PBMBPT"):
        diagrams = adg.bmbpt.diagrams_generation(commands.order,
                                                 commands.with_3NF,
                                                 commands.nbody_observable,
                                                 commands.canonical)
    elif commands.theory == "BIMSRG":
        diagrams, switch_flag = adg.bimsrg.diagrams_generation(commands.order)
    else:
        print("Invalid theory! Exiting program.")
        exit()
    print("Number of matrices produced: ", len(diagrams))

    diags = [nx.from_numpy_array(diagram, create_using=nx.MultiDiGraph(),
                                 parallel_edges=True) for diagram in diagrams]

    if commands.theory == "MBPT":
        for i_diag, diag in reversed_enumerate(diags):
            if (nx.number_weakly_connected_components(diag)) != 1:
                del diags[i_diag]

    adg.diag.label_vertices(diags,
                            commands.theory,
                            switch_flag if commands.theory == 'BIMSRG' else -1)

    if commands.theory in ('BMBPT', "PBMBPT"):
        diagrams = [adg.bmbpt.BmbptFeynmanDiagram(graph, id_generator.get())
                    for graph in diags]
    elif commands.theory == 'MBPT':
        diagrams = [adg.mbpt.MbptDiagram(graph, id_generator.get())
                    for graph in diags]
    elif commands.theory == "BIMSRG":
        diagrams = [adg.bimsrg.BimsrgDiagram(graph, id_generator.get())
                    for graph in diags]

    if commands.theory == "PBMBPT":
        for idx, diagram in reversed_enumerate(diagrams):
            new_graphs = adg.pbmbpt.generate_anomalous_diags(
                diagram,
                3 if commands.with_3NF else 2
            )
            new_diags = [adg.pbmbpt.ProjectedBmbptDiagram(diag,
                                                          id_generator.get(),
                                                          idx,
                                                          spawn_idx)
                         for spawn_idx, diag in enumerate(new_graphs)]
            del diagrams[idx]
            adg.pbmbpt.filter_new_diagrams(new_diags, diagrams)
            diagrams += new_diags
    return diagrams


def order_diagrams(diagrams, commands):
    """Return the ordered unique diagrams with a dict of numbers per type.

    Args:
        diagrams (list): The diagrams of the appropriate Class.
        commands (Namespace): Flags for the run management.

    Returns:
        tuple: First element is the list of ordered and unique diagrams.
        Second element is a dict with the number of diagrams per type. Third
        element is a dict with the identifiers of diagrams starting each output
        file section.

    """
    if commands.theory in ("BMBPT", "PBMBPT"):
        diagrams, diag_nbs, section_flags = adg.bmbpt.order_diagrams(diagrams)
    elif commands.theory == "MBPT":
        diagrams, diag_nbs, section_flags = adg.mbpt.order_diagrams(diagrams)
    elif commands.theory == "BIMSRG":
        order = max(commands.order)
        diagrams, diag_nbs, section_flags = adg.bimsrg.order_diagrams(diagrams,
                                                                      order)

    # Reattribute a number to the BMBPT diagrams
    if commands.theory == "BMBPT":
        for ind, diagram in enumerate(diagrams):
            diagram.tags[0] = ind

    return diagrams, diag_nbs, section_flags


def print_diags_numbers(commands, diags_nbs):
    """Print the number of diagrams for each major type.

    Args:
        commands (Namespace): Flags for the run management.
        diags_nbs (dict): The number of diagrams for each major type.

    """
    print("Number of connected diagrams: ", diags_nbs['nb_diags'])

    if commands.theory in ("BMBPT", "PBMBPT"):
        print(
            "\n2N valid diagrams: %i\n" % diags_nbs['nb_2']
            + "2N energy canonical diagrams: %i\n" % diags_nbs['nb_2_hf']
            + "2N canonical diagrams for a generic operator only: %i"
            % diags_nbs['nb_2_ehf']
        )
        if not commands.canonical:
            print(
                "2N non-canonical diagrams: %i" % diags_nbs['nb_2_not_hf']
            )
        if commands.with_3NF:
            print(
                "\n3N valid diagrams: %i\n" % diags_nbs['nb_3']
                + "3N energy canonical diagrams: %i\n" % diags_nbs['nb_3_hf']
                + "3N canonical diagrams for a generic operator only: %i"
                % diags_nbs['nb_3_ehf']
            )
            if not commands.canonical:
                print(
                    "3N non-canonical diagrams: %i" % diags_nbs['nb_3_not_hf']
                )
    elif commands.theory == "MBPT":
        print(
            "\nValid diagrams: %i\n\n" % diags_nbs['nb_diags']
            + "Singles: %i\n" % diags_nbs['singles']
            + "Doubles: %i\n" % diags_nbs['doubles']
            + "Triples: %i\n" % diags_nbs['triples']
            + "Quadruples: %i\n" % diags_nbs['quadruples']
            + "Quintuples and higher: %i" % diags_nbs['quintuples+']
        )
    elif commands.theory == "BIMSRG":
        print("\nValid diagrams: %i" % diags_nbs['nb_diags'])
        for n in range(1, commands.order[-1] + 1):
            print("d_max = %i diagrams: %i" % (n, diags_nbs[n]))
    print()


def prepare_drawing_instructions(directory, commands, diagrams, diagrams_time):
    """Write FeynMP files for the different diagrams.

    Args:
        directory (str): Path to the output folder.
        commands (Namespace): Flags for the run management.
        diagrams (list): All the diagrams of interest.
        diagrams_time (list): All the associated TSDs if appropriate.

    """
    create_feynmanmp_files(diagrams, commands.theory, directory, 'diag')
    if commands.draw_tsds:
        create_feynmanmp_files(diagrams_time, commands.theory,
                               directory, 'time')


def create_feynmanmp_files(diagrams, theory, directory, diag_type):
    """Create and move the appropriate feynmanmp files to the right place.

    Args:
        diagrams (list): The studied diagrams.
        theory (str): Name of the theory of interest.
        directory (str): Path to the result folder.
        diag_type (str): Type of studied diagrams used for drawing.

    """
    for diag in diagrams:
        if theory == "PBMBPT" and diag_type == 'diag':
            diag_name = '%s_%i_%i' % (diag_type, diag.tags[0], diag.tags[1])
        else:
            diag_name = '%s_%i' % (diag_type, diag.tags[0])
        adg.diag.feynmf_generator(diag.graph,
                                  'MBPT' if diag_type == 'time' else theory,
                                  diag_name)
        shutil.move('%s.tex' % diag_name,
                    "%s/Diagrams/%s.tex" % (directory, diag_name))


def write_file_header(latex_file, commands, diags_nbs):
    """Write the header of the result tex file.

    Args:
        latex_file (file): LaTeX output file of the program.
        commands (Namespace): Flags to manage the program's run.
        diags_nbs (dict): Number of diagrams per major type.

    """
    header = "\\documentclass[10pt,a4paper]{article}\n" \
        + "\\usepackage[left=2cm,right=2cm," \
        + "bottom=2.5cm,top=2.5cm]{geometry}\n" \
        + "\\usepackage[utf8]{inputenc}\n" \
        + "\\usepackage[T1]{fontenc}\n" \
        + "\\usepackage{lmodern}\n" \
        + "\\usepackage[hyperindex=true]{hyperref}" \
        + "\\usepackage[english]{babel}\n\\usepackage{amsmath}\n" \
        + "\\usepackage{amsfonts}\n\\usepackage{amssymb}\n"
    if commands.draw_diags:
        header = "%s\\usepackage{feynmp-auto}\n" % header
    if commands.theory in ('BMBPT', 'PBMBPT') and commands.order >= 3:
        header = "%s\\usepackage[landscape]{geometry}\n" % header

    if commands.theory != 'BIMSRG':
        header = header \
            + "\\title{Diagrams and algebraic expressions at order %i in %s}\n" \
            % (commands.order, commands.theory) \
            + "\\author{The ADG Dev Team}\n"
    else:
        header = header \
            + "\\allowdisplaybreaks\n\n" \
            + "\\title{Diagrams and algebraic expressions at order (%i,%i;%i) in BIMSRG}\n" \
            % commands.order \
            + "\\author{The ADG Dev Team}\n"
    latex_file.write("%s\n\\begin{document}\n\n\\maketitle\n\n" % header)

    if commands.theory in ("BMBPT", "PBMBPT"):
        adg.bmbpt.write_header(latex_file, commands, diags_nbs)
    elif commands.theory == "MBPT":
        adg.mbpt.write_header(latex_file, diags_nbs)
    elif commands.theory == "BIMSRG":
        adg.bimsrg.write_header(latex_file, commands, diags_nbs)

    latex_file.write("\n\\tableofcontents\n\n")

    if commands.theory == "BIMSRG":
        adg.bimsrg.write_permutator_section(latex_file, commands)


def compile_manager(directory):
    """Compile the program's LaTeX ouput file.

    Args:
        directory (str): Path to the ouput folder.

    """
    os.chdir(directory)
    os.system("pdflatex -shell-escape -interaction=batchmode result.tex")
    # Second compilation for table of contents and diagrams
    os.system("pdflatex -shell-escape -interaction=batchmode result.tex")
    os.chdir("../..")


def clean_folders(directory, commands):
    """Delete temporary files and folders.

    Args:
        directory (str): Path to the output folder.
        commands (Namespace): Flags to manage the program's run.

    """
    os.chdir(directory)
    if commands.draw_diags:
        # Remove the temporary Diagrams folder
        shutil.rmtree('Diagrams')
        # Get rid of undesired feynmp files to keep a clean directory
        for filename in os.listdir('.'):
            if filename.startswith("time") or filename.startswith("diag") \
                    or filename.startswith("equivalent"):
                os.unlink(filename)

    if commands.compile:
        # Remove LaTeX auxiliary files
        os.unlink("result.aux")
        os.unlink("result.log")
        os.unlink("result.out")
        os.unlink("result.toc")
    os.chdir("../..")
