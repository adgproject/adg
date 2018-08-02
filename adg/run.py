"""Routines handling the run of ADG."""

import os
import argparse
import shutil
import networkx as nx
import adg.mbpt
import adg.bmbpt
import adg.diag


def parse_command_line():
    """Return run commands from the Command Line Interface.

    Returns:
        (Namespace): Appropriate commands to manage the program's run.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="AUTOMATIC DIAGRAM GENERATOR\n\n"
        + "Generates diagrams at a given order for a "
        + "range of many-body formalisms",
        epilog="If a theory is chosen in batch mode, all flags associated\n"
        + "will automatically be deactivated to avoid conflicts.\n\n")

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
        description="Arguments available only for BMBPT calculations.\n")

    basic_args.add_argument(
        "-o", "--order", type=int, choices=range(1, 10),
        help="order of the diagrams (>=1)")
    basic_args.add_argument(
        "-t", "--theory", type=str, choices=['MBPT', 'BMBPT'],
        help="theory of interest: MBPT or BMBPT")
    basic_args.add_argument(
        "-i", "--interactive", action="store_true",
        help="execute ADG in interactive mode")

    bmbpt_args.add_argument(
        "-can", "--canonical", action="store_true",
        help="consider only canonical diagrams")
    bmbpt_args.add_argument(
        "-nobs", "--nbody_observable", type=int, choices=range(1, 4), default=2,
        help="maximal n-body character of the observable [1-3], default = 2")
    bmbpt_args.add_argument(
        "-3NF", "--with_3NF", action="store_true",
        help="use two and three-body forces for BMBPT diagrams")
    bmbpt_args.add_argument(
        "-dt", "--draw_tsds", action="store_true",
        help="draw Time-Structure Diagrams (BMBPT)")

    run_args.add_argument(
        "-d", "--draw_diags", action="store_true",
        help="draw the diagrams using FeynMF")
    run_args.add_argument(
        "-c", "--compile", action="store_true",
        help="compile the LaTeX output file with PDFLaTeX")

    mbpt_args.add_argument(
        "-cd", "--cd_output", action="store_true",
        help="produce computer-readable output for automated framework (MBPT)")

    args = parser.parse_args()

    if (not args.interactive) and ((args.order is None)
                                   or (args.theory is None)):
        print "\nPlease either run the interactive mode, or the batch mode by"
        print "providing the theory and the order for the desired diagrams.\n"
        print "Use 'adg -h' for help.\n"
        exit()

    # Avoid conflicting flags
    if args.theory != 'BMBPT' and not args.interactive:
        args.canonical = None
        args.with_3NF = None
        args.nobs = 2
        args.draw_tsds = None
    if args.theory != 'MBPT' and not args.interactive:
        args.cd_output = None

    return args


def interactive_interface(commands):
    """Run the interactive interface mode, return the appropriate commands.

    Args:
        commands (Namespace): Flags for the run management.

    Returns:
        (Namespace): Flags initialized through keyboard input.

    """
    try:
        commands.order = int(raw_input('Order of the diagrams? [1-9]\n'))
    except ValueError:
        print "Please enter an integer value! Program exiting."
        exit()
    while commands.order < 1 or commands.order > 9:
        print "Perturbative order too small or too high!"
        commands.order = int(raw_input('Order of the diagrams? [1-9]\n'))

    theories = ["BMBPT", "MBPT"]

    commands.theory = raw_input('MBPT or BMBPT?\n').upper()
    while commands.theory not in theories:
        print "Invalid theory!"
        commands.theory = raw_input('MBPT or BMBPT?\n').upper()

    if commands.theory == "BMBPT":
        commands.canonical = raw_input(
            "Consider only canonical diagrams? (y/N)").lower() == 'y'
        commands.with_3NF = raw_input(
            "Include three-body forces? (y/N)").lower() == 'y'
        try:
            commands.nbody_observable = int(raw_input(
                "Maximal n-body character of the observable? [1-3]"))
        except ValueError:
            print "Please enter an integer value! Program exiting."
            exit()
        while commands.nbody_observable < 1 or commands.nbody_observable > 3:
            print "The observable must be 1-body, 2-body or 3-body!"
            commands.nbody_observable = int(raw_input(
                "Maximal n-body character of the observable? [1-3]"))
        commands.draw_tsds = raw_input(
            "Draw time-structure diagrams? (y/N)").lower() == 'y'

    commands.draw_diags = raw_input(
        "Generate diagrams FeynMF instructions in TeX file? (y/N) "
        ).lower() == 'y'

    if commands.theory == "MBPT":
        commands.cd_output = raw_input(
            "Produce a CD output file? (y/N) ").lower() == 'y'
    commands.compile = raw_input("Compile pdf? (y/N) ").lower() == 'y'

    return commands


def attribute_directory(commands):
    """Create missing directories and return the working directory.

    Args:
        commands (Namespace): Flags for the run management.

    Returns:
        (str): Path to the result folder.

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

    """
    directory = '%s/Order-%i' % (commands.theory, commands.order)
    if commands.canonical:
        directory += '_canonical'
    if commands.theory == 'BMBPT':
        directory += '_%ibody_observable' % commands.nbody_observable
    if commands.with_3NF:
        directory += '_with3N'
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(directory+"/Diagrams"):
        os.makedirs(directory+"/Diagrams")
    return directory


def generate_diagrams(commands):
    """Return a list with diagrams of the appropriate type.

    Args:
        commands (Namespace): Flags for the run management.

    Returns:
        (list): All the diagrams of the appropriate Class and order.

    """
    if commands.theory == "MBPT":
        diagrams = adg.mbpt.diagrams_generation(commands.order)
    elif commands.theory == "BMBPT":
        diagrams = adg.bmbpt.diagrams_generation(commands.order,
                                                 commands.with_3NF,
                                                 commands.nbody_observable,
                                                 commands.canonical)
    else:
        print "Invalid theory! Exiting program."
        exit()
    print "Number of matrices produced: ", len(diagrams)

    diags = [nx.from_numpy_matrix(diagram, create_using=nx.MultiDiGraph(),
                                  parallel_edges=True) for diagram in diagrams]

    for i_diag in xrange(len(diags)-1, -1, -1):
        if (nx.number_weakly_connected_components(diags[i_diag])) != 1:
            del diags[i_diag]

    adg.diag.label_vertices(diags, commands.theory)

    if commands.theory == 'BMBPT':
        diagrams = [adg.bmbpt.BmbptFeynmanDiagram(graph, ind)
                    for ind, graph in enumerate(diags)]
    elif commands.theory == 'MBPT':
        diagrams = [adg.mbpt.MbptDiagram(graph, ind)
                    for ind, graph in enumerate(diags)]
    return diagrams


def order_diagrams(diagrams, commands):
    """Return the ordered unique diagrams with a dict of numbers per type.

    Args:
        diagrams (list): The diagrams of the appropriate Class.
        commands (Namespace): Flags for the run management.

    Returns:
        (tuple): First element is the list of ordered and unique diagrams.
        Second element is a dict with the number of diagrams per type.

    """
    if commands.theory == "BMBPT":
        diagrams, diags_per_type = adg.bmbpt.order_diagrams(diagrams)

    elif commands.theory == "MBPT":
        diagrams, diags_per_type = adg.mbpt.order_diagrams(diagrams)

    return diagrams, diags_per_type


def print_diags_numbers(commands, diags_nbs):
    """Print the number of diagrams for each major type.

    Args:
        commands (Namespace): Flags for the run management.
        diags_nbs (dict): The number of diagrams for each major type.

    """
    print "Number of connected diagrams: ", diags_nbs['nb_diags']

    if commands.theory == "BMBPT":
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
    print


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
        + "\\usepackage[utf8]{inputenc}\n" \
        + "\\usepackage[T1]{fontenc}\n" \
        + "\\usepackage{lmodern}\n" \
        + "\\usepackage[hyperindex=true]{hyperref}" \
        + "\\usepackage[english]{babel}\n\\usepackage{amsmath}\n" \
        + "\\usepackage{amsfonts}\n\\usepackage{amssymb}\n"
    if commands.draw_diags:
        header = "%s\\usepackage{feynmp-auto}\n" % header
    if commands.theory == 'BMBPT' and commands.order >= 3:
        header = "%s\\usepackage[landscape]{geometry}\n" % header

    header = header \
        + "\\title{Diagrams and algebraic expressions at order %i in %s}\n" \
        % (commands.order, commands.theory) \
        + "\\author{RDL, JR, PA, MD, AT, TD, JPE}\n"
    latex_file.write("%s\n\\begin{document}\n\n\\maketitle\n\n" % header)

    if commands.theory == "BMBPT":
        adg.bmbpt.write_header(latex_file, commands, diags_nbs)
    elif commands.theory == "MBPT":
        adg.mbpt.write_header(latex_file, diags_nbs)

    latex_file.write("\n\\tableofcontents\n\n")


def compile_manager(directory, pdiag):
    """Compile the program's LaTeX ouput file.

    Args:
        directory (str): Path to the ouput folder.
        pdiag (bool): ``True`` if one wants to draw the diagrams.

    """
    os.chdir(directory)
    os.system("pdflatex -shell-escape -interaction=batchmode result.tex")
    if pdiag:
        # Second compilation needed
        os.system("pdflatex -shell-escape -interaction=batchmode result.tex")
    os.chdir("../..")
    print "Result saved in %s/result.pdf" % directory


def clean_folders(directory, commands):
    """Delete temporary files and folders.

    Args:
        directory (str): Path to the ouput folder.
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
