"""Routines handling the run of ADG."""

import os
import argparse
import shutil
import networkx as nx
import adg.mbpt
import adg.bmbpt
import adg.generic_diag as gen


def parse_command_line():
    """Return a Namespace with the appropriate commands for the run."""
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
        "-o", "--order", type=int, choices=range(2, 10),
        help="order of the diagrams (>=2)")
    basic_args.add_argument(
        "-t", "--theory", type=str, choices=['MBPT', 'BMBPT'],
        help="theory of interest: MBPT or BMBPT")
    basic_args.add_argument(
        "-i", "--interactive", action="store_true",
        help="execute ADG in interactive mode")

    bmbpt_args.add_argument(
        "-n", "--norm", action="store_true",
        help="study norm BMBPT diagrams instead of operator ones")
    bmbpt_args.add_argument(
        "-3N", "--with_three_body", action="store_true",
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
        help="produce output for C. Drischler's framework (MBPT)")

    args = parser.parse_args()

    if (not args.interactive) and ((args.order is None)
                                   or (args.theory is None)):
        print "\nPlease either run the interactive mode, or the batch mode by"
        print "providing the theory and the order for the desired diagrams.\n"
        print "Use 'python2 adg.py -h' for help.\n"
        exit()

    # Avoid conflicting flags
    if args.theory != 'BMBPT' and not args.interactive:
        args.with_three_body = None
        args.norm = None
        args.draw_tsds = None
    if args.theory != 'MBPT' and not args.interactive:
        args.cd_output = None

    return args


def interactive_interface(commands):
    """Run the interactive interface mode, return the appropriate commands."""
    commands.order = int(raw_input('Order of the diagrams?\n'))
    while commands.order < 2:
        print "Perturbative order too small!"
        commands.order = int(raw_input('Order of the diagrams?\n'))
    commands.theory = raw_input('MBPT or BMBPT?\n').upper()

    if commands.theory == "BMBPT":
        commands.with_three_body = raw_input(
            "Include three-body forces? (y/N)").lower() == 'y'
        commands.norm = raw_input(
            "Compute norm diagrams, not operator ones? (y/N)").lower() == 'y'
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
    """Create missing directories and return the working directory."""
    directory = '%s/Order-%i' % (commands.theory, commands.order)
    if commands.with_three_body:
        directory += 'with3N'
    if commands.norm:
        directory += '_run_commands.norm'
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(directory+"/Diagrams"):
        os.makedirs(directory+"/Diagrams")
    return directory


def generate_diagrams(commands):
    """Return a list with diagrams of the appropriate type."""
    if commands.theory == "MBPT":
        diagrams = adg.mbpt.diagrams_generation(commands.order)
    elif commands.theory == "BMBPT":
        diagrams = adg.bmbpt.diagrams_generation(commands.order,
                                                 commands.with_three_body)
    else:
        print "Invalid theory!"
    print "Number of possible diagrams, ", len(diagrams)

    G = [nx.from_numpy_matrix(diagram, create_using=nx.MultiDiGraph(),
                              parallel_edges=True) for diagram in diagrams]

    for i_diag in xrange(len(G)-1, -1, -1):
        if (nx.number_weakly_connected_components(G[i_diag])) != 1:
            del G[i_diag]

    # Specific check for loop diagrams in BMBPT
    if commands.theory == "BMBPT":
        for i_diag in xrange(len(G)-1, -1, -1):
            if not nx.is_directed_acyclic_graph(G[i_diag]):
                del G[i_diag]

    gen.label_vertices(G, commands.theory, commands.norm)

    if commands.theory == 'BMBPT':
        diagrams = [adg.bmbpt.BmbptFeynmanDiagram(graph, commands.norm, ind)
                    for ind, graph in enumerate(G)]
    elif commands.theory == 'MBPT':
        diagrams = [adg.mbpt.MbptDiagram(graph, ind)
                    for ind, graph in enumerate(G)]
    return diagrams


def print_diags_numbers(commands, diags_nbs):
    """Print the number of diagrams for each major type."""
    print "Number of connected diagrams, ", diags_nbs['nb_diags']

    if commands.theory == "BMBPT":
        print "\n2N valid diagrams: %i" % diags_nbs['nb_2']
        print "2N energy canonical diagrams: %i" % diags_nbs['nb_2_hf']
        if not commands.norm:
            print "2N canonical diagrams for a generic operator only: %i" \
                % diags_nbs['nb_2_ehf']
        print "2N non-canonical diagrams: %i\n" % diags_nbs['nb_2_not_hf']
        if commands.with_three_body:
            print "3N valid diagrams: %i" % diags_nbs['nb_3']
            print "3N energy canonical diagrams: %i" % diags_nbs['nb_3_hf']
            if not commands.norm:
                print "3N canonical diagrams for a generic operator only: %i" \
                    % diags_nbs['nb_3_ehf']
            print "3N non-canonical diagrams: %i" % diags_nbs['nb_3_not_hf']
    elif commands.theory == "MBPT":
        print "\nValid diagrams: %i\n" % diags_nbs['nb_diags']
        print "Singles: %i" % diags_nbs['singles']
        print "Doubles: %i" % diags_nbs['doubles']
        print "Triples: %i" % diags_nbs['triples']
        print "Quadruples: %i" % diags_nbs['quadruples']
        print "Quintuples and higher: %i" % diags_nbs['quintuples+']


def prepare_drawing_instructions(directory, commands, diagrams, diagrams_time):
    """Write FeynMP files for the different diagrams."""
    create_feynmanmp_files(diagrams, commands.theory, directory, 'diag')
    if commands.draw_tsds:
        create_feynmanmp_files(diagrams_time, commands.theory,
                               directory, 'time')


def create_feynmanmp_files(diagrams, theory, directory, diag_type):
    """Create and move the appropriate feynmanmp files to the right place."""
    for diag in diagrams:
        diag_name = '%s_%i' % (diag_type, diag.tags[0])
        gen.feynmf_generator(diag.graph,
                             'MBPT' if diag_type == 'time' else theory,
                             diag_name)
        shutil.move('%s.tex' % diag_name,
                    "%s/Diagrams/%s.tex" % (directory, diag_name))


def write_file_header(latex_file, commands, diags_nbs):
    """Write the header of the result tex file."""
    header = "\\documentclass[10pt,a4paper]{article}\n" \
        + "\\usepackage[utf8]{inputenc}\n" \
        + "\\usepackage[T1]{fontenc}\n" \
        + "\\usepackage{lmodern}\n" \
        + "\\usepackage[hyperindex=true]{hyperref}" \
        + "\\usepackage[english]{babel}\n\\usepackage{amsmath}\n" \
        + "\\usepackage{amsfonts}\n\\usepackage{amssymb}\n"
    if commands.draw_diags:
        header = "%s\\usepackage{feynmp-auto}\n" % header
    if commands.order > 3:
        header = "%s\\usepackage[landscape]{geometry}\n" % header

    header = header \
        + "\\title{Diagrams and algebraic expressions at order %i in %s}\n" \
        % (commands.order, commands.theory) \
        + "\\author{RDL, JR, PA, MD, AT, TD, JPE}\n"
    latex_file.write("%s\n\\begin{document}\n\n\\maketitle\n\n" % header)

    if commands.theory == "BMBPT":
        adg.bmbpt.write_header(latex_file, commands.with_three_body,
                               commands.norm, diags_nbs)
    elif commands.theory == "MBPT":
        adg.mbpt.write_header(latex_file, diags_nbs)

    latex_file.write("\n\\tableofcontents\n\n")


def compile_and_clean(directory, pdiag):
    """Compile result.pdf and delete useless files."""
    os.chdir(directory)
    os.system("pdflatex -shell-escape -interaction=batchmode result.tex")
    if pdiag:
        # Second compilation needed
        os.system("pdflatex -shell-escape -interaction=batchmode result.tex")
        # Get rid of undesired feynmp files to keep a clean directory
        for filename in os.listdir('.'):
            if filename.startswith("time") or filename.startswith("diag") \
                    or filename.startswith("equivalent"):
                os.unlink(filename)
    print "Result saved in %s/result.pdf" % directory
