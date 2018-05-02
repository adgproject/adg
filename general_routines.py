"""Module containg methods to be called by ADG."""

import os
import shutil
import argparse
import networkx as nx


def parse_command_line():
    """Return a Namespace with the appropriate commands for the run."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="AUTOMATIC DIAGRAM GENERATOR\n\n"
        + "Generates diagrams at a given order for a "
        + "range of many-body formalisms")
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
    import mbpt
    import bmbpt
    if commands.theory == "MBPT":
        diagrams = mbpt.diagrams_generation(commands.order)
    elif commands.theory == "BMBPT":
        diagrams = bmbpt.diagrams_generation(commands.order,
                                             commands.with_three_body)
    else:
        print "Invalid theory!"
    print "Number of possible diagrams, ", len(diagrams)
    return diagrams


def no_trace(matrices):
    """Select matrices with full 0 diagonal."""
    traceless_matrices = []
    for matrix in matrices:
        test_traceless = True
        for ind_i, line in enumerate(matrix):
            if line[ind_i] != 0:
                test_traceless = False
                break
        if test_traceless:
            traceless_matrices.append(matrix)
    return traceless_matrices


def no_loop(matrices):
    """Select out matrices with loops between two vertices."""
    for i_mat in xrange(len(matrices)-1, -1, -1):
        test_no_loop = True
        matrix = matrices[i_mat]
        for ind_i in xrange(len(matrix[0])):
            for ind_j in xrange(ind_i+1):
                if (matrix[ind_i][ind_j] != 0) and (matrix[ind_j][ind_i] != 0):
                    test_no_loop = False
                    break
        if not test_no_loop:
            del matrices[i_mat]


def check_degree(matrices, three_body_use):
    """Discard matrices with wrong N-body character."""
    for i_mat in xrange(len(matrices)-1, -1, -1):
        matrix = matrices[i_mat]
        for ind_i in xrange(len(matrix[0])):
            degree = sum(matrix[ind_i][ind_j]
                         + matrix[ind_j][ind_i]
                         for ind_j in xrange(len(matrix[0])))
            if (degree != 2) and (degree != 4):
                if (not three_body_use) or (degree != 6):
                    del matrices[i_mat]
                    break


def check_vertex_degree(matrices, three_body_use, vertex_id):
    """Check the degree of a specific vertex in a set of matrices."""
    for i_mat in xrange(len(matrices)-1, -1, -1):
        matrix = matrices[i_mat]
        vertex_degree = sum(matrix[index][vertex_id] + matrix[vertex_id][index]
                            for index in xrange(len(matrix[0])))
        if (vertex_degree != 2) and (vertex_degree != 4):
            if (not three_body_use) or (vertex_degree != 6):
                del matrices[i_mat]


def topologically_distinct_diagrams(diagrams):
    """Return a list of diagrams all topologically distinct."""
    import time_structure as tsd
    iso = nx.algorithms.isomorphism
    op_nm = iso.categorical_node_match('operator', False)
    for i_diag in xrange(len(diagrams)-1, -1, -1):
        graph = diagrams[i_diag].graph
        diag_io_degrees = diagrams[i_diag].io_degrees
        for i_comp_diag in xrange(i_diag+1, len(diagrams), 1):
            if diag_io_degrees == diagrams[i_comp_diag].io_degrees:
                matcher = iso.DiGraphMatcher(graph,
                                             diagrams[i_comp_diag].graph,
                                             node_match=op_nm)
                if matcher.is_isomorphic():
                    diagrams[i_diag].tags += diagrams[i_comp_diag].tags
                    if isinstance(diagrams[i_diag], tsd.TimeStructureDiagram):
                        diagrams[i_diag].perms.update(
                            diagrams[i_comp_diag].perms)
                        diagrams[i_diag].perms[diagrams[i_comp_diag].tags[0]] \
                            = matcher.mapping
                    del diagrams[i_comp_diag]
                    break
    return diagrams


def label_vertices(graphs_list, theory_type, study_norm):
    """Account for different status of vertices in operator diagrams."""
    for graph in graphs_list:
        for node in graph:
            graph.node[node]['operator'] = False
        if (theory_type == "BMBPT") and not study_norm:
            graph.node[0]['operator'] = True


def feynmf_generator(graph, theory_type, diagram_name):
    """Generate the feynmanmp instructions corresponding to the diagram."""
    p_order = graph.number_of_nodes()
    diag_size = 20*p_order

    theories = ["MBPT", "BMBPT", "SCGF"]
    prop_types = ["half_prop", "prop_pm", "double_arrow"]
    propa = prop_types[theories.index(theory_type)]

    fmf_file = open(diagram_name + ".tex", 'w')
    fmf_file.write("\\parbox{%ipt}{\\begin{fmffile}{%s}\n" % (diag_size,
                                                              diagram_name)
                   + "\\begin{fmfgraph*}(%i,%i)\n" % (diag_size, diag_size))

    # Set the position of the vertices
    fmf_file.write("\\fmftop{v%i}\\fmfbottom{v0}\n" % (p_order-1))
    for vert in xrange(p_order-1):
        fmf_file.write("\\fmf{phantom}{v%i,v%i}\n" % (vert, (vert+1)))
        fmf_file.write(
            "\\fmfv{d.shape=square,d.filled=full,d.size=3thick"
            if graph.node[vert]['operator']
            else "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick")
        fmf_file.write("}{v%i}\n" % vert)
    fmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n"
                   % (p_order-1))
    fmf_file.write("\\fmffreeze\n")

    # Loop over all elements of the graph to draw associated propagators
    for vert_i in graph:
        for vert_j in graph:
            props_left_to_draw = graph.number_of_edges(vert_i, vert_j)
            # Special config for consecutive vertices
            if (props_left_to_draw % 2 == 1) and (abs(vert_i-vert_j) == 1):
                fmf_file.write("\\fmf{%s" % propa)
                # Check for specific MBPT configuration
                if graph.number_of_edges(vert_j, vert_i) == 1:
                    fmf_file.write(",right=0.5")
                fmf_file.write("}{v%i,v%i}\n" % (vert_i, vert_j))
                props_left_to_draw -= 1
            while props_left_to_draw > 0:
                fmf_file.write("\\fmf{%s," % propa)
                fmf_file.write("right=" if props_left_to_draw % 2 == 1
                               else "left=")
                if (props_left_to_draw == 6) or (props_left_to_draw == 5):
                    fmf_file.write("0.9")
                elif (props_left_to_draw == 4) or (props_left_to_draw == 3) \
                    or ((props_left_to_draw == 1) and
                        (graph.number_of_edges(vert_j, vert_i) == 2)):
                    fmf_file.write("0.75")
                elif (props_left_to_draw == 2) or (props_left_to_draw == 1):
                    fmf_file.write("0.5" if abs(vert_i-vert_j) == 1 else "0.6")
                fmf_file.write("}{v%i,v%i}\n" % (vert_i, vert_j))
                props_left_to_draw -= 1
    fmf_file.write("\\end{fmfgraph*}\n\\end{fmffile}}\n")
    fmf_file.close()


def create_feynmanmp_files(diagrams_list, theory_type, directory, diag_type):
    """Create and move the appropriate feynmanmp files to the right place."""
    for diag in diagrams_list:
        diag_name = '%s_%i' % (diag_type, diag.tags[0])
        feynmf_generator(diag.graph,
                         'MBPT' if diag_type == 'time' else theory_type,
                         diag_name)
        shutil.move('%s.tex' % diag_name,
                    "%s/Diagrams/%s.tex" % (directory, diag_name))


def write_file_header(latex_file, pdiag, norder, theory):
    """Write the header of the result tex file."""
    header = "\\documentclass[10pt,a4paper]{article}\n" \
        + "\\usepackage[utf8]{inputenc}\n" \
        + "\\usepackage[hyperindex=true]{hyperref}" \
        + "\\usepackage{braket}\n\\usepackage{graphicx}\n" \
        + "\\usepackage[english]{babel}\n\\usepackage{amsmath}\n" \
        + "\\usepackage{amsfonts}\n\\usepackage{amssymb}\n"
    if pdiag:
        header = "%s\\usepackage{feynmp-auto}\n" % header
    if norder > 3:
        header = "%s\\usepackage[landscape]{geometry}\n" % header

    header = header \
        + "\\title{Diagrams and algebraic expressions at order %i" % norder \
        + " in %s}\n" % theory \
        + "\\author{RDL, JR, PA, MD, AT, TD, JPE}\n"
    latex_file.write("%s\\begin{document}\n\\maketitle\n" % header)


def draw_diagram(directory, result_file, diagram_index, diag_type):
    """Copy the diagram feynmanmp instructions in the result file."""
    diag_file = open(directory+"/Diagrams/%s_%i.tex" % (diag_type,
                                                        diagram_index))
    result_file.write(diag_file.read())
    diag_file.close()


def compile_and_clean(directory, pdiag):
    """Compile result.pdf and delete useless files."""
    os.chdir(directory)
    os.system("pdflatex -shell-escape result.tex")
    if pdiag:
        # Second compilation needed
        os.system("pdflatex -shell-escape result.tex")
        # Get rid of undesired feynmp files to keep a clean directory
        for filename in os.listdir('.'):
            if filename.startswith("time") or filename.startswith("diag") \
                    or filename.startswith("equivalent"):
                os.unlink(filename)
    print "Result saved in %s/result.pdf" % directory


def to_skeleton(graph):
    """Return the bare skeleton of a graph, i.e. only non-redundant links."""
    for vertex_a in graph:
        for vertex_b in graph:
            while graph.number_of_edges(vertex_a, vertex_b) > 1:
                graph.remove_edge(vertex_a, vertex_b)
            if len(list(nx.all_simple_paths(graph, vertex_a, vertex_b))) > 1:
                while len(nx.shortest_path(graph, vertex_a, vertex_b)) == 2:
                    graph.remove_edge(vertex_a, vertex_b)
    return graph


def extract_denom(start_graph, subgraph):
    """Extract the appropriate denominator using the subgraph rule."""
    denomin = r"\epsilon^{" \
        + "".join("%s"
                  % start_graph.adj[propa[0]][propa[1]][propa[2]]['qp_state']
                  for propa
                  in start_graph.in_edges(subgraph, keys=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2])) \
        + "}_{" \
        + "".join("%s"
                  % start_graph.adj[propa[0]][propa[1]][propa[2]]['qp_state']
                  for propa
                  in start_graph.out_edges(subgraph, keys=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2])) \
        + "}"
    return denomin


class Diagram(object):
    """Describes a diagram with its related properties."""

    def __init__(self, nx_graph):
        """Generate a Diagram object starting from the NetworkX graph."""
        self.graph = nx_graph
        self.degrees = sorted([nx_graph.degree(node) for node in nx_graph])
        self.unsort_io_degrees = tuple((nx_graph.in_degree(node),
                                        nx_graph.out_degree(node))
                                       for node in nx_graph)
        self.io_degrees = sorted(self.unsort_io_degrees)
        self.max_degree = self.degrees[-1]
        self.tags = [0]

    def write_graph(self, latex_file, directory, write_time):
        """Write the graph of the diagram to the LaTeX file."""
        latex_file.write('\n\\begin{center}\n')
        draw_diagram(directory, latex_file, self.tags[0], 'diag')
        latex_file.write('\n\\end{center}\n\n')
