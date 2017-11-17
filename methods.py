"""Module containg methods to be called by ADG."""

import os
import shutil
import networkx as nx


def no_trace(matrices):
    """Select matrices with full 0 diagonal."""
    traceless_matrices = []
    for matrix in matrices:
        test_traceless = True
        for ind_i, n in enumerate(matrix):
            if n[ind_i] == 1:
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


def check_degree(matrices, three_N_use):
    """Discard matrices with wrong N-body character."""
    for i_mat in xrange(len(matrices)-1, -1, -1):
        matrix = matrices[i_mat]
        for ind_i in xrange(len(matrix[0])):
            degree = sum(matrix[ind_i][ind_j]
                         + matrix[ind_j][ind_i]
                         for ind_j in xrange(len(matrix[0])))
            if (degree != 2) and (degree != 4):
                if (not three_N_use) or (degree != 6):
                    del matrices[i_mat]
                    break


def check_vertex_degree(matrices, three_N_use, vertex_id):
    """Check the degree of a specific vertex in a set of matrices."""
    for i_mat in xrange(len(matrices)-1, -1, -1):
        matrix = matrices[i_mat]
        vertex_degree = sum(matrix[index][vertex_id] + matrix[vertex_id][index]
                            for index in xrange(len(matrix[0])))
        if (vertex_degree != 2) and (vertex_degree != 4):
            if (not three_N_use) or (vertex_degree != 6):
                del matrices[i_mat]


def topologically_distinct_graphs(graphs):
    """Return a list of graphs all topologically distinct."""
    nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
    for i_graph in xrange(len(graphs)-1, -1, -1):
        graph = graphs[i_graph]
        vert_degrees = sorted([degree for node, degree in graph.degree_iter()])
        for i_comp_graph in xrange(i_graph+1, len(graphs), 1):
            if vert_degrees == sorted([
                    degree for node, degree
                    in graphs[i_comp_graph].degree_iter()]):
                if nx.is_isomorphic(graph, graphs[i_comp_graph],
                                    node_match=nm):
                    del graphs[i_comp_graph]
                    break
    return graphs


def topologically_distinct_diagrams(diagrams):
    """Return a list of diagrams all topologically distinct."""
    nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
    for i_diag in xrange(len(diagrams)-1, -1, -1):
        graph = diagrams[i_diag].graph
        diag_io_degrees = diagrams[i_diag].io_degrees
        for i_comp_diag in xrange(i_diag+1, len(diagrams), 1):
            if diag_io_degrees == diagrams[i_comp_diag].io_degrees:
                if nx.is_isomorphic(graph, diagrams[i_comp_diag].graph,
                                    node_match=nm):
                    diagrams[i_diag].tags += diagrams[i_comp_diag].tags
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


def has_only_branch_operator_subgraphs(graph):
    """Return True if graph has operator subgraphs that are all branches."""
    has_branch_subgraphs = True
    for connected_subgraph in nx.weakly_connected_component_subgraphs(graph):
        if len(connected_subgraph) > 1:
            if nx.dag_longest_path_length(connected_subgraph) != \
                    (len(connected_subgraph)-1):
                has_branch_subgraphs = False
    return has_branch_subgraphs


def number_of_sinks(graph):
    """Return the number of vertices in the graph with no edges going out."""
    return sum(1 for vertex in graph if graph.out_degree(vertex) == 0)


def feynmf_generator(start_graph, theory_type, diagram_name):
    """Generate the feynmanmp instructions corresponding to the diagram."""
    p_order = start_graph.number_of_nodes()
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
        fmf_file.write("\\fmfv{d.shape=square,d.filled=full,d.size=3thick"
                       if start_graph.node[vert]['operator']
                       else "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick")
        fmf_file.write("}{v%i}\n" % vert)
    fmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n"
                   % (p_order-1))
    fmf_file.write("\\fmffreeze\n")

    # Loop over all elements of the graph to draw associated propagators
    for vert_i in start_graph:
        for vert_j in start_graph:
            props_left_to_draw = start_graph.number_of_edges(vert_i, vert_j)
            # Special config for consecutive vertices
            if (props_left_to_draw % 2 == 1) and (abs(vert_i-vert_j) == 1):
                fmf_file.write("\\fmf{%s" % propa)
                # Check for specific MBPT configuration
                if start_graph.number_of_edges(vert_j, vert_i) == 1:
                    fmf_file.write(",right=0.5")
                fmf_file.write("}{v%i,v%i}\n" % (vert_i, vert_j))
                props_left_to_draw -= 1
            while props_left_to_draw > 0:
                fmf_file.write("\\fmf{%s," % propa)
                fmf_file.write("right=" if props_left_to_draw % 2 == 1
                               else "left=")
                if (props_left_to_draw == 6) or (props_left_to_draw == 5):
                    fmf_file.write("0.9")
                elif (props_left_to_draw == 4) or (props_left_to_draw == 3):
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
        diag_name = '%s%i' % (diag_type, diag.tags[0])
        feynmf_generator(diag.graph,
                         'MBPT' if diag_type == 'time' else theory_type,
                         diag_name)
        shutil.move('%s.tex' % diag_name,
                    "%s/Diagrams/%s.tex" % (directory, diag_name))


def write_file_header(directory, latex_file, pdiag, norder, theory):
    """Write the header of the result tex file."""
    header = "\\documentclass[10pt,a4paper]{article}\n" \
        + "\\usepackage[utf8]{inputenc}\n" \
        + "\\usepackage[hyperindex=true]{hyperref}" \
        + "\\usepackage{braket}\n\\usepackage{graphicx}\n" \
        + "\\usepackage[english]{babel}\n\\usepackage{amsmath}\n" \
        + "\\usepackage{amsfonts}\n\\usepackage{amssymb}\n"
    if pdiag:
        header = "%s\\usepackage{feynmp-auto}\n" % header
    land = False
    if norder > 3:
        msg = 'Expressions may be long, rotate pdf?'
        land = raw_input("%s (y/N) " % msg).lower() == 'y'
    if land:
        header = "%s\\usepackage[landscape]{geometry}\n" % header

    header = header \
        + "\\title{Diagrams and algebraic expressions at order %i" % norder \
        + " in %s}\n" % theory \
        + "\\author{RDL, JR, PA, MD, AT, TD, JPE}\n"
    latex_file.write("%s\\begin{document}\n" % header)
    latex_file.write("\\maketitle\n")
    latex_file.write("\\graphicspath{{Diagrams/}}")


def draw_diagram(directory, result_file, diagram_index, diag_type):
    """Copy the diagram feynmanmp instructions in the result file."""
    diag_file = open(directory+"/Diagrams/%s%i.tex" % (diag_type,
                                                       diagram_index))
    result_file.write(diag_file.read())


def compile_and_clean(directory, pdiag, diagrams, write_time, time_diagrams):
    """Compile result.pdf and delete useless files."""
    os.chdir(directory)
    os.system("pdflatex -shell-escape result.tex")
    if pdiag:
        # Second compilation needed
        os.system("pdflatex -shell-escape result.tex")
        # Get rid of undesired feynmp files to keep a clean directory
        for diagram in diagrams:
            os.unlink("diag_%i.1" % diagram.tags[0])
            os.unlink("diag_%i.mp" % diagram.tags[0])
            os.unlink("diag_%i.log" % diagram.tags[0])
        if write_time:
            for tdiag in time_diagrams:
                os.unlink("time_%i.1" % tdiag.tags[0])
                os.unlink("time_%i.mp" % tdiag.tags[0])
                os.unlink("time_%i.log" % tdiag.tags[0])
        for filename in os.listdir('.'):
            if filename.startswith('equivalent'):
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


class Diagram(object):
    """Describes a diagram with its related properties."""

    def __init__(self, nx_graph):
        self.graph = nx_graph
        self.degrees = sorted([degree for node, degree
                               in nx_graph.degree_iter()])
        self.unsort_io_degrees = tuple((nx_graph.in_degree(node),
                                        nx_graph.out_degree(node))
                                       for node in nx_graph)
        self.io_degrees = sorted(self.unsort_io_degrees)
        self.max_degree = self.degrees[-1]
