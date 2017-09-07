"""Module containg methods to be called by ADG."""

import os
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
    no_loop_matrices = []
    for matrix in matrices:
        test_no_loop = True
        for ind_i in range(len(matrix[0])):
            for ind_j in range(ind_i+1):
                if (matrix[ind_i][ind_j] != 0) and (matrix[ind_j][ind_i] != 0):
                    test_no_loop = False
                    break
        if test_no_loop:
            no_loop_matrices.append(matrix)
    return no_loop_matrices


def check_degree(matrices, three_N_use):
    """Check the degrees of the vertices
    (i.e. its effective one-, two- or three-body structure).
    """
    deg_ok = []
    for matrix in matrices:
        test_degree = True
        for ind_i in range(len(matrix[0])):
            degree = 0
            for ind_j in range(len(matrix[0])):
                degree += matrix[ind_i][ind_j] + matrix[ind_j][ind_i]
            if (degree != 2) and (degree != 4):
                if (not three_N_use) or (degree != 6):
                    test_degree = False
                    break
        if test_degree:
            deg_ok.append(matrix)
    return deg_ok


def check_vertex_degree(matrices, three_N_use, vertex_id):
    """Check the degree of a specific vertex in a set of matrices."""
    good_matrices = []
    for matrix in matrices:
        vertex_degree_OK = True
        vertex_degree = 0
        for index in range(len(matrix[0])):
            vertex_degree += matrix[index][vertex_id] \
                + matrix[vertex_id][index]
        if (vertex_degree != 2) and (vertex_degree != 4):
            if (not three_N_use) or (vertex_degree != 6):
                vertex_degree_OK = False
        if vertex_degree_OK:
            good_matrices.append(matrix)
    return good_matrices


def empty_matrix_generation(size):
    """Generate an empty matrix of size (size,size)."""
    empty_matrix = []
    for line in range(size):
        empty_matrix.append([])
        for element in range(size):
            empty_matrix[line].append(0)
    return empty_matrix


def topologically_distinct_diags(diagrams):
    """Return a list of diagrams all topologically distinct."""
    distinct_diagrams = []
    nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
    for diag in diagrams:
        if distinct_diagrams == []:
            distinct_diagrams.append(diag)
        else:
            test = True
            for good_diag in distinct_diagrams:
                if nx.is_isomorphic(diag, good_diag, node_match=nm):
                    test = False
                    break
            if test:
                distinct_diagrams.append(diag)
    return distinct_diagrams


def has_only_branch_operator_subgraphs(diagram):
    """Return True if diagram has operator subgraphs that are all branches."""
    has_branch_subgraphs = True
    for connected_subgraph in nx.weakly_connected_component_subgraphs(diagram):
        if len(connected_subgraph) > 1:
            if nx.dag_longest_path_length(connected_subgraph) != (len(connected_subgraph)-1):
                has_branch_subgraphs = False
    return has_branch_subgraphs


def number_of_sinks(diagram):
    """Return the number of vertices in the graph with no edges going out."""
    nb = 0
    for vertex in diagram:
        if diagram.out_degree(vertex) == 0:
            nb += 1
    return nb


def feynmf_generator(start_diag, theory_type, diagram_name):
    """Generate the feynmanmp instructions corresponding to the diagram."""
    p_order = start_diag.number_of_nodes()
    diag_size = 20*p_order

    theories = ["MBPT", "BMBPT", "SCGF"]
    prop_types = ["half_prop", "prop_pm", "double_arrow"]
    propa = prop_types[theories.index(theory_type)]

    fmf_file = open(diagram_name + ".tex", 'w')
    begin_file = "\\parbox{%i" % diag_size + "pt}{\\begin{fmffile}{" \
        + diagram_name + "}\n\\begin{fmfgraph*}(%i" % diag_size \
        + ",%i)\n" % diag_size
    fmf_file.write(begin_file)

    # Set the position of the vertices
    fmf_file.write("\\fmftop{v%i}\\fmfbottom{v0}\n" % (p_order-1))
    for vert in range(p_order-1):
        fmf_file.write("\\fmf{phantom}{v%i" % vert + ",v%i}\n" % (vert+1))
        if start_diag.node[vert]['operator']:
            fmf_file.write("\\fmfv{d.shape=square,d.filled=full,d.size=3thick")
        else:
            fmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick")
        fmf_file.write("}{v%i}\n" % vert)
    fmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}")
    fmf_file.write("{v%i}\n" % (p_order-1))
    fmf_file.write("\\fmffreeze\n")

    # Loop over all elements of the graph to draw associated propagators
    for vert_i in start_diag:
        for vert_j in start_diag:
            props_left_to_draw = start_diag.number_of_edges(vert_i, vert_j)
            # Special config for consecutive vertices
            if (props_left_to_draw % 2 == 1) and (abs(vert_i-vert_j) == 1):
                fmf_file.write("\\fmf{" + propa)
                # Check for specific MBPT configuration
                if start_diag.number_of_edges(vert_j, vert_i) == 1:
                    fmf_file.write(",right=0.5")
                fmf_file.write("}{v%i," % vert_i + "v%i}\n" % vert_j)
                props_left_to_draw -= 1
            while props_left_to_draw > 0:
                fmf_file.write("\\fmf{" + propa + ",")
                if props_left_to_draw % 2 == 1:
                    fmf_file.write("right=")
                else:
                    fmf_file.write("left=")
                if (props_left_to_draw == 6) or (props_left_to_draw == 5):
                    fmf_file.write("0.9")
                elif (props_left_to_draw == 4) or (props_left_to_draw == 3):
                    fmf_file.write("0.75")
                elif (props_left_to_draw == 2) or (props_left_to_draw == 1):
                    if abs(vert_i-vert_j) == 1:
                        fmf_file.write("0.5")
                    else:
                        fmf_file.write("0.6")
                fmf_file.write("}{v%i," % vert_i + "v%i}\n" % vert_j)
                props_left_to_draw -= 1
    fmf_file.write("\\end{fmfgraph*}\n\\end{fmffile}}\n")
    fmf_file.close()


def write_file_header(directory, latex_file, pdiag, norder, theory):
    """Write the header of the result tex file."""
    header = "\\documentclass[10pt,a4paper]{article}\n" \
        + "\\usepackage[utf8]{inputenc}\n" \
        + "\\usepackage[hyperindex=true]{hyperref}" \
        + "\\usepackage{braket}\n\\usepackage{graphicx}\n" \
        + "\\usepackage[english]{babel}\n\\usepackage{amsmath}\n" \
        + "\\usepackage{amsfonts}\n\\usepackage{amssymb}\n"
    if pdiag:
        header = header + "\\usepackage[force]{feynmp-auto}\n"
    land = False
    if norder > 3:
        msg = 'Expressions may be long, rotate pdf?'
        land = raw_input("%s (y/N) " % msg).lower() == 'y'
    if land:
        header = header + "\\usepackage[landscape]{geometry}\n"

    header = header \
        + "\\title{Diagrams and algebraic expressions at order %i" % norder \
        + " in " + theory + "}\n" \
        + "\\author{RDL, JR, PA, MD, AT, TD, JPE}\n"
    latex_file.write(header)
    begdoc = "\\begin{document}\n"
    latex_file.write(begdoc)
    latex_file.write("\\maketitle\n")
    latex_file.write("\\graphicspath{{Diagrams/}}")


def draw_diagram(directory, result_file, diagram_index, diag_type):
    """Copy the diagram feynmanmp instructions in the result file."""
    if diag_type == 'diag':
        diag_file = open(directory+"/Diagrams/diag_%i.tex" % diagram_index)
    elif diag_type == 'time':
        diag_file = open(directory+"/Diagrams/time_%i.tex" % diagram_index)
    result_file.write(diag_file.read())


def compile_and_clean(directory, pdiag, numdiag, write_time, nb_time_diags):
    """Compile result.pdf and delete useless files."""
    os.chdir(directory)
    os.system("pdflatex -shell-escape result.tex")
    if pdiag:
        # Second compilation needed
        os.system("pdflatex -shell-escape result.tex")
        # Get rid of undesired feynmp files to keep a clean directory
        for diagram in range(0, numdiag):
            os.unlink("diag_%i.1" % diagram)
            os.unlink("diag_%i.mp" % diagram)
            os.unlink("diag_%i.log" % diagram)
        if write_time:
            for i_tdiag in range(nb_time_diags):
                os.unlink("time_%i.1" % i_tdiag)
                os.unlink("time_%i.mp" % i_tdiag)
                os.unlink("time_%i.log" % i_tdiag)
    print "Result saved in "+directory + '/result.pdf'
