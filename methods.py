"""Module containg methods to be called by ADG."""

import os
import copy
import itertools
import string
import numpy as np
import networkx as nx


def seed(n):
    """Generate all 1-magic square of dimension n."""
    return [k for k in itertools.permutations(range(n), n)]


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


def diagram_generation(n):
    """Generate the diagrams for the MBPT case."""
    seeds = seed(n)
    all = [[[0 if i != j else 1 for i in range(n)] for j in k] for k in seeds]
    traceless = no_trace(all)
    coeffs = [i for i in itertools.combinations_with_replacement(
        range(len(traceless)), 2)]
    double = []

    for coef in coeffs:
        matrix = copy.deepcopy(traceless[coef[0]])
        for i, line in enumerate(traceless[coef[1]]):
            for j, elem in enumerate(line):
                matrix[i][j] += elem
        double.append(matrix)
    doubleUniq = []
    for matrix in double:
        if matrix not in doubleUniq:
            doubleUniq.append(matrix)
    doubleUniq.sort(reverse=True)
    mbpt_diagrams = []
    for matrix in doubleUniq:
        mbpt_diagrams.append(np.array(matrix))
    return mbpt_diagrams


def BMBPT_generation(p_order, three_N_use, norm_diagrams):
    """Generate diagrams for BMBPT from bottom up."""
    deg_max = 4
    if three_N_use:
        deg_max = 6

    # Create a null oriented adjacency matrix of dimension (p_order,p_order)
    temp_matrices = []
    temp_matrices.append(empty_matrix_generation(p_order))

    # Generate oriented adjacency matrices going vertex-wise
    for vertex in range(p_order):
        for sum_index in range(vertex+1, p_order):
            matrices = []
            for mat in temp_matrices:
                matrices.append(mat)
                if mat[sum_index][vertex] == 0:
                    vert_degree = 0
                    for k in range(0, p_order):
                        vert_degree += mat[k][vertex] + mat[vertex][k]
                    elem = 1
                    while (elem + vert_degree) <= deg_max:
                        temp_mat = copy.deepcopy(mat)
                        temp_mat[vertex][sum_index] = elem
                        matrices.append(temp_mat)
                        elem += 1
            temp_matrices = copy.deepcopy(matrices)
        temp_matrices = check_vertex_degree(matrices, three_N_use, vertex)
        matrices = copy.deepcopy(temp_matrices)

    # Checks to exclude non-conform matrices
    matrices = check_degree(matrices, three_N_use)
    matrices = no_loop(matrices)
    matricesUniq = []
    for mat in matrices:
        if mat not in matricesUniq:
            matricesUniq.append(mat)
    matricesUniq.sort(reverse=True)
    bmbpt_diagrams = []
    for mat in matricesUniq:
        bmbpt_diagrams.append(np.array(mat))
    return bmbpt_diagrams


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


def order_2B_or_3B(diagrams, TwoB_diags, ThreeB_diags):
    """Order the diagrams depending on their 2B or 3B effective status."""
    for diag in diagrams:
        max_deg = 0
        for node in diag:
            max_deg = max(max_deg, diag.degree(node))
        if max_deg == 6:
            ThreeB_diags.append(diag)
        else:
            TwoB_diags.append(diag)


def order_HF_or_not(diagrams, HF_diags, EHF_diags, noHF_diags, norm):
    """Order the diagrams depending on their HF status."""
    for diag in diagrams:
        test_HF = True
        test_EHF = True
        for node in diag:
            if diag.degree(node) == 2:
                test_HF = False
                if node != 0:
                    test_EHF = False
        if test_HF:
            HF_diags.append(diag)
        elif (not test_EHF) or norm:
            noHF_diags.append(diag)
        else:
            EHF_diags.append(diag)


def line_label_h(n):
    """Select appropriate label for hole line."""
    labels = list(string.ascii_lowercase)
    labels = labels[0:15]
    return labels[n]


def line_label_p(n):
    """Select appropriate label for particle line."""
    labels = list(string.ascii_lowercase)
    labels = labels[15:-1]
    return labels[n]


def attribute_qp_labels(diagram):
    """Attribute the appropriate qp labels to the graph's propagators."""
    i = 1
    for prop in diagram.edges_iter(keys=True):
        diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state'] = "k_{%i}" % i
        i += 1


def time_structure_graph(diagram):
    """Return the time-structure graph associated to the diagram."""
    time_diag = diagram.to_directed()
    if time_diag.node[0]['operator']:
        for vertex in range(1, len(time_diag)):
            time_diag.add_edge(0, vertex)
    for vertex_a in time_diag:
        for vertex_b in time_diag:
            while time_diag.number_of_edges(vertex_a, vertex_b) > 1:
                time_diag.remove_edge(vertex_a, vertex_b)
            if len(list(nx.all_simple_paths(time_diag, vertex_a, vertex_b))) > 1:
                while len(nx.shortest_path(time_diag, vertex_a, vertex_b)) == 2:
                    time_diag.remove_edge(vertex_a, vertex_b)
    return time_diag


def omega_subgraph(diagram):
    """Return the graph without any operator vertex."""
    subgraph_stack = []
    for vertex in diagram:
        if diagram.node[vertex]['operator'] is False:
            subgraph_stack.append(vertex)
    return diagram.subgraph(subgraph_stack)


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


def extract_numerator(diagram):
    """"Return the numerator associated to a BMBPT diagram."""
    numerator = ""
    for vertex in diagram:
        # Attribute the correct operator to each vertex
        if diagram.node[vertex]['operator']:
            numerator += "O"
        else:
            numerator += "\\Omega"
        # Attribute the good "type number" to each vertex
        numerator = numerator + "^{%i" % diagram.out_degree(vertex) \
            + "%i}_{" % diagram.in_degree(vertex)
        # First add the qp states corresponding to propagators going out
        for prop in diagram.out_edges_iter(vertex, keys=True):
            numerator += diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state']
        # Add the qp states corresponding to propagators coming in
        previous_vertex = vertex - 1
        while previous_vertex >= 0:
            for prop in diagram.in_edges_iter(vertex, keys=True):
                if prop[0] == previous_vertex:
                    numerator += diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state']
            previous_vertex -= 1
        numerator += "} "
    return numerator


def extract_denom(start_diag, subdiagram):
    """Extract the appropriate denominator using the subdiagram rule."""
    denomin = ""
    for propa in start_diag.in_edges_iter(subdiagram, keys=True):
        if subdiagram.has_edge(propa[0], propa[1], propa[2]) is False:
            denomin += " + E_{%s}" \
                % start_diag.edge[propa[0]][propa[1]][propa[2]]['qp_state']
    for propa in start_diag.out_edges_iter(subdiagram, keys=True):
        if subdiagram.has_edge(propa[0], propa[1], propa[2]) is False:
            denomin += " - E_{%s}" \
                % start_diag.edge[propa[0]][propa[1]][propa[2]]['qp_state']
    return denomin


def time_tree_denominator(diagram, time_diagram, denominator):
    """Add the denominator for a time-tree diagram."""
    for vertex_i in range(1, len(time_diagram)):
        subgraph_stack = []
        subgraph_stack.append(diagram.nodes()[vertex_i])
        for vertex_j in nx.descendants(time_diagram, vertex_i):
            subgraph_stack.append(diagram.nodes()[vertex_j])
        subdiag = diagram.subgraph(subgraph_stack)
        denominator += "(" + extract_denom(diagram, subdiag) + ")"
    return denominator


def extract_integral(diagram):
    """Return the integral part of the Feynman expression of the diagram."""
    integral = ""
    norder = diagram.number_of_nodes()
    for vertex in range(1, norder):
        integral += "\\mathrm{d}\\tau_%i" % vertex
    if norder > 2:
        for vertex_i in range(1, norder):
            for vertex_j in range(1, norder):
                if diagram.has_edge(vertex_i, vertex_j):
                    integral += "\\theta(\\tau_%i" % vertex_j \
                        + "-\\tau_%i) " % vertex_i
    for vertex in range(1, norder):
        integral += "e^{-\\tau_%i (" % vertex
        for prop in diagram.in_edges_iter(vertex, keys=True):
            integral += " + E_{%s}" \
                % diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state']
        for prop in diagram.out_edges_iter(vertex, keys=True):
            integral += " - E_{%s}" \
                % diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state']
        integral += ")}"
    return integral


def extract_BMBPT_crossing_sign(diagram):
    """Return True if there's a sign factor associated with crossing propagators.

    Use the fact that all lines propagate upwards and the
    canonical representation of the diagrams and vertices.
    """
    nb_crossings = 0
    for vertex in diagram:
        for propagator in diagram.out_edges_iter(vertex, keys=True):
            for vertex_ante in range(propagator[0]):
                for vertex_post in range(propagator[0]+1, propagator[1]):
                    nb_crossings += diagram.number_of_edges(vertex_ante,
                                                            vertex_post)
    return nb_crossings % 2 == 1


def multiplicity_symmetry_factor(diagram):
    """Return the symmetry factor associated with propagators multiplicity."""
    factor = ""
    prop_multiplicity = []
    for i in range(6):
        prop_multiplicity.append(0)
    for vertex_i in diagram:
        for vertex_j in diagram:
            if diagram.number_of_edges(vertex_i, vertex_j) >= 2:
                prop_multiplicity[diagram.number_of_edges(
                    vertex_i, vertex_j) - 1] += 1

    for prop_id, multiplicity in enumerate(prop_multiplicity):
        if multiplicity == 1:
            factor += "(%i!)" % (prop_id+1)
        elif multiplicity >= 2:
            factor += "(%i!)" % (prop_id+1) + "^%i" % multiplicity
    return factor


def vertex_exchange_sym_factor(diagram):
    """Return the symmetry factor associated with vertex exchange."""
    # Starts at -2 as the identity belongs to the set of permutations
    factor = -2
    non_op_vertices = []
    for vertex in diagram:
        if diagram.node[vertex]['operator'] is False:
            non_op_vertices.append(vertex)
    for permutation in itertools.permutations(non_op_vertices, len(non_op_vertices)):
        mapping = dict(zip(non_op_vertices, permutation))
        permuted_diag = nx.relabel_nodes(diagram, mapping, copy=True)
        if nx.is_isomorphic(diagram, nx.intersection(diagram, permuted_diag)):
            factor += 2
    if factor != 0:
        return "%i" % factor
    else:
        return ""


def has_tree_time_structure(diagram):
    """Return True if the time structure of the diagram is a tree."""
    diag_copy = diagram.to_directed()
    for vertex in range(1, len(diag_copy)):
        if diag_copy.in_degree(vertex) == 0:
            diag_copy.add_edge(0, vertex)
    time_diag = nx.DiGraph()
    for vertex_i in diag_copy:
        for vertex_j in diag_copy:
            lgst_path = []
            for path in nx.all_simple_paths(diag_copy,
                                            source=vertex_i,
                                            target=vertex_j):
                if len(path) > len(lgst_path):
                    lgst_path = path
            time_diag.add_path(lgst_path)
    return nx.is_arborescence(time_diag)


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


def write_BMBPT_header(tex_file, numdiag, three_N, norm, nb_2_HF,
                       nb_2_EHF, nb_2_noHF, nb_3_HF, nb_3_EHF, nb_3_noHF):
    """Write overall header for BMBPT result file."""
    tex_file.write("Valid diagrams: %i\n\n" % numdiag)
    tex_file.write("2N valid diagrams: %i\n\n" %
                   (nb_2_HF + nb_2_EHF + nb_2_noHF))
    tex_file.write("2N canonical diagrams for the energy: %i\n\n" % nb_2_HF)
    if not norm:
        tex_file.write(
            "2N canonical diagrams for a generic operator only: %i\n\n" % nb_2_EHF)
    tex_file.write("2N non-canonical diagrams: %i\n\n" % nb_2_noHF)
    if three_N:
        tex_file.write("3N valid diagrams: %i\n\n" %
                       (nb_3_HF + nb_3_EHF + nb_3_noHF))
        tex_file.write(
            "3N canonical diagrams for the energy: %i\n\n" % nb_3_HF)
        if not norm:
            tex_file.write(
                "3N canonical diagrams for a generic operator only: %i\n\n" % nb_3_EHF)
        tex_file.write("3N non-canonical diagrams: %i\n\n" % nb_3_noHF)


def write_BMBPT_section(result, diag_index, three_N, norm,
                        nb_2, nb_2_HF, nb_2_EHF, nb_3_HF, nb_3_EHF):
    """Write section and subsections for BMBPT result file."""
    if (diag_index == nb_2_HF) and (not norm):
        result.write("\\subsection{Two-body canonical diagrams for a generic operator only}\n\n")
    elif diag_index == nb_2_HF + nb_2_EHF:
        result.write("\\subsection{Two-body non-canonical diagrams}\n\n")
    if three_N:
        if diag_index == nb_2:
            result.write("\\section{Three-body diagrams}\n\n")
            result.write("\\subsection{Three-body energy canonical diagrams}\n\n")
        elif (diag_index == nb_2 + nb_3_HF) and (not norm):
            result.write("\\subsection{Three-body canonical diagrams for a generic operator only}\n\n")
        elif diag_index == nb_2 + nb_3_HF + nb_3_EHF:
            result.write("\\subsection{Three-body non-canonical diagrams}\n\n")


def draw_diagram(directory, result_file, diagram_index, type):
    """Copy the diagram feynmanmp instructions in the result file."""
    if type == 'diag':
        diag_file = open(directory+"/Diagrams/diag_%i.tex" % diagram_index)
    elif type == 'time':
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


def tree_time_structure_den(time_diagram):
    """Return the denominator associated to a tree time-structure diagram."""
    denominator = ""
    labels = list(string.ascii_lowercase)
    i = 0
    for vertex in time_diagram:
        if not time_diagram.node[vertex]['operator']:
            time_diagram.node[vertex]['label'] = labels[i]
            i += 1
    for vertex in time_diagram:
        if not time_diagram.node[vertex]['operator']:
            if time_diagram.out_degree(vertex) == 0:
                denominator += time_diagram.node[vertex]['label']
            else:
                denominator += "(" + time_diagram.node[vertex]['label']
                for descendant in nx.descendants(time_diagram, vertex):
                    denominator += "+" + time_diagram.node[descendant]['label']
                denominator += ")"
    return denominator
