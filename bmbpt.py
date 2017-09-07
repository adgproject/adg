"""Module containg methods relative to BMBPT, to be called by ADG."""

import copy
import itertools
import numpy as np
import methods as mth
import networkx as nx


def BMBPT_generation(p_order, three_N_use, norm_diagrams):
    """Generate diagrams for BMBPT from bottom up."""
    deg_max = 4
    if three_N_use:
        deg_max = 6

    # Create a null oriented adjacency matrix of dimension (p_order,p_order)
    temp_matrices = []
    temp_matrices.append(mth.empty_matrix_generation(p_order))

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
        temp_matrices = mth.check_vertex_degree(matrices, three_N_use, vertex)
        matrices = copy.deepcopy(temp_matrices)

    # Checks to exclude non-conform matrices
    matrices = mth.check_degree(matrices, three_N_use)
    matrices = mth.no_loop(matrices)
    matricesUniq = []
    for mat in matrices:
        if mat not in matricesUniq:
            matricesUniq.append(mat)
    matricesUniq.sort(reverse=True)
    bmbpt_diagrams = []
    for mat in matricesUniq:
        bmbpt_diagrams.append(np.array(mat))
    return bmbpt_diagrams


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


def attribute_qp_labels(diagram):
    """Attribute the appropriate qp labels to the graph's propagators."""
    i = 1
    for prop in diagram.edges_iter(keys=True):
        diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state'] = "k_{%i}" % i
        i += 1


def omega_subgraph(diagram):
    """Return the graph without any operator vertex."""
    subgraph_stack = []
    for vertex in diagram:
        if diagram.node[vertex]['operator'] is False:
            subgraph_stack.append(vertex)
    return diagram.subgraph(subgraph_stack)


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
                    numerator += \
                        diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state']
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
    for permutation in itertools.permutations(non_op_vertices,
                                              len(non_op_vertices)):
        mapping = dict(zip(non_op_vertices, permutation))
        permuted_diag = nx.relabel_nodes(diagram, mapping, copy=True)
        if nx.is_isomorphic(diagram, nx.intersection(diagram, permuted_diag)):
            factor += 2
    if factor != 0:
        return "%i" % factor
    else:
        return ""


def write_BMBPT_header(tex_file, numdiag, three_N, norm, nb_2_HF,
                       nb_2_EHF, nb_2_noHF, nb_3_HF, nb_3_EHF, nb_3_noHF):
    """Write overall header for BMBPT result file."""
    tex_file.write("Valid diagrams: %i\n\n" % numdiag)
    tex_file.write("2N valid diagrams: %i\n\n" %
                   (nb_2_HF + nb_2_EHF + nb_2_noHF))
    tex_file.write("2N canonical diagrams for the energy: %i\n\n" % nb_2_HF)
    if not norm:
        tex_file.write(
            "2N canonical diagrams for a generic operator only: %i\n\n"
            % nb_2_EHF)
    tex_file.write("2N non-canonical diagrams: %i\n\n" % nb_2_noHF)
    if three_N:
        tex_file.write("3N valid diagrams: %i\n\n" %
                       (nb_3_HF + nb_3_EHF + nb_3_noHF))
        tex_file.write(
            "3N canonical diagrams for the energy: %i\n\n" % nb_3_HF)
        if not norm:
            tex_file.write(
                "3N canonical diagrams for a generic operator only: %i\n\n"
                % nb_3_EHF)
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
            result.write(
                "\\subsection{Three-body energy canonical diagrams}\n\n")
        elif (diag_index == nb_2 + nb_3_HF) and (not norm):
            result.write("\\subsection{Three-body canonical diagrams for a generic operator only}\n\n")
        elif diag_index == nb_2 + nb_3_HF + nb_3_EHF:
            result.write("\\subsection{Three-body non-canonical diagrams}\n\n")
