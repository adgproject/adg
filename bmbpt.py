"""Module containg methods relative to BMBPT, to be called by ADG."""

import copy
import itertools
import string
import numpy as np
import methods as mth
import networkx as nx


def BMBPT_generation(p_order, three_N_use, norm_diagrams):
    """Generate diagrams for BMBPT from bottom up."""
    deg_max = 6 if three_N_use else 4

    # Create a null oriented adjacency matrix of dimension (p_order,p_order)
    matrices = [mth.empty_matrix_generation(p_order)]

    # Generate oriented adjacency matrices going vertex-wise
    vertices = range(p_order)
    for vertex in vertices:
        for sum_index in xrange(vertex+1, p_order):
            for mat_indx in xrange(len(matrices)-1, -1, -1):
                mat = matrices[mat_indx]
                if mat[sum_index][vertex] == 0:
                    vert_degree = sum(mat[k][vertex]
                                      + mat[vertex][k] for k in vertices)
                    elem = 1
                    while (elem + vert_degree) <= deg_max:
                        temp_mat = copy.deepcopy(mat)
                        temp_mat[vertex][sum_index] = elem
                        matrices.append(temp_mat)
                        elem += 1
        mth.check_vertex_degree(matrices, three_N_use, vertex)

    # Checks to exclude non-conform matrices
    mth.check_degree(matrices, three_N_use)
    mth.no_loop(matrices)
    matricesUniq = []
    for mat in matrices:
        if mat not in matricesUniq:
            matricesUniq.append(mat)
    matricesUniq.sort(reverse=True)
    return [np.array(mat) for mat in matricesUniq]


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
    subgraph_stack = [vertex for vertex in diagram
                      if diagram.node[vertex]['operator'] is False]
    return diagram.subgraph(subgraph_stack)


def extract_numerator(diagram):
    """Return the numerator associated to a BMBPT diagram."""
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
        subgraph_stack = [diagram.nodes()[vertex_i]]
        for vertex_j in nx.descendants(time_diagram, vertex_i):
            subgraph_stack.append(diagram.nodes()[vertex_j])
        subdiag = diagram.subgraph(subgraph_stack)
        denominator += "(" + extract_denom(diagram, subdiag) + ")"
    return denominator


def extract_integral(diagram):
    """Return the integral part of the Feynman expression of the diagram."""
    integral = ""
    norder = diagram.number_of_nodes()
    pert_vertex_indices = range(1, norder)
    for vertex in pert_vertex_indices:
        integral += "\\mathrm{d}\\tau_%i" % vertex
    if norder > 2:
        for vertex_i in pert_vertex_indices:
            for vertex_j in pert_vertex_indices:
                if diagram.has_edge(vertex_i, vertex_j):
                    integral += "\\theta(\\tau_%i" % vertex_j \
                        + "-\\tau_%i) " % vertex_i
    for vertex in pert_vertex_indices:
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
            for vertex_ante in xrange(propagator[0]):
                for vertex_post in xrange(propagator[0]+1, propagator[1]):
                    nb_crossings += diagram.number_of_edges(vertex_ante,
                                                            vertex_post)
    return nb_crossings % 2 == 1


def multiplicity_symmetry_factor(diagram):
    """Return the symmetry factor associated with propagators multiplicity."""
    factor = ""
    prop_multiplicity = [0 for i in xrange(6)]
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
    non_op_vertices = [vertex for vertex in diagram
                       if diagram.node[vertex]['operator'] is False]
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


def write_diag_exps(latex_file, bmbpt_diag, norder):
    """Write the expressions associated to a diagram in the LaTeX file."""
    latex_file.write("\\begin{align}\n\\text{PO}%i" % norder
                     + ".%i\n" % (bmbpt_diag.tags[0] + 1))
    latex_file.write("&= " + bmbpt_diag.feynman_exp + r" \nonumber \\" + "\n")
    latex_file.write("&= " + bmbpt_diag.diag_exp)
    latex_file.write("\\end{align}\n")


def write_vertices_values(latex_file, diagram):
    """Write the qp energies associated to each vertex of the diagram."""
    latex_file.write("\\begin{align*}\n")
    labels = list(string.ascii_lowercase)
    for vertex in xrange(1, len(diagram)):
        latex_file.write(labels[vertex-1] + " &= ")
        for prop in diagram.in_edges_iter(vertex, keys=True):
            latex_file.write(" + E_{%s}"
                             % diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state'])
        for prop in diagram.out_edges_iter(vertex, keys=True):
            latex_file.write(" - E_{%s}"
                             % diagram.edge[prop[0]][prop[1]][prop[2]]['qp_state'])
        if vertex != len(diagram)-1:
            latex_file.write(r"\\")
        latex_file.write('\n')
    latex_file.write("\\end{align*}\n")


class BmbptFeynmanDiagram(mth.Diagram):
    """Describes a BMBPT Feynman diagram with its related properties."""

    def __init__(self, nx_graph, use_norm, tag_num):
        mth.Diagram.__init__(self, nx_graph)
        self.two_or_three_body = 3 if self.max_degree == 6 else 2
        self.tags = [tag_num]
        self.time_tag = -1
        self.tst_is_tree = False
        self.feynman_exp = ""
        self.diag_exp = ""
        if 2 not in self.degrees:
            self.HF_type = "HF"
        elif use_norm:
            self.HF_type = "noHF"
        else:
            for node in xrange(1, len(self.graph)):
                if self.graph.degree(node) == 2:
                    self.HF_type = "noHF"
                    break
                if self.graph.degree(len(self.graph)-1) != 2:
                    self.HF_type = "EHF"

    def attribute_expressions(self, time_diags):
        """Attribute the correct Feynman and Goldstone expressions."""
        norder = len(self.graph)
        numerator = extract_numerator(self.graph)
        denominator = ""
        extra_factor = ""
        if self.tst_is_tree:
            for tdiag in time_diags:
                if tdiag.tags[0] == self.time_tag:
                    time_graph = tdiag.graph
                    break
            denominator = time_tree_denominator(self.graph,
                                                time_graph, denominator)

        elif (norder == 4) and (mth.number_of_sinks(self.graph) == 1):
            testdiag = omega_subgraph(self.graph)
            for i in range(2):
                subgraph_stack = [nx.dag_longest_path(testdiag)[1]]
                if i == 0:
                    subgraph_stack.append(nx.dag_longest_path(testdiag)[0])
                else:
                    for vertex_1 in testdiag:
                        test_vertex = True
                        for vertex_2 in nx.dag_longest_path(testdiag):
                            if vertex_1 == vertex_2:
                                test_vertex = False
                        if test_vertex:
                            subgraph_stack.append(vertex_1)
                subdiag = testdiag.subgraph(subgraph_stack)
                denominator += "(" + extract_denom(self.graph, subdiag) + ")"
            for vertex in self.graph:
                if self.graph.out_degree(vertex) == 0:
                    subdiag = self.graph.subgraph(vertex)
            denominator_a = extract_denom(self.graph, subdiag)
            denominator_abc = extract_denom(self.graph, testdiag)
            extra_factor += "\\left[ \\frac{1}{" + denominator_a \
                + "} + \\frac{1}{" + denominator_abc + "} \\right]"
        # Determine the pre-factor
        prefactor = "(-1)^%i " % (norder - 1)
        if extract_BMBPT_crossing_sign(self.graph):
            prefactor = "-" + prefactor
        sym_fact = vertex_exchange_sym_factor(self.graph) \
            + multiplicity_symmetry_factor(self.graph)
        if sym_fact != "":
            prefactor = "\\frac{" + prefactor + "}{" \
                + sym_fact + "}\\sum_{k_i}"
        else:
            prefactor = prefactor + "\\sum_{k_i}"
        self.feynman_exp = prefactor + numerator + "\\int_{0}^{\\tau}" \
            + extract_integral(self.graph) + "\n"
        if denominator != "":
            self.diag_exp = prefactor + "\\frac{ " + numerator \
                + " }{ " + denominator + " }" + extra_factor + "\n"
        else:
            self.diag_exp = prefactor + numerator + extra_factor + "\n"
