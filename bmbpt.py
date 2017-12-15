"""Module containg methods relative to BMBPT, to be called by ADG."""

import copy
import itertools
import string
import numpy as np
import methods as mth
import networkx as nx


def BMBPT_generation(p_order, three_N_use):
    """Generate diagrams for BMBPT from bottom up."""
    deg_max = 6 if three_N_use else 4

    # Create a null oriented adjacency matrix of dimension (p_order,p_order)
    matrices = [[[0 for element in range(p_order)] for line in range(p_order)]]

    # Generate oriented adjacency matrices going vertex-wise
    vertices = range(p_order)
    add = matrices.append
    for vertex in vertices:
        for sum_index in xrange(vertex+1, p_order):
            for mat_indx in xrange(len(matrices)-1, -1, -1):
                mat = matrices[mat_indx]
                elem_max = deg_max - sum(mat[k][vertex] + mat[vertex][k]
                                         for k in vertices)
                for elem in xrange(1, elem_max + 1, 1):
                    temp_mat = copy.deepcopy(mat)
                    temp_mat[vertex][sum_index] = elem
                    add(temp_mat)
        mth.check_vertex_degree(matrices, three_N_use, vertex)
        if 0 < vertex < p_order-1:
            check_unconnected_spawn(matrices, vertex, p_order)

    # Checks to exclude non-conform matrices
    mth.check_degree(matrices, three_N_use)
    mth.no_loop(matrices)
    matricesUniq = []
    for mat in matrices:
        if mat not in matricesUniq:
            matricesUniq.append(mat)
    matricesUniq.sort(reverse=True)
    return [np.array(mat) for mat in matricesUniq]


def check_unconnected_spawn(matrices, max_filled_vertex, length_mat):
    """Exclude some matrices that would spawn unconnected diagrams."""
    empty_block = [0 for i in range(length_mat - max_filled_vertex - 1)]
    for ind_mat in xrange(len(matrices)-1, -1, -1):
        mat = matrices[ind_mat]
        is_disconnected = True
        empty_lines = [index for index, line
                       in enumerate(mat[0:max_filled_vertex + 1])
                       if line[max_filled_vertex + 1:length_mat] == empty_block]
        test_block = [0 for i in range(length_mat - len(empty_lines))]
        for index in empty_lines:
            test_line = copy.deepcopy(mat[index])
            for index2 in empty_lines:
                test_line.remove(mat[index][index2])
            if test_line != test_block:
                is_disconnected = False
                break
        if is_disconnected and empty_lines != []:
            for index, line in enumerate(mat[0:max_filled_vertex + 1]):
                if index not in empty_lines:
                    for index2 in empty_lines:
                        if line[index2] != 0:
                            is_disconnected = False
                            break
            if is_disconnected:
                del matrices[ind_mat]


def order_2B_or_3B(graphs, TwoB_diags, ThreeB_diags):
    """Order the graphs depending on their 2B or 3B effective status."""
    for graph in graphs:
        max_deg = 0
        max_deg = (max(max_deg, graph.degree(node)) for node in graph)
        if max_deg == 6:
            ThreeB_diags.append(graph)
        else:
            TwoB_diags.append(graph)


def order_HF_or_not(graphs, HF_graphs, EHF_graphs, noHF_graphs, norm):
    """Order the graphs depending on their HF status."""
    for graph in graphs:
        test_HF = True
        test_EHF = True
        for node in graph:
            if graph.degree(node) == 2:
                test_HF = False
                if node != 0:
                    test_EHF = False
        if test_HF:
            HF_graphs.append(graph)
        elif (not test_EHF) or norm:
            noHF_graphs.append(graph)
        else:
            EHF_graphs.append(graph)


def attribute_qp_labels(graph):
    """Attribute the appropriate qp labels to the graph's propagators."""
    i = 1
    for prop in graph.edges(keys=True):
        graph.adj[prop[0]][prop[1]][prop[2]]['qp_state'] = "k_{%i}" % i
        i += 1


def omega_subgraph(graph):
    """Return the graph without any operator vertex."""
    subgraph_stack = [vertex for vertex in graph
                      if graph.node[vertex]['operator'] is False]
    return graph.subgraph(subgraph_stack)


def extract_numerator(graph):
    """Return the numerator associated to a BMBPT graph."""
    numerator = ""
    for vertex in graph:
        # Attribute the correct operator to each vertex
        numerator += "O" if graph.node[vertex]['operator'] else "\\Omega"
        # Attribute the good "type number" to each vertex
        numerator += "^{%i%i}_{" % (graph.out_degree(vertex),
                                    graph.in_degree(vertex))
        # First add the qp states corresponding to propagators going out
        numerator += "".join(graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                             for prop in graph.out_edges(vertex, keys=True))
        # Add the qp states corresponding to propagators coming in
        previous_vertex = vertex - 1
        while previous_vertex >= 0:
            numerator += "".join(graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                                 for prop in graph.in_edges(vertex, keys=True)
                                 if prop[0] == previous_vertex)
            previous_vertex -= 1
        numerator += "} "
    return numerator


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


def time_tree_denominator(graph, time_graph):
    """Return the denominator for a time-tree graph."""
    denominator = ""
    for vertex_i in range(1, len(time_graph)):
        subgraph_stack = [vertex_j
                          for vertex_j in nx.descendants(time_graph, vertex_i)]
        subgraph_stack.append(vertex_i)
        subdiag = graph.subgraph(subgraph_stack)
        denominator += "%s\\ " % extract_denom(graph, subdiag)
    return denominator


def extract_integral(diag):
    """Return the integral part of the Feynman expression of the diag."""
    graph = diag.graph
    pert_vertex_indices = range(1, len(graph))
    integral = "".join("\\mathrm{d}\\tau_%i" % vertex
                       for vertex in pert_vertex_indices)
    if len(graph) > 2:
        for vertex_i in pert_vertex_indices:
            integral += "".join("\\theta(\\tau_%i-\\tau_%i) " % (vertex_j,
                                                                 vertex_i)
                                for vertex_j in pert_vertex_indices
                                if graph.has_edge(vertex_i, vertex_j))
    integral += "".join("e^{-\\tau_%i %s}" % (vertex, diag.vert_exp[vertex])
                        for vertex in pert_vertex_indices)
    return integral


def extract_BMBPT_crossing_sign(graph):
    """Return True if there's a minus sign associated with crossing propagators.

    Use the fact that all lines propagate upwards and the
    canonical representation of the diagrams and vertices.
    """
    nb_crossings = 0
    for vertex in graph:
        for propagator in graph.out_edges(vertex, keys=True):
            for vertex_ante in xrange(propagator[0]):
                for vertex_post in xrange(propagator[0]+1, propagator[1]):
                    nb_crossings += graph.number_of_edges(vertex_ante,
                                                          vertex_post)
    return nb_crossings % 2 == 1


def multiplicity_symmetry_factor(graph):
    """Return the symmetry factor associated with propagators multiplicity."""
    factor = ""
    prop_multiplicity = [0 for i in xrange(6)]
    for vertex_i in graph:
        for vertex_j in graph:
            if graph.number_of_edges(vertex_i, vertex_j) >= 2:
                prop_multiplicity[graph.number_of_edges(
                    vertex_i, vertex_j) - 1] += 1

    for prop_id, multiplicity in enumerate(prop_multiplicity):
        if multiplicity == 1:
            factor += "(%i!)" % (prop_id+1)
        elif multiplicity >= 2:
            factor += "(%i!)" % (prop_id+1) + "^%i" % multiplicity
    return factor


def vertex_exchange_sym_factor(diag):
    """Return the symmetry factor associated with vertex exchange."""
    # Starts at -2 as the identity belongs to the set of permutations
    factor = -2
    graph = diag.graph
    perm_vertices = [vertex for vertex, degrees
                     in enumerate(diag.unsort_io_degrees)
                     if graph.node[vertex]['operator'] is False
                     and diag.unsort_io_degrees.count(degrees) >= 2]
    for permutation in itertools.permutations(perm_vertices):
        permuted_graph = nx.relabel_nodes(graph,
                                          dict(zip(perm_vertices,
                                                   permutation)),
                                          copy=True)
        if nx.is_isomorphic(graph, nx.intersection(graph, permuted_graph)):
            factor += 2
    return "%i" % factor if factor != 0 else ""


def write_BMBPT_header(tex_file, numdiag, three_N, norm, nb_2_HF,
                       nb_2_EHF, nb_2_noHF, nb_3_HF, nb_3_EHF, nb_3_noHF):
    """Write overall header for BMBPT result file."""
    tex_file.write("Valid diagrams: %i\n\n" % numdiag
                   + "2N valid diagrams: %i\n\n" % (nb_2_HF
                                                    + nb_2_EHF + nb_2_noHF)
                   + "2N canonical diagrams for the energy: %i\n\n" % nb_2_HF)
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
    if diag_index == 0:
        result.write("\\section{Two-body diagrams}\n\n"
                     + "\\subsection{Two-body energy canonical diagrams}\n\n")
    elif (diag_index == nb_2_HF) and (not norm):
        result.write("\\subsection{Two-body canonical diagrams for a generic operator only}\n\n")
    elif diag_index == nb_2_HF + nb_2_EHF:
        result.write("\\subsection{Two-body non-canonical diagrams}\n\n")
    if three_N:
        if diag_index == nb_2:
            result.write("\\section{Three-body diagrams}\n\n"
                         + "\\subsection{Three-body energy canonical diagrams}\n\n")
        elif (diag_index == nb_2 + nb_3_HF) and (not norm):
            result.write("\\subsection{Three-body canonical diagrams for a generic operator only}\n\n")
        elif diag_index == nb_2 + nb_3_HF + nb_3_EHF:
            result.write("\\subsection{Three-body non-canonical diagrams}\n\n")


def write_diag_exps(latex_file, bmbpt_diag, norder):
    """Write the expressions associated to a diagram in the LaTeX file."""
    latex_file.write("\\begin{align}\n\\text{PO}%i.%i\n" % (norder,
                                                            (bmbpt_diag.tags[0]
                                                             + 1))
                     + "&= %s" % bmbpt_diag.feynman_exp
                     + r" \nonumber \\" + "\n"
                     + "&= %s\\end{align}\n" % bmbpt_diag.diag_exp)


def write_vertices_values(latex_file, diag):
    """Write the qp energies associated to each vertex of the diag."""
    latex_file.write("\\begin{align*}\n")
    labels = list(string.ascii_lowercase)
    for ind in range(1, len(diag.vert_exp)):
        latex_file.write("%s &= %s" % (labels[ind-1], diag.vert_exp[ind]))
        if ind != len(diag.vert_exp)-1:
            latex_file.write(r"\\")
        latex_file.write('\n')
    latex_file.write("\\end{align*}\n")


class BmbptFeynmanDiagram(mth.Diagram):
    """Describes a BMBPT Feynman diagram with its related properties."""

    def __init__(self, nx_graph, use_norm, tag_num):
        """Generate a BMBPT diagrams using a NetworkX graph."""
        mth.Diagram.__init__(self, nx_graph)
        self.type = 'BMBPT'
        self.two_or_three_body = 3 if self.max_degree == 6 else 2
        self.tags = [tag_num]
        self.time_tag = -1
        self.tst_is_tree = False
        self.feynman_exp = ""
        self.diag_exp = ""
        self.vert_exp = []
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
        self.vert_exp = [self.vertex_expression(vertex)
                         for vertex in self.graph]
        numerator = extract_numerator(self.graph)
        denominator = time_tree_denominator(self.graph,
                                            time_diags[self.time_tag].graph) \
            if self.tst_ist_tree else ""
        extra_factor = "" if self.tst_is_tree \
            else "\\left[" \
            + " + ".join("\\frac{1}{%s}"
                         % time_tree_denominator(self.graph, equi_t_graph)
                         for equi_t_graph
                         in time_diags[self.time_tag].equivalent_trees) \
            + " \\right]"
        # Determine the pre-factor
        prefactor = "(-1)^%i " % (len(self.graph) - 1)
        if extract_BMBPT_crossing_sign(self.graph):
            prefactor = "-%s" % prefactor
        sym_fact = ""
        for vertex_degrees in self.unsort_io_degrees:
            if self.unsort_io_degrees.count(vertex_degrees) >= 2:
                sym_fact += vertex_exchange_sym_factor(self)
                break
        sym_fact += multiplicity_symmetry_factor(self.graph)
        prefactor = "\\frac{%s}{%s}\\sum_{k_i}" % (prefactor, sym_fact) \
            if sym_fact != "" else "%s\\sum_{k_i}" % prefactor
        # Set the Feynman and Goldstone expressions
        self.feynman_exp = "%s%s\\int_{0}^{\\tau}%s\n" \
                           % (prefactor, numerator, extract_integral(self))
        self.diag_exp = "%s\\frac{%s}{%s} %s\n" % (prefactor, numerator,
                                                   denominator, extra_factor) \
            if denominator != "" \
            else "%s%s%s\n" % (prefactor, numerator, extra_factor)

    def vertex_expression(self, vertex):
        """Return the expression associated to a given vertex."""
        expression = r"\epsilon^{" \
            + "".join("%s"
                      % self.graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                      for prop in self.graph.in_edges(vertex, keys=True)) \
            + "}_{" \
            + "".join("%s"
                      % self.graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                      for prop in self.graph.out_edges(vertex, keys=True)) \
            + "}"
        return expression
