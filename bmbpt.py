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
    matrices = [mth.empty_matrix_generation(p_order)]

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

    # Checks to exclude non-conform matrices
    mth.check_degree(matrices, three_N_use)
    mth.no_loop(matrices)
    matricesUniq = []
    for mat in matrices:
        if mat not in matricesUniq:
            matricesUniq.append(mat)
    matricesUniq.sort(reverse=True)
    return [np.array(mat) for mat in matricesUniq]


def order_2B_or_3B(graphs, TwoB_diags, ThreeB_diags):
    """Order the graphs depending on their 2B or 3B effective status."""
    for graph in graphs:
        max_deg = 0
        for node in graph:
            max_deg = max(max_deg, graph.degree(node))
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
    for prop in graph.edges_iter(keys=True):
        graph.edge[prop[0]][prop[1]][prop[2]]['qp_state'] = "k_{%i}" % i
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
        if graph.node[vertex]['operator']:
            numerator += "O"
        else:
            numerator += "\\Omega"
        # Attribute the good "type number" to each vertex
        numerator = numerator + "^{%i" % graph.out_degree(vertex) \
            + "%i}_{" % graph.in_degree(vertex)
        # First add the qp states corresponding to propagators going out
        numerator += "".join(graph.edge[prop[0]][prop[1]][prop[2]]['qp_state']
                             for prop in graph.out_edges_iter(vertex,
                                                              keys=True))
        # Add the qp states corresponding to propagators coming in
        previous_vertex = vertex - 1
        while previous_vertex >= 0:
            numerator += "".join(graph.edge[prop[0]][prop[1]][prop[2]]['qp_state']
                                 for prop in graph.in_edges_iter(vertex,
                                                                 keys=True)
                                 if prop[0] == previous_vertex)
            previous_vertex -= 1
        numerator += "} "
    return numerator


def extract_denom(start_graph, subgraph):
    """Extract the appropriate denominator using the subgraph rule."""
    denomin = "".join(" + E_{%s}"
                      % start_graph.edge[propa[0]][propa[1]][propa[2]]['qp_state']
                      for propa
                      in start_graph.in_edges_iter(subgraph, keys=True)
                      if not subgraph.has_edge(propa[0], propa[1], propa[2]))
    denomin += "".join(" - E_{%s}"
                       % start_graph.edge[propa[0]][propa[1]][propa[2]]['qp_state']
                       for propa
                       in start_graph.out_edges_iter(subgraph, keys=True)
                       if not subgraph.has_edge(propa[0], propa[1], propa[2]))
    return denomin


def time_tree_denominator(graph, time_graph, denominator):
    """Add the denominator for a time-tree graph."""
    for vertex_i in range(1, len(time_graph)):
        subgraph_stack = [graph.nodes()[vertex_j]
                          for vertex_j in nx.descendants(time_graph, vertex_i)]
        subgraph_stack.append(graph.nodes()[vertex_i])
        subdiag = graph.subgraph(subgraph_stack)
        denominator += "(%s)" % extract_denom(graph, subdiag)
    return denominator.replace("( +", "(")


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
    """Return True if there's a sign factor associated with crossing propagators.

    Use the fact that all lines propagate upwards and the
    canonical representation of the diagrams and vertices.
    """
    nb_crossings = 0
    for vertex in graph:
        for propagator in graph.out_edges_iter(vertex, keys=True):
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
                                          dict(zip(perm_vertices, permutation)),
                                          copy=True)
        if nx.is_isomorphic(graph, nx.intersection(graph, permuted_graph)):
            factor += 2
    return "%i" % factor if factor != 0 else ""


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


def write_vertices_values(latex_file, diag):
    """Write the qp energies associated to each vertex of the diag."""
    latex_file.write("\\begin{align*}\n")
    labels = list(string.ascii_lowercase)
    for ind in range(1, len(diag.vert_exp)):
        latex_file.write("%s &= " % labels[ind-1])
        latex_file.write(diag.vert_exp[ind].replace("(", "").replace(")", ""))
        if ind != len(diag.vert_exp)-1:
            latex_file.write(r"\\")
        latex_file.write('\n')
    latex_file.write("\\end{align*}\n")


class BmbptFeynmanDiagram(mth.Diagram):
    """Describes a BMBPT Feynman diagram with its related properties."""

    def __init__(self, nx_graph, use_norm, tag_num):
        """Generate a BMBPT diagrams using a NetworkX graph."""
        mth.Diagram.__init__(self, nx_graph)
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
        self.attribute_vertices_expressions()
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
                denominator += "(%s)" % extract_denom(self.graph, subdiag)
            for vertex in self.graph:
                if self.graph.out_degree(vertex) == 0:
                    subdiag = self.graph.subgraph(vertex)
            extra_factor = "\\left[ \\frac{1}{%s} + \\frac{1}{%s} \\right]" \
                % (extract_denom(self.graph, subdiag),
                   extract_denom(self.graph, testdiag))
        # Determine the pre-factor
        prefactor = "(-1)^%i " % (norder - 1)
        if extract_BMBPT_crossing_sign(self.graph):
            prefactor = "-%s" % prefactor
        sym_fact = ""
        for vertex_degrees in self.unsort_io_degrees:
            if self.unsort_io_degrees.count(vertex_degrees) >= 2:
                sym_fact += vertex_exchange_sym_factor(self)
                break
        sym_fact += multiplicity_symmetry_factor(self.graph)
        if sym_fact != "":
            prefactor = "\\frac{%s}{%s}\\sum_{k_i}" % (prefactor, sym_fact)
        else:
            prefactor = "%s\\sum_{k_i}" % prefactor
        self.feynman_exp = "%s%s\\int_{0}^{\\tau}%s\n" % (prefactor, numerator,
                                                          extract_integral(self))
        if denominator != "":
            self.diag_exp = "%s\\frac{%s}{%s} %s\n" % (prefactor,
                                                       numerator,
                                                       denominator,
                                                       extra_factor)
        else:
            self.diag_exp = "%s%s%s\n" % (prefactor, numerator, extra_factor)

    def attribute_vertices_expressions(self):
        """Attribute the appropriate expression to each vertex."""
        vertices_expressions = []
        for vertex in self.graph:
            v_exp = "(" \
                + "".join(" + E_{%s}" % self.graph.edge[prop[0]][prop[1]][prop[2]]['qp_state']
                          for prop
                          in self.graph.in_edges_iter(vertex, keys=True)) \
                + "".join(" - E_{%s}" % self.graph.edge[prop[0]][prop[1]][prop[2]]['qp_state']
                          for prop
                          in self.graph.out_edges_iter(vertex, keys=True)) \
                + ")"
            if "( +" in v_exp:
                v_exp = v_exp.replace("( +", "(")
            vertices_expressions.append(v_exp)
        self.vert_exp = vertices_expressions
