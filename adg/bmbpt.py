"""Routines and class for Bogoliubov MBPT diagrams."""

import copy
import itertools
import numpy as np
import networkx as nx
import adg.tsd
import adg.diag


def diagrams_generation(p_order, three_body_use):
    """Generate diagrams for BMBPT from bottom up.

    Args:
        p_order (int): The conventional order of the studied diagrams.
        three_body_use (bool): Flag for the use of three-body operators.

    Returns:
        (list): NumPy arrays encoding the adjacency matrices of the graphs.

    """
    deg_max = 6 if three_body_use else 4

    # Create a null oriented adjacency matrix of dimension (p_order,p_order)
    matrices = [[[0 for _ in range(p_order)] for _ in range(p_order)]]

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
        adg.diag.check_vertex_degree(matrices, three_body_use, vertex)
        if 0 < vertex < p_order-1:
            check_unconnected_spawn(matrices, vertex, p_order)

    # Checks to exclude non-conform matrices
    adg.diag.check_degree(matrices, three_body_use)
    adg.diag.no_loop(matrices)
    matrices_uniq = []
    for mat in matrices:
        if mat not in matrices_uniq:
            matrices_uniq.append(mat)
    matrices_uniq.sort(reverse=True)
    return [np.array(mat) for mat in matrices_uniq]


def check_unconnected_spawn(matrices, max_filled_vertex, length_mat):
    """Exclude some matrices that would spawn unconnected diagrams.

    Args:
        matrices (list): The adjacency matrices to be checked.
        max_filled_vertex (int): The furthest vertex until which the matrices
            have been filled.
        length_mat (int): The size of the square matrices.

    """
    empty_block = [0 for _ in range(length_mat - max_filled_vertex - 1)]
    for ind_mat in xrange(len(matrices)-1, -1, -1):
        mat = matrices[ind_mat]
        is_disconnected = True
        empty_lines = [index for index, line
                       in enumerate(mat[0:max_filled_vertex + 1])
                       if line[max_filled_vertex + 1:length_mat]
                       == empty_block]
        test_block = [0 for _ in range(length_mat - len(empty_lines))]
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
                    for _ in (idx for idx in empty_lines if line[idx] != 0):
                        is_disconnected = False
                        break
            if is_disconnected:
                del matrices[ind_mat]


def attribute_qp_labels(graph):
    """Attribute the appropriate qp labels to the graph's propagators.

    Args:
        graph(NetworkX MultiDiGraph): The graph of interest.

    """
    for idx, prop in enumerate(graph.edges(keys=True, data=True)):
        prop[3]['qp_state'] = "k_{%i}" % (idx+1)


def extract_numerator(graph):
    """Return the numerator associated to a BMBPT graph.

    Args:
        graph(NetworkX MultiDiGraph): The graph of interest.

    Returns:
        (str): The numerator of the graph.

    """
    numerator = ""
    for vertex in graph:
        # Attribute the correct operator to each vertex
        numerator += "O" if graph.node[vertex]['operator'] else "\\Omega"
        # Attribute the good "type number" to each vertex
        numerator += "^{%i%i}_{" % (graph.out_degree(vertex),
                                    graph.in_degree(vertex))
        # First add the qp states corresponding to propagators going out
        numerator += "".join(prop[3]['qp_state']
                             for prop
                             in graph.out_edges(vertex, keys=True, data=True))
        # Add the qp states corresponding to propagators coming in
        previous_vertex = vertex - 1
        while previous_vertex >= 0:
            numerator += "".join(
                prop[3]['qp_state']
                for prop in graph.in_edges(vertex, keys=True, data=True)
                if prop[0] == previous_vertex)
            previous_vertex -= 1
        numerator += "} "
    return numerator


def time_tree_denominator(graph, time_graph):
    """Return the denominator for a time-tree graph.

    Args:
        graph (NetworkX MultiDiGraph): The graph of interest.
        time_graph (NetworkX MultiDiGraph): Its associated time-structure
            graph.

    Returns:
        (str): The denominator of the graph.

    """
    denominator = ""
    for vertex_i in range(1, len(time_graph)):
        subgraph_stack = [vertex_j
                          for vertex_j in nx.descendants(time_graph, vertex_i)]
        subgraph_stack.append(vertex_i)
        subdiag = graph.subgraph(subgraph_stack)
        denominator += "%s\\ " % adg.diag.extract_denom(graph, subdiag)
    return denominator


def extract_integral(diag):
    """Return the integral part of the Feynman expression of the diag.

    Args:
        diag (BmbptFeynmanDiagram): The diagram of interest.

    Returns:
        (str): The integral part of its Feynman expression.

    """
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


def has_crossing_sign(graph):
    """Return True if there's a minus sign associated with crossing propagators.

    Use the fact that all lines propagate upwards and the
    canonical representation of the diagrams and vertices.

    Args:
        graph (NetworkX MultiDiGraph): The graph of interest.

    Returns:
        (bool): Encode for the sign factor associated with crossing
            propagators.

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
    """Return the symmetry factor associated with propagators multiplicity.

    Args:
        graph (NetworkX MultiDiGraph): The graph of interest.

    Returns:
        (str): The symmetry factor associated with equivalent lines.

    """
    factor = ""
    prop_multiplicity = [0 for _ in xrange(6)]
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
    """Return the symmetry factor associated with vertex exchange.

    Args:
        diag (BmbptFeynmanDiagram): The diagram of interest.

    Returns:
        (str): The symmetry factor for vertex exchange.

    """
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


def write_header(tex_file, three_body_use, norm, diags_nbs):
    """Write overall header for BMBPT result file.

    Args:
        tex_file (file): The ouput LaTeX file of the program.
        three_body_use (bool): True if one uses three-body operators.
        diags_nbs (dict): The number of diagrams per type.

    """
    tex_file.write(
        "Valid diagrams: %i\n\n" % diags_nbs['nb_diags']
        + "2N valid diagrams: %i\n\n" % diags_nbs['nb_2']
        + "2N canonical diagrams for the energy: %i\n\n"
        % diags_nbs['nb_2_hf'])
    if not norm:
        tex_file.write(
            "2N canonical diagrams for a generic operator only: %i\n\n"
            % diags_nbs['nb_2_ehf'])
    tex_file.write(
        "2N non-canonical diagrams: %i\n\n" % diags_nbs['nb_2_not_hf'])
    if three_body_use:
        tex_file.write(
            "3N valid diagrams: %i\n\n" % diags_nbs['nb_3_hf'])
        tex_file.write(
            "3N canonical diagrams for the energy: %i\n\n"
            % diags_nbs['nb_3_hf'])
        if not norm:
            tex_file.write(
                "3N canonical diagrams for a generic operator only: %i\n\n"
                % diags_nbs['nb_3_ehf'])
        tex_file.write(
            "3N non-canonical diagrams: %i\n\n" % diags_nbs['nb_3_not_hf'])


def write_section(result, diag, commands, diags_nbs):
    """Write section and subsections for BMBPT result file.

    Args:
        result (file): The LaTeX output file of the program.
        diag (BmbptFeynmanDiagram): The diagram of interest.
        commands (dict): The flags associated with run management.
        diags_nbs (dict): The number od diagrams per type.

    """
    if diag.tags[0] == 0:
        result.write("\\section{Two-body diagrams}\n\n"
                     + "\\subsection{Two-body energy canonical diagrams}\n\n")
    elif (diag.tags[0] == diags_nbs['nb_2_hf']) and (not commands.norm):
        result.write("\\subsection{Two-body canonical diagrams " +
                     "for a generic operator only}\n\n")
    elif diag.tags[0] == diags_nbs['nb_2_hf'] + diags_nbs['nb_2_ehf']:
        result.write("\\subsection{Two-body non-canonical diagrams}\n\n")
    if commands.with_three_body:
        if diag.tags[0] == diags_nbs['nb_2']:
            result.write(
                "\\section{Three-body diagrams}\n\n"
                + "\\subsection{Three-body energy canonical diagrams}\n\n")
        elif (diag.tags[0] == diags_nbs['nb_2'] + diags_nbs['nb_3_hf']) \
                and (not commands.norm):
            result.write("\\subsection{Three-body canonical diagrams " +
                         "for a generic operator only}\n\n")
        elif diag.tags[0] == diags_nbs['nb_2'] + diags_nbs['nb_3_hf'] \
                + diags_nbs['nb_3_ehf']:
            result.write("\\subsection{Three-body non-canonical diagrams}\n\n")
    result.write("\\paragraph{Diagram %i:}\n" % (diag.tags[0] + 1))
    if not commands.norm:
        write_diag_exps(result, diag, commands.order)


def write_diag_exps(latex_file, bmbpt_diag, norder):
    """Write the expressions associated to a diagram in the LaTeX file.

    Args:
        latex_file (file): The LaTeX outputfile of the program.
        bmbpt_diag (BmbptFeynmanDiagram): The diagram of interest.
        norder (int): The order in BMBPT formalism.

    """
    latex_file.write("\\begin{align}\n\\text{PO}%i.%i\n" % (norder,
                                                            (bmbpt_diag.tags[0]
                                                             + 1))
                     + "&= %s" % bmbpt_diag.feynman_exp
                     + r" \nonumber \\" + "\n"
                     + "&= %s\\end{align}\n" % bmbpt_diag.diag_exp)


def write_vertices_values(latex_file, diag, mapping):
    """Write the qp energies associated to each vertex of the diag.

    Args:
        latex_file (file): The LaTeX output file of the program.
        diag (BmbptFeynmanDiagram): The diagram of interest.
        mapping (dict): A mapping between the vertices in the diagram and the
            vertices in its euivalent TSD, since permutations between vertices
            are possible.

    """
    latex_file.write("\\begin{align*}\n")
    for ind in range(1, len(diag.vert_exp)):
        latex_file.write("a_%i &= %s" % (ind, diag.vert_exp[mapping[ind]]))
        if ind != len(diag.vert_exp)-1:
            latex_file.write(r"\\")
        latex_file.write('\n')
    latex_file.write("\\end{align*}\n")


def produce_expressions(diagrams, diagrams_time):
    """Produce and store the expressions associated to the BMBPT diagrams.

    Args:
        diagrams (list): The list of all BmbptFeynmanDiagrams.
        diagrams_time (list): Their associates TSDs.

    """
    for diag in diagrams:
        attribute_qp_labels(diag.graph)
        for t_diag in diagrams_time:
            if diag.tags[0] in t_diag.tags[1:]:
                diag.time_tag = t_diag.tags[0]
                diag.tsd_is_tree = t_diag.is_tree
                break
        diag.attribute_expressions(diagrams_time[diag.time_tag])


def treat_tsds(diagrams_time):
    """Order TSDs, produce their expressions, return also number of trees.

    Args:
        diagrams_time (list): All the associated TSDs.

    """
    tree_tsds = []
    for i_diag in xrange(len(diagrams_time)-1, -1, -1):
        if diagrams_time[i_diag].is_tree:
            tree_tsds.append(diagrams_time[i_diag])
            del diagrams_time[i_diag]

    adg.diag.topologically_distinct_diagrams(tree_tsds)
    adg.diag.topologically_distinct_diagrams(diagrams_time)

    diagrams_time = tree_tsds + diagrams_time

    for index, t_diag in enumerate(diagrams_time):
        t_diag.tags.insert(0, index)
        if not t_diag.is_tree:
            t_diag.equivalent_trees = adg.tsd.treat_cycles(t_diag.graph)
            t_diag.expr = " + ".join("\\frac{1}{%s}"
                                     % adg.tsd.tree_time_structure_den(graph)
                                     for graph
                                     in t_diag.equivalent_trees)
    return diagrams_time, len(tree_tsds)


def order_diagrams(diagrams):
    """Order the BMBPT diagrams and return number of diags for each type.

    Args:
        diagrams (list): Possibly redundant BmbptFeynmanDiagrams.

    Returns:
        (tuple): First element is the list of topologically unique, ordered
            diagrams. Second element is a dict with the number of diagrams
            for each major type.

    """
    diagrams_2_hf = []
    diagrams_2_ehf = []
    diagrams_2_not_hf = []
    diagrams_3_hf = []
    diagrams_3_ehf = []
    diagrams_3_not_hf = []

    for i_diag in xrange(len(diagrams)-1, -1, -1):
        if diagrams[i_diag].two_or_three_body == 2:
            if diagrams[i_diag].hf_type == "HF":
                diagrams_2_hf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "EHF":
                diagrams_2_ehf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "noHF":
                diagrams_2_not_hf.append(diagrams[i_diag])
        elif diagrams[i_diag].two_or_three_body == 3:
            if diagrams[i_diag].hf_type == "HF":
                diagrams_3_hf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "EHF":
                diagrams_3_ehf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "noHF":
                diagrams_3_not_hf.append(diagrams[i_diag])
        del diagrams[i_diag]

    adg.diag.topologically_distinct_diagrams(diagrams_2_hf)
    adg.diag.topologically_distinct_diagrams(diagrams_2_ehf)
    adg.diag.topologically_distinct_diagrams(diagrams_2_not_hf)
    adg.diag.topologically_distinct_diagrams(diagrams_3_hf)
    adg.diag.topologically_distinct_diagrams(diagrams_3_ehf)
    adg.diag.topologically_distinct_diagrams(diagrams_3_not_hf)

    diagrams = diagrams_2_hf + diagrams_2_ehf + diagrams_2_not_hf \
        + diagrams_3_hf + diagrams_3_ehf + diagrams_3_not_hf
    for ind, diagram in enumerate(diagrams):
        diagram.tags[0] = ind

    diags_nb_per_type = {
        'nb_2_hf': len(diagrams_2_hf),
        'nb_2_ehf': len(diagrams_2_ehf),
        'nb_2_not_hf': len(diagrams_2_not_hf),
        'nb_3_hf': len(diagrams_3_hf),
        'nb_3_ehf': len(diagrams_3_ehf),
        'nb_3_not_hf': len(diagrams_3_not_hf),
        'nb_diags': len(diagrams),
        'nb_2': (len(diagrams_2_hf) + len(diagrams_2_ehf)
                 + len(diagrams_2_not_hf)),
        'nb_3': (len(diagrams_3_hf) + len(diagrams_3_ehf)
                 + len(diagrams_3_not_hf))
        }

    return diagrams, diags_nb_per_type


class BmbptFeynmanDiagram(adg.diag.Diagram):
    """Describes a BMBPT Feynman diagram with its related properties."""

    def __init__(self, nx_graph, use_norm, tag_num):
        """Generate a BMBPT diagrams using a NetworkX graph.

        Args:
            nx_graph (NetworkX MultiDiGraph): The graph of interest.
            use_norm (bool): ``True`` if we are sturdying norm diagrams.
            tag_num (int): The tag number associated to the diagram.

        """
        adg.diag.Diagram.__init__(self, nx_graph)
        self.two_or_three_body = 3 if self.max_degree == 6 else 2
        self.tags = [tag_num]
        self.time_tag = -1
        self.tsd_is_tree = False
        self.feynman_exp = ""
        self.diag_exp = ""
        self.vert_exp = []
        if 2 not in self.degrees:
            self.hf_type = "HF"
        elif use_norm:
            self.hf_type = "noHF"
        else:
            for node in xrange(1, len(self.graph)):
                if self.graph.degree(node) == 2:
                    self.hf_type = "noHF"
                    break
                if self.graph.degree(len(self.graph)-1) != 2:
                    self.hf_type = "EHF"

    def attribute_expressions(self, time_diag):
        """Attribute the correct Feynman and Goldstone expressions.

        Args:
            time_diag (TimeStructureDiagram): The associated TSD.

        """
        self.vert_exp = [self.vertex_expression(vertex)
                         for vertex in self.graph]
        numerator = extract_numerator(self.graph)
        denominator = time_tree_denominator(
            self.graph,
            nx.relabel_nodes(time_diag.graph,
                             time_diag.perms[self.tags[0]])
        ) if self.tsd_is_tree else ""
        extra_factor = "" if self.tsd_is_tree \
            else "\\left[" \
            + " + ".join("\\frac{1}{%s}"
                         % time_tree_denominator(
                             self.graph,
                             nx.relabel_nodes(
                                 equi_t_graph,
                                 time_diag.perms[self.tags[0]]))
                         for equi_t_graph in time_diag.equivalent_trees) \
            + " \\right]"
        # Determine the pre-factor
        prefactor = "(-1)^%i " % (len(self.graph) - 1)
        if has_crossing_sign(self.graph):
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
        self.feynman_exp = \
            "\\lim\\limits_{\\tau \\to \\infty}%s%s\\int_{0}^{\\tau}%s\n" \
            % (prefactor, numerator, extract_integral(self))
        self.diag_exp = \
            "%s\\frac{%s}{%s} %s\n" % (prefactor, numerator,
                                       denominator, extra_factor) \
            if denominator != "" \
            else "%s%s%s\n" % (prefactor, numerator, extra_factor)

    def vertex_expression(self, vertex):
        """Return the expression associated to a given vertex.

        Args:
            vertex (int): The vertex of interest in the graph.

        """
        expression = r"\epsilon^{" \
            + "".join("%s"
                      % prop[3]['qp_state']
                      for prop
                      in self.graph.out_edges(vertex, keys=True, data=True)) \
            + "}_{" \
            + "".join("%s"
                      % prop[3]['qp_state']
                      for prop
                      in self.graph.in_edges(vertex, keys=True, data=True)) \
            + "}"
        return expression

    def write_graph(self, latex_file, directory, write_time):
        """Write the BMBPT graph and its associated TSD to the LaTeX file.

        Args:
            latex_file (file): The LaTeX output file of the program.
            directory (str): The path to the result folder.
            write_time (bool): ``True`` if we want informations on the
                associated TSDs.

        """
        latex_file.write('\n\\begin{center}\n')
        adg.diag.draw_diagram(directory, latex_file, self.tags[0], 'diag')
        if write_time:
            latex_file.write(
                '\\hspace{10pt} $\\rightarrow$ \\hspace{10pt} T%i:'
                % (self.time_tag + 1))
            adg.diag.draw_diagram(directory, latex_file, self.time_tag, 'time')
        latex_file.write('\n\\end{center}\n\n')

    def write_tsd_info(self, diagrams_time, latex_file):
        """Write info related to the BMBPT associated TSD to the LaTeX file.

        Args:
            diagrams_time (list): The associated TSDs.
            latex_file (file): The LaTeX output file of the program.

        """
        for tdiag in diagrams_time:
            if self.time_tag == tdiag.tags[0]:
                time_diag = tdiag
                break
        latex_file.write(
            "\\begin{equation}\n\\text{T}%i = " % (self.time_tag + 1)
            + "%s\\end{equation}\n" % time_diag.expr)
        write_vertices_values(latex_file, self, time_diag.perms[self.tags[0]])
