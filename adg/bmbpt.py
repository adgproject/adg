"""Routines and class for Bogoliubov MBPT diagrams."""

import copy
import itertools
import numpy as np
import networkx as nx
import adg.tsd
import adg.diag


def diagrams_generation(p_order, three_body_use, nbody_obs, canonical):
    """Generate diagrams for BMBPT from bottom up.

    Args:
        p_order (int): The BMBPT perturbative order of the studied diagrams.
        three_body_use (bool): Flag for the use of three-body forces.
        nbody_obs (int): N-body character of the obervable of interest.
        canonical (bool): ``True`` if one draws only canonical diagrams.

    Returns:
        (list): NumPy arrays encoding the adjacency matrices of the graphs.

    >>> diagrams_generation(1, False, 2, False) #doctest: +NORMALIZE_WHITESPACE
    [array([[0, 4], [0, 0]]), array([[0, 2], [0, 0]])]
    >>> diagrams_generation(1, True, 3, False)  #doctest: +NORMALIZE_WHITESPACE
    [array([[0, 6], [0, 0]]), array([[0, 4], [0, 0]]), array([[0, 2], [0, 0]])]
    >>> diagrams_generation(2, False, 2, True)  #doctest: +NORMALIZE_WHITESPACE
    [array([[0, 2, 2], [0, 0, 2], [0, 0, 0]]),
     array([[0, 1, 1], [0, 0, 3], [0, 0, 0]])]

    """
    # Matrices contain operator vertex + p_order perturbative vertices
    order = p_order + 1

    # Create a null oriented adjacency matrix of dimension (p_order,p_order)
    matrices = [[[0 for _ in range(order)] for _ in range(order)]]

    # Generate oriented adjacency matrices going vertex-wise
    vertices = range(order)
    add = matrices.append
    for vertex in vertices:
        if vertex == 0:
            deg_max = 2*nbody_obs
        else:
            deg_max = 6 if three_body_use else 4
        for sum_index in xrange(vertex+1, order):
            for mat_indx in xrange(len(matrices)-1, -1, -1):
                mat = matrices[mat_indx]
                elem_max = deg_max - sum(mat[k][vertex] + mat[vertex][k]
                                         for k in vertices)
                for elem in xrange(1, elem_max + 1, 1):
                    temp_mat = copy.deepcopy(mat)
                    temp_mat[vertex][sum_index] = elem
                    add(temp_mat)
        adg.diag.check_vertex_degree(
            matrices, three_body_use, nbody_obs, canonical, vertex
        )
        if 0 < vertex < order-1:
            check_unconnected_spawn(matrices, vertex, order)

    matrices.sort(reverse=True)
    return [np.array(matrix) for matrix in matrices]


def check_unconnected_spawn(matrices, max_filled_vertex, length_mat):
    """Exclude some matrices that would spawn unconnected diagrams.

    Args:
        matrices (list): The adjacency matrices to be checked.
        max_filled_vertex (int): The furthest vertex until which the matrices
            have been filled.
        length_mat (int): The size of the square matrices.

    >>> mats = [[[0, 2, 0], [2, 0, 0], [0, 0, 0]], \
                [[0, 2, 1], [2, 0, 1], [0, 0, 0]]]
    >>>
    >>> check_unconnected_spawn(mats, 1, 3)
    >>> mats
    [[[0, 2, 1], [2, 0, 1], [0, 0, 0]]]

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


def write_header(tex_file, commands, diags_nbs):
    """Write overall header for BMBPT result file.

    Args:
        tex_file (file): The ouput LaTeX file of the program.
        commands (Namespace): Flags for the program run.
        diags_nbs (dict): The number of diagrams per type.

    """
    tex_file.write(
        "Valid diagrams: %i\n\n" % diags_nbs['nb_diags']
        + "2N valid diagrams: %i\n\n" % diags_nbs['nb_2']
        + "2N canonical diagrams for the energy: %i\n\n" % diags_nbs['nb_2_hf']
        + "2N canonical diagrams for a generic operator only: %i\n\n"
        % diags_nbs['nb_2_ehf']
    )
    if not commands.canonical:
        tex_file.write(
            "2N non-canonical diagrams: %i\n\n" % diags_nbs['nb_2_not_hf']
        )
    if commands.with_3NF:
        tex_file.write(
            "3N valid diagrams: %i\n\n" % diags_nbs['nb_3_hf']
            + "3N canonical diagrams for the energy: %i\n\n"
            % diags_nbs['nb_3_hf']
            + "3N canonical diagrams for a generic operator only: %i\n\n"
            % diags_nbs['nb_3_ehf']
        )
        if not commands.canonical:
            tex_file.write(
                "3N non-canonical diagrams: %i\n\n" % diags_nbs['nb_3_not_hf']
            )


def produce_expressions(diagrams, diagrams_time):
    """Produce and store the expressions associated to the BMBPT diagrams.

    Args:
        diagrams (list): The list of all BmbptFeynmanDiagrams.
        diagrams_time (list): Their associates TSDs.

    """
    for diag in diagrams:
        diag.attribute_qp_labels()
        for t_diag in diagrams_time:
            if diag.tags[0] in t_diag.tags[1:]:
                diag.time_tag = t_diag.tags[0]
                diag.tsd_is_tree = t_diag.is_tree
                break
        diag.attribute_expressions(diagrams_time[diag.time_tag])


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
    """Describes a BMBPT Feynman diagram with its related properties.

    Attributes:
        two_or_three_body (int): The 2 or 3-body characted of the vertices.
        time_tag (int): The tag number associated to the diagram's
            associated TSD.
        tsd_is_tree (bool): The tree or non-tree character of the
            associated TSD.
        feynman_exp (str): The Feynman expression associated to the diagram.
        diag_exp (str): The Goldstone expression associated to the diagram.
        vert_exp (list): The expression associated to the vertices.
        hf_type (str): The Hartree-Fock, non-Hartree-Fock or Hartree-Fock for
            the energy operator only character of the graph.

        """

    def __init__(self, nx_graph, tag_num):
        """Generate a BMBPT diagrams using a NetworkX graph.

        Args:
            nx_graph (NetworkX MultiDiGraph): The graph of interest.
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
        elif 2 not in self.unsort_degrees[1:]:
            self.hf_type = "EHF"
        else:
            self.hf_type = "noHF"

    def attribute_expressions(self, time_diag):
        """Attribute the correct Feynman and Goldstone expressions.

        Args:
            time_diag (TimeStructureDiagram): The associated TSD.

        """
        self.vert_exp = [self.vertex_expression(vertex)
                         for vertex in self.graph]
        numerator = self.extract_numerator()
        denominator = self.time_tree_denominator(
            nx.relabel_nodes(time_diag.graph, time_diag.perms[self.tags[0]])
        ) if self.tsd_is_tree else ""

        extra_factor = "" if self.tsd_is_tree \
            else "\\left[" \
            + " + ".join("\\frac{1}{%s}"
                         % self.time_tree_denominator(
                             nx.relabel_nodes(equi_t_graph,
                                              time_diag.perms[self.tags[0]]))
                         for equi_t_graph in time_diag.equivalent_trees) \
            + " \\right]"

        # Determine the pre-factor
        prefactor = "(-1)^%i " % (len(self.graph) - 1)
        if self.has_crossing_sign():
            prefactor = "-%s" % prefactor
        sym_fact = ""
        for vertex_degrees in self.unsort_io_degrees:
            if self.unsort_io_degrees.count(vertex_degrees) >= 2:
                sym_fact += self.vertex_exchange_sym_factor()
                break
        sym_fact += self.multiplicity_symmetry_factor()
        prefactor = "\\frac{%s}{%s}\\sum_{k_i}" % (prefactor, sym_fact) \
            if sym_fact != "" else "%s\\sum_{k_i}" % prefactor

        # Set the Feynman and Goldstone expressions
        self.feynman_exp = \
            "\\lim\\limits_{\\tau \\to \\infty}%s%s\\int_{0}^{\\tau}%s\n" \
            % (prefactor, numerator, self.extract_integral())
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
        self.write_vertices_values(latex_file, time_diag.perms[self.tags[0]])

    def write_section(self, result, commands, diags_nbs):
        """Write section and subsections for BMBPT result file.

        Args:
            result (file): The LaTeX output file of the program.
            commands (dict): The flags associated with run management.
            diags_nbs (dict): The number of diagrams per type.

        """
        if self.tags[0] == 0:
            result.write(
                "\\section{Two-body diagrams}\n\n"
                + "\\subsection{Two-body energy canonical diagrams}\n\n")
        elif self.tags[0] == diags_nbs['nb_2_hf']:
            result.write("\\subsection{Two-body canonical diagrams " +
                         "for a generic operator only}\n\n")
        elif self.tags[0] == diags_nbs['nb_2_hf'] + diags_nbs['nb_2_ehf']:
            result.write("\\subsection{Two-body non-canonical diagrams}\n\n")
        if commands.with_3NF:
            if self.tags[0] == diags_nbs['nb_2']:
                result.write(
                    "\\section{Three-body diagrams}\n\n"
                    + "\\subsection{Three-body energy canonical diagrams}\n\n")
            elif self.tags[0] == diags_nbs['nb_2'] + diags_nbs['nb_3_hf']:
                result.write("\\subsection{Three-body canonical diagrams " +
                             "for a generic operator only}\n\n")
            elif self.tags[0] == diags_nbs['nb_2'] + diags_nbs['nb_3_hf'] \
                    + diags_nbs['nb_3_ehf']:
                result.write(
                    "\\subsection{Three-body non-canonical diagrams}\n\n")
        result.write("\\paragraph{Diagram %i:}\n" % (self.tags[0] + 1))
        self.write_diag_exps(result, commands.order)

    def write_vertices_values(self, latex_file, mapping):
        """Write the qp energies associated to each vertex of the diag.

        Args:
            latex_file (file): The LaTeX output file of the program.
            mapping (dict): A mapping between the vertices in the diagram and
                the vertices in its euivalent TSD, since permutations between
                vertices are possible.

        """
        latex_file.write("\\begin{align*}\n")
        for ind in range(1, len(self.vert_exp)):
            latex_file.write("a_%i &= %s" % (ind, self.vert_exp[mapping[ind]]))
            if ind != len(self.vert_exp)-1:
                latex_file.write(r"\\")
            latex_file.write('\n')
        latex_file.write("\\end{align*}\n")

    def write_diag_exps(self, latex_file, norder):
        """Write the expressions associated to a diagram in the LaTeX file.

        Args:
            latex_file (file): The LaTeX outputfile of the program.
            norder (int): The order in BMBPT formalism.

        """
        latex_file.write(
            "\\begin{align}\n\\text{PO}%i.%i\n" % (norder, (self.tags[0] + 1))
            + "&= %s" % self.feynman_exp
            + r" \nonumber \\" + "\n"
            + "&= %s\\end{align}\n" % self.diag_exp)

    def vertex_exchange_sym_factor(self):
        """Return the symmetry factor associated with vertex exchange.

        Returns:
            (str): The symmetry factor for vertex exchange.

        """
        # Starts at -2 as the identity belongs to the set of permutations
        factor = -2
        graph = self.graph
        perm_vertices = [vertex for vertex, degrees
                         in enumerate(self.unsort_io_degrees)
                         if graph.node[vertex]['operator'] is False
                         and self.unsort_io_degrees.count(degrees) >= 2]
        for permutation in itertools.permutations(perm_vertices):
            permuted_graph = nx.relabel_nodes(graph,
                                              dict(zip(perm_vertices,
                                                       permutation)),
                                              copy=True)
            if nx.is_isomorphic(graph, nx.intersection(graph, permuted_graph)):
                factor += 2
        return "%i" % factor if factor != 0 else ""

    def extract_integral(self):
        """Return the integral part of the Feynman expression of the diag.

        Returns:
            (str): The integral part of its Feynman expression.

        """
        pert_vertex_indices = range(1, len(self.graph))
        integral = "".join("\\mathrm{d}\\tau_%i" % vertex
                           for vertex in pert_vertex_indices)
        if len(pert_vertex_indices) > 1:
            for vertex_i in pert_vertex_indices:
                integral += "".join("\\theta(\\tau_%i-\\tau_%i) " % (vertex_j,
                                                                     vertex_i)
                                    for vertex_j in pert_vertex_indices
                                    if self.graph.has_edge(vertex_i, vertex_j))
        integral += "".join("e^{-\\tau_%i %s}"
                            % (vertex, self.vert_exp[vertex])
                            for vertex in pert_vertex_indices)
        return integral

    def attribute_qp_labels(self):
        """Attribute the appropriate qp labels to the graph's propagators."""
        for idx, prop in enumerate(self.graph.edges(keys=True, data=True)):
            prop[3]['qp_state'] = "k_{%i}" % (idx+1)

    def extract_numerator(self):
        """Return the numerator associated to a BMBPT graph.

        Returns:
            (str): The numerator of the graph.

        """
        graph = self.graph
        numerator = ""
        for vertex in graph:
            # Attribute the correct operator to each vertex
            numerator += "O" if graph.node[vertex]['operator'] else "\\Omega"
            # Attribute the good "type number" to each vertex
            numerator += "^{%i%i}_{" % (self.unsort_io_degrees[vertex][1],
                                        self.unsort_io_degrees[vertex][0])
            # First add the qp states corresponding to propagators going out
            numerator += "".join(prop[3]['qp_state']
                                 for prop
                                 in graph.out_edges(vertex,
                                                    keys=True, data=True))
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

    def has_crossing_sign(self):
        """Return True for a minus sign associated with crossing propagators.

        Use the fact that all lines propagate upwards and the
        canonical representation of the diagrams and vertices.

        Returns:
            (bool): Encode for the sign factor associated with crossing
                propagators.

        """
        nb_crossings = 0
        for vertex in self.graph:
            for propagator in self.graph.out_edges(vertex, keys=True):
                for vertex_ante in xrange(propagator[0]):
                    for vertex_post in xrange(propagator[0]+1, propagator[1]):
                        nb_crossings += self.graph.number_of_edges(vertex_ante,
                                                                   vertex_post)
        return nb_crossings % 2 == 1

    def multiplicity_symmetry_factor(self):
        """Return the symmetry factor associated with propagators multiplicity.

        Returns:
            (str): The symmetry factor associated with equivalent lines.

        """
        factor = ""
        # Account for up to three-body operators
        prop_multiplicity = [0 for _ in xrange(6)]
        for vertex_i in self.graph:
            for vertex_j in self.graph:
                if self.graph.number_of_edges(vertex_i, vertex_j) >= 2:
                    prop_multiplicity[self.graph.number_of_edges(
                        vertex_i, vertex_j) - 1] += 1

        for prop_id, multiplicity in enumerate(prop_multiplicity):
            if multiplicity == 1:
                factor += "(%i!)" % (prop_id+1)
            elif multiplicity >= 2:
                factor += "(%i!)" % (prop_id+1) + "^%i" % multiplicity
        return factor

    def time_tree_denominator(self, time_graph):
        """Return the denominator for a time-tree graph.

        Args:
            time_graph (NetworkX MultiDiGraph): Its associated time-structure
                graph.

        Returns:
            (str): The denominator of the graph.

        """
        denominator = ""
        for vertex_i in range(1, len(time_graph)):
            subgraph_stack = [vertex_j for vertex_j
                              in nx.descendants(time_graph, vertex_i)]
            subgraph_stack.append(vertex_i)
            subdiag = self.graph.subgraph(subgraph_stack)
            denominator += "%s\\ " % adg.diag.extract_denom(self.graph,
                                                            subdiag)
        return denominator
