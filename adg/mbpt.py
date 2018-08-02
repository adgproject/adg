"""Routines and class for Many-Body Perturbation Theory diagrams."""

import copy
import itertools
import string
import numpy as np
import networkx as nx
import adg.diag


def diagrams_generation(order):
    """Generate the diagrams for the MBPT case.

    Args:
        order (int): The perturbative order of interest.

    Returns:
        (list): A list of NumPy arrays with the diagrams adjacency matrices.

    >>> diagrams_generation(2) # doctest: +NORMALIZE_WHITESPACE
    [array([[0, 2], [2, 0]])]
    >>> diagrams_generation(3) # doctest: +NORMALIZE_WHITESPACE
    [array([[0, 2, 0], [0, 0, 2], [2, 0, 0]]),
     array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]),
     array([[0, 0, 2], [2, 0, 0], [0, 2, 0]])]
    >>> diagrams_generation(1)
    []

    """
    # Generate all 1-magic square of dimension order
    seeds = [k for k in itertools.permutations(range(order), order)]
    all_matrices = [[[0 if i != j else 1 for i in range(order)]
                     for j in k]
                    for k in seeds]
    traceless = adg.diag.no_trace(all_matrices)
    coeffs = [i for i in itertools.combinations_with_replacement(
        range(len(traceless)), 2)]
    double = []

    for coef in coeffs:
        matrix = copy.deepcopy(traceless[coef[0]])
        for i, line in enumerate(traceless[coef[1]]):
            for j, elem in enumerate(line):
                matrix[i][j] += elem
        double.append(matrix)
    double_uniq = []
    for matrix in double:
        if matrix not in double_uniq:
            double_uniq.append(matrix)
    double_uniq.sort(reverse=True)
    return [np.array(matrix) for matrix in double_uniq]


def write_diag_exp(latex_file, mbpt_diag):
    """Write the expression associated to a diagram in the LaTeX file.

    Args:
        latex_file (file): The LaTeX output file to be written in.
        mbpt_diag (MbptDiagram): The diagram which expression is being written.

    """
    latex_file.write("\\begin{equation}\n")
    latex_file.write(mbpt_diag.expr)
    latex_file.write("\\end{equation}\n")


def write_header(tex_file, diags_nbs):
    """Write tha appropriate header for the LaTeX file for MBPT diagrams.

    Args:
        tex_file (file): The LaTeX ouput file to be written in.
        diags_nbs (dict): A dict with the number of diagrams per
            excitation level type.

    """
    tex_file.write("Valid diagrams: %i\n\n" % diags_nbs['nb_diags']
                   + "Singles: %i\n\n" % diags_nbs['singles']
                   + "Doubles: %i\n\n" % diags_nbs['doubles']
                   + "Triples: %i\n\n" % diags_nbs['triples']
                   + "Quadruples: %i\n\n" % diags_nbs['quadruples']
                   + "Quintuples and higher excitation levels: %i\n\n"
                   % diags_nbs['quintuples+'])


def print_cd_output(directory, diagrams):
    """Print a computer-readable file for automated frameworks.

    Args:
        directory (str): The path to the output directory.
        diagrams (list): All the MbptDiagrams.

    """
    cd_file = open(directory + '/CD_output.txt', 'w')
    conjug_file = open(directory + '/CD_conjug_pairs.list', 'w')
    for diag in diagrams:
        cd_file.write('config[%i] = %s\n' % (diag.tags[0] + 1, diag.cd_expr))
        if (diag.complex_conjugate != -1) \
                and (diag.complex_conjugate > diag.tags[0]):
            conjug_file.write("%i\t%i\n" % (diag.tags[0] + 1,
                                            diag.complex_conjugate + 1))
    cd_file.write('\n')
    cd_file.close()
    conjug_file.close()


def order_diagrams(diagrams):
    """Order the MBPT diagrams and return the number of diags for each type.

    Args:
        diagrams (list): The unordered MbptDiagrams.

    Returns:
        (tuple): First element are the ordered MbptDiagrams. Second element is
        the number of diagrams for each excitation level type.

    """
    singles = []
    doubles = []
    triples = []
    quadruples = []
    quintuples_and_higher = []

    for i_diag in xrange(len(diagrams)-1, -1, -1):
        if diagrams[i_diag].excitation_level == 1:
            singles.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level == 2:
            doubles.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level == 3:
            triples.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level == 4:
            quadruples.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level >= 5:
            quintuples_and_higher.append(diagrams[i_diag])
        else:
            print "Zero or negative excitation level!\n"
            exit()
        del diagrams[i_diag]

    diagrams = singles + doubles + triples + quadruples \
        + quintuples_and_higher

    for ind, diagram in enumerate(diagrams):
        diagram.tags[0] = ind

    attribute_conjugate(diagrams)

    diags_nb_per_type = {
        'nb_diags': len(diagrams),
        'singles': len(singles),
        'doubles': len(doubles),
        'triples': len(triples),
        'quadruples': len(quadruples),
        'quintuples+': len(quintuples_and_higher)
        }

    return diagrams, diags_nb_per_type


def attribute_conjugate(diagrams):
    """Attribute to each diagram its complex conjugate.

    The diagrams involved in conjugate pairs receive the tag associated to
    their partner in the ``complex_conjugate`` attribute.

    Args:
        diagrams (list): The topologically unique MbptDiagrams.
    """
    for idx, diag1 in enumerate(diagrams):
        if diag1.complex_conjugate == -1:
            for diag2 in diagrams[idx+1:]:
                if diag2.complex_conjugate == -1:
                    if diag1.is_complex_conjug_of(diag2):
                        diag1.complex_conjugate = diag2.tags[0]
                        diag2.complex_conjugate = diag1.tags[0]
                        break


def extract_cd_denom(start_graph, subgraph):
    """Extract the computer-readable denominator using the subgraph rule.

    Args:
        start_graph (NetworkX MultiDiGraph): The studied graph.
        subgraph (NetworkX MultiDiGraph): The subgaph for this particular
            factor.

    Returns:
        (str): The denominator factor associated to this subgraph.

    """
    denomin = "{" \
        + "".join("%s, "
                  % propa[3]['qp_state']
                  for propa
                  in start_graph.in_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2])) \
        + "".join("%s, "
                  % propa[3]['qp_state']
                  for propa
                  in start_graph.out_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2]))
    denomin = denomin.strip(', ') + '}'
    return denomin


class MbptDiagram(adg.diag.Diagram):
    """Describes a MBPT diagram with its related properties.

    Attributes:
        incidence (NumPy array): The incidence matrix of the graph.
        excitation_level (int): The single, double, etc., excitation character.
        complex_conjugate (int): The tag number of the diagram's complex
            conjugate. -1 is the graph has none.
        expr (str): The MBPT expression associated to the diagram.
        cd_expr (str): The expression associated to the diagram in a
            computer-readable format.

    """

    def __init__(self, mbpt_graph, tag_num):
        """Generate a MBPT diagram using the appropriate NetworkX graph.

        Args:
            mbpt_graph (NetworkX MultiDiGraph): The actual diagram.
            tag_num (int): The tag number associated to the graph.

        """
        adg.diag.Diagram.__init__(self, mbpt_graph)
        self.tags = [tag_num]
        # Beware of the sign convention !!!
        self.incidence = - nx.incidence_matrix(self.graph,
                                               oriented=True).todense()
        self.attribute_ph_labels()
        self.attribute_expression()
        self.excitation_level = self.calc_excitation()
        self.complex_conjugate = -1
        self.loops_number()

    def attribute_expression(self):
        """Initialize the expression associated to the diagram."""
        sign = "-" if (self.count_hole_lines()
                       - self.loops_number()) % 2 == 1 else ""
        eq_lines = np.array(self.incidence.transpose())
        neq_lines = np.asarray(list(i for i in set(map(tuple, eq_lines))))
        nedges_eq = 2**(len(eq_lines)-len(neq_lines))

        self.expr = sign \
            + ("\\dfrac{1}{%i}" % nedges_eq if nedges_eq != 1 else "") \
            + "\\sum{\\dfrac{%s}{%s}}\n" % (self.extract_numerator(),
                                            self.extract_denominator())
        self.cd_expr = "{%i, {%s}, {%s}};" % (nedges_eq,
                                              self.cd_numerator(),
                                              self.cd_denominator())

    def calc_excitation(self):
        """Return an integer coding for the excitation level of the diag.

        Returns:
            (int): The singles / doubles / etc. character of the graph.

        """
        max_excited_states = 0
        for row in xrange(1, self.graph.number_of_nodes()):
            nb_excited_states = 0
            for col in range(self.graph.number_of_edges()):
                if self.incidence[0:row, col].sum() == 1:
                    nb_excited_states += 1
                elif self.incidence[0:row, col].sum() == -1:
                    nb_excited_states += 1
            if nb_excited_states > max_excited_states \
                    and nb_excited_states != 4:
                max_excited_states = nb_excited_states
        return max_excited_states / 2 if max_excited_states != 0 else 2

    def count_hole_lines(self):
        """Return an integer for the number of hole lines in the graph.

        Returns:
            (int): The number of holes in the diagram.

        """
        return sum(1 for edge in self.graph.edges() if edge[0] > edge[1])

    def is_complex_conjug_of(self, test_diagram):
        """Return True if self and test_diagram are complex conjugate.

        Args:
            test_diagram (MbptDiagram): A diagram to compare with.

        Return:
            (bool): The complex conjugate status of the pair of diagrams.

        """
        # Check the adjacency mat against the anti-transposed one of test_diag
        return np.array_equal(self.adjacency_mat,
                              test_diagram.adjacency_mat[::-1, ::-1].T)

    def attribute_ph_labels(self):
        """Attribute the appropriate qp labels to the graph's propagators."""
        labels = list(string.ascii_lowercase)
        # Labelling needs to be shifted for higher orders
        if len(self.graph) < 6:
            h_labels = labels[0:8]
            p_labels = labels[8:]
        else:
            h_labels = labels[0:13]
            p_labels = labels[13:]
        for prop in self.graph.edges(keys=True, data=True):
            if prop[0] < prop[1]:
                prop[3]['qp_state'] = h_labels.pop(0)
            else:
                prop[3]['qp_state'] = p_labels.pop(0)

    def extract_denominator(self):
        """Return the denominator for a MBPT graph.

        Returns:
            (str): The denominator of the diagram.

        """
        denominator = ""
        graph = self.graph
        vertices = list(range(1, len(graph)))
        while len(vertices) >= 1:
            denominator += "%s\\ " % adg.diag.extract_denom(
                graph, graph.subgraph(vertices))
            del vertices[0]
        return denominator

    def cd_denominator(self):
        """Return the computer-readable denominator of the graph.

        Return:
            (str): The graph denominator tailored for automated frameworks.

        """
        denominator = ""
        graph = self.graph
        vertices = list(range(1, len(graph)))
        while len(vertices) >= 1:
            denominator += "%s, " % extract_cd_denom(graph,
                                                     graph.subgraph(vertices))
            del vertices[0]
        return denominator.strip(', ')

    def extract_numerator(self):
        """Return the numerator associated to a MBPT graph.

        Returns:
            (str): The numerator of the diagram.

        """
        graph = self.graph
        numerator = ""
        for vertex in graph:
            # First add the qp states corresponding to propagators going out
            numerator += "v_{" + "".join(
                prop[3]['qp_state']
                for prop in graph.out_edges(vertex, keys=True, data=True))
            # Add the qp states corresponding to propagators coming in
            numerator += "".join(
                prop[3]['qp_state']
                for prop in graph.in_edges(vertex, keys=True, data=True)) \
                + "} "
        return numerator

    def cd_numerator(self):
        """Return the computer-readable numerator.

        Returns:
            (str): The graph numerator tailored for automated frameworks.

        """
        graph = self.graph
        numerator = ""
        for vertex in graph:
            # Attribute the correct operator to each vertex
            numerator += "{"
            # First add the qp states corresponding to propagators going out
            numerator += ", ".join(
                prop[3]['qp_state']
                for prop in graph.out_edges(vertex, keys=True, data=True))
            numerator += ', '
            # Add the qp states corresponding to propagators coming in
            numerator += ", ".join(
                prop[3]['qp_state']
                for prop in graph.in_edges(vertex, keys=True, data=True))
            numerator = numerator.strip(', ') + "}, "
        return numerator.strip(', ')

    def loops_number(self):
        """Return the number of loops in the diagram as an integer.

        Returns:
            (int): The number of loops in the graph.

        """
        nb_loops = 0
        nb_checked_props = 0
        diag = self.graph
        nx.set_edge_attributes(diag, False, 'checked')
        while nb_checked_props < diag.number_of_edges():
            prop = list(edge for edge in diag.edges(keys=True, data=True)
                        if edge[3]['checked'] is False)[0]
            while prop[3]['checked'] is False:
                prop[3]['checked'] = True
                nb_checked_props += 1
                left_right_label = list(diag.in_edges(prop[1],
                                                      keys=True,
                                                      data=True)).index(prop)
                prop = list(diag.out_edges(prop[1],
                                           keys=True,
                                           data=True))[left_right_label]
            nb_loops += 1
        return nb_loops

    def write_section(self, result, commands, diags_nbs):
        """Write sections for MBPT result file.

        Args:
            result (file): The LaTeX output file to be written in.
            commands (dict): The flags associated with run management.
            diags_nbs (dict): A dict with the number of diagrams per
                excitation level type.

        """
        if self.tags[0] == 0 and diags_nbs['singles'] != 0:
            result.write("\\section{Singles}\n\n")
        elif self.tags[0] == diags_nbs['singles']:
            result.write("\\section{Doubles}\n\n")
        elif self.tags[0] == diags_nbs['singles'] + diags_nbs['doubles']:
            result.write("\\section{Triples}\n\n")
        elif self.tags[0] == (diags_nbs['singles'] + diags_nbs['doubles']
                              + diags_nbs['triples']):
            result.write("\\section{Quadruples}\n\n")
        elif self.tags[0] == (diags_nbs['singles'] + diags_nbs['doubles']
                              + diags_nbs['triples']
                              + diags_nbs['quadruples']):
            result.write("\\section{Quintuples and higher}\n\n")
        result.write("\\paragraph{Diagram %i:}\n" % (self.tags[0] + 1))
        if self.complex_conjugate >= 0:
            result.write("Complex conjugate diagram: %i\n"
                         % (self.complex_conjugate + 1))
        write_diag_exp(result, self)
