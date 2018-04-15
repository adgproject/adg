"""Module containg methods relative to MBPT, to be called by ADG."""

import copy
import itertools
import string
import numpy as np
import networkx as nx
import general_routines as gen


def seed(n):
    """Generate all 1-magic square of dimension n."""
    return [k for k in itertools.permutations(range(n), n)]


def diagram_generation(n):
    """Generate the diagrams for the MBPT case."""
    seeds = seed(n)
    all_matrices = [[[0 if i != j else 1 for i in range(n)]
                     for j in k]
                    for k in seeds]
    traceless = gen.no_trace(all_matrices)
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
    return [np.array(matrix) for matrix in doubleUniq]


def line_label_h(n, order):
    """Select appropriate label for hole line."""
    labels = list(string.ascii_lowercase)
    if order < 5:
        labels = labels[0:6]
    else:
        labels = labels[0:13]
    return labels[n]


def line_label_p(n, order):
    """Select appropriate label for particle line."""
    labels = list(string.ascii_lowercase)
    if order < 5:
        labels = labels[6:-1]
    else:
        labels = labels[13:-1]
    return labels[n]


def write_diag_exp(latex_file, mbpt_diag):
    """Write the expression associated to a diagram in the LaTeX file."""
    latex_file.write("\\begin{equation}\n")
    latex_file.write(mbpt_diag.expr)
    latex_file.write("\\end{equation}\n")


def write_MBPT_header(tex_file, numdiag, nb_singles, nb_doubles, nb_triples,
                      nb_quadruples, nb_quintuples_and_higher):
    """Write tha appropriate header for the LaTeX file for MBPT diagrams."""
    tex_file.write("Valid diagrams: %i\n\n" % numdiag
                   + "Singles: %i\n\n" % nb_singles
                   + "Doubles: %i\n\n" % nb_doubles
                   + "Triples: %i\n\n" % nb_triples
                   + "Quadruples: %i\n\n" % nb_quadruples
                   + "Quintuples and higher excitation levels: %i\n\n"
                   % nb_quintuples_and_higher)


def write_MBPT_section(result, diag_index, nb_singles, nb_doubles, nb_triples,
                       nb_quadruples, nb_quintuples_and_higher):
    """Write sections for MBPT result file."""
    if diag_index == 0 and nb_singles != 0:
        result.write("\\section{Singles}\n\n")
    elif diag_index == nb_singles:
        result.write("\\section{Doubles}\n\n")
    elif diag_index == nb_singles + nb_doubles:
        result.write("\\section{Triples}\n\n")
    elif diag_index == nb_singles + nb_doubles + nb_triples:
        result.write("\\section{Quadruples}\n\n")
    elif diag_index == nb_singles + nb_doubles + nb_triples + nb_quadruples:
        result.write("\\section{Quintuples and higher}\n\n")


def print_CD_output(directory, diagrams):
    """Print a computer-readable file for C. Drischler's framework."""
    CD_file = open(directory + '/CD_output.txt', 'w')
    conjug_file = open(directory + '/CD_conjug_pairs.list', 'w')
    for diag in diagrams:
        CD_file.write('config[%i] = %s\n' % (diag.tags[0] + 1, diag.CD))
        if (diag.complex_conjugate != -1) \
                and (diag.complex_conjugate > diag.tags[0]):
            conjug_file.write("%i\t%i\n" % (diag.tags[0] + 1,
                                            diag.complex_conjugate + 1))
    CD_file.write('\n')
    CD_file.close()
    conjug_file.close()
    with open(directory+"/CD_adj_matrices.list", "w") as f:
        for idx, diagram in enumerate(diagrams):
            f.write("Diagram n: %i\n" % (idx + 1))
            np.savetxt(f, diagram.adjacency_mat, fmt='%d')
            f.write("\n")


def order_diagrams(diagrams):
    """Order the MBPT diagrams and return the number of diags for each type."""
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

    return diagrams, len(singles), len(doubles), len(triples), \
        len(quadruples), len(quintuples_and_higher)


def attribute_conjugate(diagrams):
    """Attribute to each diagram its complex conjugate."""
    for idx, diag1 in enumerate(diagrams):
        if diag1.complex_conjugate == -1:
            for diag2 in diagrams[idx+1:]:
                if diag2.complex_conjugate == -1:
                    if diag1.is_complex_conjug_of(diag2):
                        diag1.complex_conjugate = diag2.tags[0]
                        diag2.complex_conjugate = diag1.tags[0]
                        break


def extract_CD_denom(start_graph, subgraph):
    """Extract the appropriate CD denominator using the subgraph rule."""
    denomin = "{" \
        + "".join("%s, "
                  % start_graph.adj[propa[0]][propa[1]][propa[2]]['qp_state']
                  for propa
                  in start_graph.in_edges(subgraph, keys=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2])) \
        + "".join("%s, "
                  % start_graph.adj[propa[0]][propa[1]][propa[2]]['qp_state']
                  for propa
                  in start_graph.out_edges(subgraph, keys=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2]))
    denomin = denomin.strip(', ') + '}'
    return denomin


class MbptDiagram(gen.Diagram):
    """Describes a MBPT diagram with its related properties."""

    def __init__(self, mbpt_graph, tag_num):
        """Generate a MBPT diagram using the appropriate NetworkX graph."""
        gen.Diagram.__init__(self, mbpt_graph)
        self.type = 'MBPT'
        self.tags = [tag_num]
        # Beware of the sign convention !!!
        self.incidence = - nx.incidence_matrix(self.graph,
                                               oriented=True).todense()
        self.attribute_ph_labels()
        self.attribute_expression()
        self.excitation_level = self.calc_excitation()
        self.complex_conjugate = -1
        self.adjacency_mat = nx.to_numpy_matrix(self.graph, dtype=int)

    def attribute_expression(self):
        """Initialize the expression associated to the diagram."""
        phases = '(-1)^{%i+l}' % self.count_hole_lines()
        eq_lines = np.array(self.incidence.transpose())
        neq_lines = np.asarray(list(i for i in set(map(tuple, eq_lines))))
        nedges_eq = 2**(len(eq_lines)-len(neq_lines))

        self.expr = "\\dfrac{1}{%i}%s" % (nedges_eq, phases) \
            + "\\sum{\\dfrac{%s}{%s}}\n" % (self.extract_numerator(),
                                            self.extract_denominator())
        self.CD = "{%i, {%s}, {%s}};" % (nedges_eq,
                                         self.CD_numerator(),
                                         self.CD_denominator())

    def calc_excitation(self):
        """Return an integer coding for the excitation level of the diag."""
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
        """Return an integer for the number of hole lines in the graph."""
        nb_holes = 0
        for edge in self.graph.edges():
            if edge[0] > edge[1]:
                nb_holes += 1
        return nb_holes

    def edges_type(self):
        """Return a list with the types of the graph edges."""
        ncol = self.graph.number_of_edges()
        type_edg = []
        for col in range(ncol):
            flat = list(self.incidence[:, col].A1)
            if flat.index(1) < flat.index(-1):
                type_edg.append('h')
            else:
                type_edg.append('p')
        return type_edg

    def is_complex_conjug_of(self, test_diagram):
        """Return True if self and test_diagram are complex conjugate."""
        is_conjug = True
        # Check the adjacency mat against the anti-transposed one of test_diag
        if not np.array_equal(self.adjacency_mat,
                              test_diagram.adjacency_mat[::-1, ::-1].T):
            is_conjug = False
        return is_conjug

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
        for prop in self.graph.edges(keys=True):
            if prop[0] < prop[1]:
                self.graph.adj[prop[0]][prop[1]][prop[2]]['qp_state'] \
                    = h_labels.pop(0)
            else:
                self.graph.adj[prop[0]][prop[1]][prop[2]]['qp_state'] \
                    = p_labels.pop(0)

    def extract_denominator(self):
        """Return the denominator for a MBPT graph."""
        denominator = ""
        graph = self.graph
        for vertex_i in range(1, len(graph)):
            stack = [vertex_j for vertex_j in graph if vertex_j >= vertex_i]
            denominator += "%s\\ " % gen.extract_denom(graph,
                                                       graph.subgraph(stack))
        return denominator

    def CD_denominator(self):
        """Return the CD-formatted denominator of the graph."""
        denominator = ""
        graph = self.graph
        for vertex_i in range(1, len(graph)):
            substack = [vertex_j
                        for vertex_j in graph
                        if vertex_j >= vertex_i]
            denominator += "%s, " % extract_CD_denom(graph,
                                                     graph.subgraph(substack))
        return denominator.strip(', ')

    def extract_numerator(self):
        """Return the numerator associated to a MBPT graph."""
        graph = self.graph
        numerator = ""
        for vertex in graph:
            # Attribute the correct operator to each vertex
            numerator += "O_{" if graph.node[vertex]['operator'] \
                               else "\\Omega_{"
            # First add the qp states corresponding to propagators going out
            numerator += "".join(
                graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                for prop in graph.out_edges(vertex, keys=True))
            # Add the qp states corresponding to propagators coming in
            numerator += "".join(
                graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                for prop in graph.in_edges(vertex, keys=True)) \
                + "} "
        return numerator

    def CD_numerator(self):
        """Return the numerator under CD format."""
        graph = self.graph
        numerator = ""
        for vertex in graph:
            # Attribute the correct operator to each vertex
            numerator += "{"
            # First add the qp states corresponding to propagators going out
            numerator += ", ".join(
                graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                for prop in graph.out_edges(vertex, keys=True))
            numerator += ', '
            # Add the qp states corresponding to propagators coming in
            numerator += ", ".join(
                graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                for prop in graph.in_edges(vertex, keys=True))
            numerator = numerator.strip(', ') + "}, "
        return numerator.strip(', ')
