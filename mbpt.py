"""Module containg methods relative to MBPT, to be called by ADG."""

import copy
import itertools
import string
import numpy as np
import networkx as nx
import generic_diag as gen


def diagrams_generation(order):
    """Generate the diagrams for the MBPT case."""
    # Generate all 1-magic square of dimension order
    seeds = [k for k in itertools.permutations(range(order), order)]
    all_matrices = [[[0 if i != j else 1 for i in range(order)]
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
    double_uniq = []
    for matrix in double:
        if matrix not in double_uniq:
            double_uniq.append(matrix)
    double_uniq.sort(reverse=True)
    return [np.array(matrix) for matrix in double_uniq]


def write_diag_exp(latex_file, mbpt_diag):
    """Write the expression associated to a diagram in the LaTeX file."""
    latex_file.write("\\begin{equation}\n")
    latex_file.write(mbpt_diag.expr)
    latex_file.write("\\end{equation}\n")


def write_header(tex_file, diags_nbs):
    """Write tha appropriate header for the LaTeX file for MBPT diagrams."""
    tex_file.write("Valid diagrams: %i\n\n" % diags_nbs['nb_diags']
                   + "Singles: %i\n\n" % diags_nbs['singles']
                   + "Doubles: %i\n\n" % diags_nbs['doubles']
                   + "Triples: %i\n\n" % diags_nbs['triples']
                   + "Quadruples: %i\n\n" % diags_nbs['quadruples']
                   + "Quintuples and higher excitation levels: %i\n\n"
                   % diags_nbs['quintuples+'])


def write_section(result, diag_index, diags_nbs):
    """Write sections for MBPT result file."""
    if diag_index == 0 and diags_nbs['singles'] != 0:
        result.write("\\section{Singles}\n\n")
    elif diag_index == diags_nbs['singles']:
        result.write("\\section{Doubles}\n\n")
    elif diag_index == diags_nbs['singles'] + diags_nbs['doubles']:
        result.write("\\section{Triples}\n\n")
    elif diag_index == (diags_nbs['singles'] + diags_nbs['doubles']
                        + diags_nbs['triples']):
        result.write("\\section{Quadruples}\n\n")
    elif diag_index == (diags_nbs['singles'] + diags_nbs['doubles']
                        + diags_nbs['triples'] + diags_nbs['quadruples']):
        result.write("\\section{Quintuples and higher}\n\n")


def print_cd_output(directory, diagrams):
    """Print a computer-readable file for C. Drischler's framework."""
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
    with open(directory+"/CD_adj_matrices.list", "w") as mat_file:
        for idx, diagram in enumerate(diagrams):
            mat_file.write("Diagram n: %i\n" % (idx + 1))
            np.savetxt(mat_file, diagram.adjacency_mat, fmt='%d')
            mat_file.write("\n")


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
    """Attribute to each diagram its complex conjugate."""
    for idx, diag1 in enumerate(diagrams):
        if diag1.complex_conjugate == -1:
            for diag2 in diagrams[idx+1:]:
                if diag2.complex_conjugate == -1:
                    if diag1.is_complex_conjug_of(diag2):
                        diag1.complex_conjugate = diag2.tags[0]
                        diag2.complex_conjugate = diag1.tags[0]
                        break


def extract_cd_denom(start_graph, subgraph):
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
        self.cd_expr = "{%i, {%s}, {%s}};" % (nedges_eq,
                                              self.cd_numerator(),
                                              self.cd_denominator())

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

    def cd_denominator(self):
        """Return the CD-formatted denominator of the graph."""
        denominator = ""
        graph = self.graph
        for vertex_i in range(1, len(graph)):
            substack = [vertex_j
                        for vertex_j in graph
                        if vertex_j >= vertex_i]
            denominator += "%s, " % extract_cd_denom(graph,
                                                     graph.subgraph(substack))
        return denominator.strip(', ')

    def extract_numerator(self):
        """Return the numerator associated to a MBPT graph."""
        graph = self.graph
        numerator = ""
        for vertex in graph:
            # First add the qp states corresponding to propagators going out
            numerator += "v_{" + "".join(
                graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                for prop in graph.out_edges(vertex, keys=True))
            # Add the qp states corresponding to propagators coming in
            numerator += "".join(
                graph.adj[prop[0]][prop[1]][prop[2]]['qp_state']
                for prop in graph.in_edges(vertex, keys=True)) \
                + "} "
        return numerator

    def cd_numerator(self):
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
