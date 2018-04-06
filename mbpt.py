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


def line_label_h(n):
    """Select appropriate label for hole line."""
    labels = list(string.ascii_lowercase)
    labels = labels[0:13]
    return labels[n]


def line_label_p(n):
    """Select appropriate label for particle line."""
    labels = list(string.ascii_lowercase)
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
    for diag in diagrams:
        CD_file.write('config[%i] = %s\n' % (diag.tags[0] + 1, diag.CD))
    CD_file.write('\n')
    CD_file.close()


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
        for diag2 in diagrams[idx+1:]:
            if diag1.is_complex_conjug_of(diag2):
                diag1.complex_conjugate = diag2.tags[0]
                diag2.complex_conjugate = diag1.tags[0]
                break


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
        self.attribute_expression()
        self.excitation_level = self.calc_excitation()
        self.complex_conjugate = -1

    def attribute_expression(self):
        """Initialize the expression associated to the diagram."""
        type_edg = self.edges_type()
        braket = ''
        braket_CD = ''
        nrow = self.graph.number_of_nodes()
        ncol = self.graph.number_of_edges()
        for row in xrange(nrow):
            ket = ''
            bra = ''
            ket_CD = ''
            bra_CD = ''
            for col in xrange(ncol):
                if self.incidence[row, col] == 1:
                    bra += line_label_h(col) if type_edg[col] == 'h' \
                        else line_label_p(col)
                    bra_CD += line_label_h(col) if type_edg[col] == 'h' \
                        else line_label_p(col)
                    bra_CD += ", "
                elif self.incidence[row, col] == -1:
                    ket += line_label_h(col) if type_edg[col] == 'h' \
                        else line_label_p(col)
                    ket_CD += line_label_h(col) if type_edg[col] == 'h' \
                        else line_label_p(col)
                    ket_CD += ", "
            braket += '\\braket{%s|H|%s}' % (bra, ket)
            braket_CD += '{%s, %s}, ' % (bra_CD.strip(', '),
                                         ket_CD.strip(', '))
        braket_CD = braket_CD.strip(', ')
        denom = ''
        denom_CD = ''
        for row in xrange(1, nrow):
            denom += '('
            denom_CD += '{'
            for col in range(ncol):
                if self.incidence[0:row, col].sum() == 1:
                    denom += '+E_' + line_label_h(col) if type_edg[col] == 'h'\
                        else '+E_' + line_label_p(col)
                    denom_CD += line_label_h(col) if type_edg[col] == 'h'\
                        else line_label_p(col)
                    denom_CD += ', '
                elif self.incidence[0:row, col].sum() == -1:
                    denom += '-E_' + line_label_h(col) if type_edg[col] == 'h'\
                        else '-E_' + line_label_p(col)
                    denom_CD += line_label_h(col) if type_edg[col] == 'h'\
                        else line_label_p(col)
                    denom_CD += ', '
            denom += ')'
            denom_CD = denom_CD.strip(', ') + '}, '
        denom_CD = denom_CD.strip(', ')
        denom = denom.replace('( +', '(').strip(' ')
        phases = '(-1)^{%i+l}' % self.count_hole_lines()
        eq_lines = np.array(self.incidence.transpose())
        neq_lines = np.asarray(list(i for i in set(map(tuple, eq_lines))))
        nedges_eq = 2**(len(eq_lines)-len(neq_lines))
        self.expr = "\\dfrac{1}{%i}%s" % (nedges_eq, phases) \
            + "\\sum{\\dfrac{%s}{%s}}\n" % (braket, denom)
        self.CD = "{%i, {%s}, {%s}};" % (nedges_eq, braket_CD, denom_CD)

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
        """Return True if the diagram and test_diagram are complex conjugate."""
        is_conjug = True
        if (self.excitation_level != test_diagram.excitation_level) \
                or (self.graph.number_of_edges()
                    != test_diagram.graph.number_of_edges()):
            is_conjug = False
            return is_conjug
        vertex_max = self.graph.number_of_nodes()
        for vertex_a in self.graph:
            for vertex_b in self.graph:
                if self.graph.number_of_edges(vertex_a, vertex_b) \
                        != test_diagram.graph.reverse().number_of_edges(
                                vertex_max-1-vertex_a, vertex_max-1-vertex_b):
                    is_conjug = False
                    break
        return is_conjug
