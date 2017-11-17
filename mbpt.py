"""Module containg methods relative to MBPT, to be called by ADG."""

import copy
import itertools
import string
import numpy as np
import networkx as nx
import methods as mth


def seed(n):
    """Generate all 1-magic square of dimension n."""
    return [k for k in itertools.permutations(range(n), n)]


def diagram_generation(n):
    """Generate the diagrams for the MBPT case."""
    seeds = seed(n)
    all_matrices = [[[0 if i != j else 1 for i in range(n)]
                     for j in k]
                    for k in seeds]
    traceless = mth.no_trace(all_matrices)
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
    labels = labels[0:15]
    return labels[n]


def line_label_p(n):
    """Select appropriate label for particle line."""
    labels = list(string.ascii_lowercase)
    labels = labels[15:-1]
    return labels[n]


def write_diag_exp(latex_file, mbpt_diag):
    """Write the expression associated to a diagram in the LaTeX file."""
    latex_file.write("\\begin{equation}\n")
    latex_file.write(mbpt_diag.expr)
    latex_file.write("\\end{equation}\n")


class MbptDiagram(mth.Diagram):
    """Describes a MBPT diagram with its related properties."""

    def __init__(self, mbpt_graph, tag_num):
        """Generate a MBPT diagram using the appropriate NetworkX graph."""
        mth.Diagram.__init__(self, mbpt_graph)
        self.tags = [tag_num]
        self.attribute_expression()

    def attribute_expression(self):
        """Initialize the expression associated to the diagram."""
        type_edg = []
        braket = ''
        # Beware of the sign convention !!!
        incidence = - nx.incidence_matrix(self.graph, oriented=True).todense()
        nrow = self.graph.number_of_nodes()
        ncol = self.graph.number_of_edges()
        n_holes = 0
        diffcols = set()
        for col in range(ncol):
            flat = list(incidence[:, col].A1)
            if flat.index(1) < flat.index(-1):
                n_holes += 1
                type_edg.append('h')
            else:
                type_edg.append('p')
            diffcols.add(repr(flat))

        for row in xrange(nrow):
            ket = ''
            bra = ''
            for col in xrange(ncol):
                if (incidence[row, col] == 1):
                    bra += line_label_h(col) if type_edg[col] == 'h' \
                        else line_label_p(col)
                elif (incidence[row, col] == -1):
                    ket += line_label_h(col) if type_edg[col] == 'h' \
                        else line_label_p(col)
            braket += '\\braket{%s|H|%s}' % (bra, ket)
        denom = ''
        for row in xrange(1, nrow):
            denom += '('
            for col in range(ncol):
                if incidence[0:row, col].sum() == 1:
                    denom += '+E_' + line_label_h(col) if type_edg[col] == 'h'\
                        else '+E_' + line_label_p(col)
                elif incidence[0:row, col].sum() == -1:
                    denom += '-E_' + line_label_h(col) if type_edg[col] == 'h'\
                        else '-E_' + line_label_p(col)
            denom += ')'
        denom = denom.replace('( +', '(').strip(' ')
        phases = '(-1)^{%i+l}' % n_holes
        eq_lines = np.array(incidence.transpose())
        neq_lines = np.asarray(list(i for i in set(map(tuple, eq_lines))))
        nedges_eq = 2**(len(eq_lines)-len(neq_lines))
        self.expr = "\\dfrac{1}{%i}%s" % (nedges_eq, phases) \
            + "\\sum{\\dfrac{%s}{%s}}\n" % (braket, denom)
