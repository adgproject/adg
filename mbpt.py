"""Module containg methods relative to MBPT, to be called by ADG."""

import copy
import itertools
import string
import numpy as np
import methods as mth


def seed(n):
    """Generate all 1-magic square of dimension n."""
    return [k for k in itertools.permutations(range(n), n)]


def diagram_generation(n):
    """Generate the diagrams for the MBPT case."""
    seeds = seed(n)
    all = [[[0 if i != j else 1 for i in range(n)] for j in k] for k in seeds]
    traceless = mth.no_trace(all)
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
    mbpt_diagrams = []
    for matrix in doubleUniq:
        mbpt_diagrams.append(np.array(matrix))
    return mbpt_diagrams


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


def write_diag_exp(latex_file, nedges_eq, phases, matrix_elements,
                   denominators):
    """Write the expression associated to a diagram in the LaTeX file."""
    diag_exp = "\\dfrac{1}{%i}" % nedges_eq + phases \
        + "\\sum{\\dfrac{" + matrix_elements + "}{" \
        + denominators + "}}\n"
    latex_file.write("\\begin{equation}\n")
    latex_file.write(diag_exp)
    latex_file.write("\\end{equation}\n")
