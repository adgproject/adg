"""Tests for the adg.bmbpt module."""

import argparse
import os
import numpy as np
import networkx as nx
import adg.diag
import adg.bimsrg


def test_attribute_expressions():
    """Test the attribute_expressions method."""
    # Test for a simple diagram
    diagram = np.array([[0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BIMSRG", 2)
    diag = adg.bimsrg.BimsrgDiagram(graph, 0)

    diag_ref = ' C^{20}_{k_{1}k_{2}}(11,20) = P(k_{1}/k_{2}) ' \
        + '\\sum_{p_{1}} A^{11}_{k_{1} p_{1}} B^{20}_{p_{1} k_{2}}'

    assert diag.expr == diag_ref


def test_write_header():
    """Test the write_header function."""
    com = argparse.Namespace()
    com.theory, com.with_3NF, com.canonical = 'BIMSRG', False, False
    com.order = (2, 2, 2)

    diags_nbs = {'nb_diags': 82, 1: 10, 2: 72}

    header_ref = '$C=\\left[A,B\\right]$ with $N_A = 2$, $N_B = 2$ and $N_C = 2$\n\n' \
        + '$d_\\mathrm{max} \\equiv \\mathrm{max}(d_A, d_B, d_C)\n\n' \
        + '\\par \\vspace{\\baselineskip}\nValid diagrams: 82\n\n' \
        + '$d_\\mathrm{max} = 1$ diagrams: 10\n\n' \
        + '$d_\\mathrm{max} = 2$ diagrams: 72\n\n'

    with open('header.tex', 'w') as header_file:
        adg.bimsrg.write_header(header_file, com, diags_nbs)

    with open('header.tex', 'r') as header_file:
        assert header_file.read() == header_ref
    os.unlink('header.tex')


def test_write_permutator_section():
    """Test the write_permutator_section function."""
    com = argparse.Namespace()
    com.theory, com.with_3NF, com.canonical = 'BIMSRG', False, False
    com.order = (2, 2, 2)

    perm_ref = '\\section{Permutators definitions}\n\n\\begin{align*}\n' \
        + 'P(k_{1}/k_{2}) &= 1 - P_{k_{1} k_{2}} \\\\\n' \
        + 'P(k_{1}/k_{2}k_{3}) &= 1 - P_{k_{1} k_{2}} - P_{k_{1} k_{3}} \\\\\n' \
        + 'P(k_{1}/k_{2}k_{3}k_{4}) &= 1 - P_{k_{1} k_{2}} - P_{k_{1} k_{3}} - P_{k_{1} k_{4}} \\\\\n' \
        + 'P(k_{1}k_{2}/k_{3}k_{4}) &= 1 - P_{k_{1} k_{3}} - P_{k_{1} k_{4}} ' \
        + '- P_{k_{2} k_{3}} - P_{k_{2} k_{4}} + P_{k_{1} k_{3}} P_{k_{2} k_{4}} ' \
        + '+ P_{k_{2} k_{3}} P_{k_{1} k_{4}} \\\\\n\\end{align*}\n'

    with open('perm.tex', 'w') as perm_file:
        adg.bimsrg.write_permutator_section(perm_file, com)

    with open('perm.tex', 'r') as perm_file:
        assert perm_file.read() == perm_ref
    os.unlink('perm.tex')


def test_write_section():
    """Test the write_section method."""

    diagram = np.array([[0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BIMSRG", 2)
    diag = adg.bimsrg.BimsrgDiagram(graph, 0)

    com = argparse.Namespace()
    com.theory, com.with_3NF, com.canonical = 'BIMSRG', False, False
    com.order = (2, 2, 2)

    section_flags = {'new_op_struct': [0], 1: -1, 2: 0}

    section_ref = '\\section{$d_\\mathrm{max} = 2$}\n\n' \
        + '\\subsection{$C^{20}$}\n\n' \
        + '\\paragraph{Diagram 1 ($+AB$):}\n' \
        + '\\begin{equation}\n' \
        + ' C^{20}_{k_{1}k_{2}}(11,20) = P(k_{1}/k_{2}) ' \
        + '\\sum_{p_{1}} A^{11}_{k_{1} p_{1}} B^{20}_{p_{1} k_{2}}\n' \
        + '\\end{equation}\n'

    with open('section_file.tex', 'w+') as result_file:
        diag.write_section(result_file, com, section_flags)

    with open('section_file.tex', "r") as result_file:
        assert result_file.read() == section_ref
    os.unlink('section_file.tex')
