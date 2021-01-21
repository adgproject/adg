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
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], "BIMSRG", 2)
    diag = adg.bimsrg.BimsrgDiagram(graph, 0)

    diag_ref = ' C^{20}(11,20) = P(k_{1}/k_{2}) ' \
        + '\\sum_{p_{1}} A^{11}_{k_{1} p_{1}} B^{20}_{p_{1} k_{2}}'

    assert diag.expr == diag_ref


def test_write_section():
    """Test the write_section method."""

    diagram = np.array([[0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
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
        + ' C^{20}(11,20) = P(k_{1}/k_{2}) ' \
        + '\\sum_{p_{1}} A^{11}_{k_{1} p_{1}} B^{20}_{p_{1} k_{2}}\n' \
        + '\\end{equation}\n'

    with open('section_file.tex', 'w+') as result_file:
        diag.write_section(result_file, com, section_flags)

    with open('section_file.tex', "r") as result_file:
        assert result_file.read() == section_ref
    os.unlink('section_file.tex')
