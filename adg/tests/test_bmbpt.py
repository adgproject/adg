"""Tests for the adg.bmbpt module."""

import numpy as np
import networkx as nx
import adg.bmbpt
import adg.diag


def test_vertex_expression():
    """Test the vertex_expression method."""
    # Set up the test case
    diagram = np.array([[0, 2, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    # attribute_qp_labels is used and thus tested in the process
    diag.attribute_qp_labels()

    assert diag.vertex_expression(0) == '\\epsilon^{k_{1}k_{2}k_{3}k_{4}}_{}'
    assert diag.vertex_expression(1) == '\\epsilon^{k_{5}k_{6}}_{k_{1}k_{2}}'
    assert diag.vertex_expression(2) == '\\epsilon^{}_{k_{3}k_{4}k_{5}k_{6}}'


def test_vertex_exchange_sym_factor():
    """Test the vertex_exchange_sym_factor method."""
    # Test case with symmetry factor
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT')
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.vertex_exchange_sym_factor() == '2'

    # Test case without symmetry factor
    diagram = np.array([[0, 2, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT')
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.vertex_exchange_sym_factor() == ''


def test_has_crossing_sign():
    """Test for the has_crossing_sign method."""
    # Test case without crossing sign
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert not diag.has_crossing_sign()

    # Test case with a crossing sign
    diagram = np.array([[0, 2, 1, 1],
                        [0, 0, 1, 1],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.has_crossing_sign()


def test_multiplicity_symmetry_factor():
    """Test the multiplicity_symmetry_factor method."""
    # Test a simple case with twice two lines
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.multiplicity_symmetry_factor() == '(2!)^2'

    # Test a case with once three lines
    diagram = np.array([[0, 3, 1], [0, 0, 1], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.multiplicity_symmetry_factor() == '(3!)'
