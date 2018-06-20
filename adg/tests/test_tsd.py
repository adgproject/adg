"""Tests for the adg.tsd module."""

import numpy as np
import networkx as nx
import adg.tsd
import adg.bmbpt
import adg.diag


def test_time_structure_graph():
    """Test the time_structure_graph routine."""
    # Test a simple case
    diagram = np.array([[0, 2, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT')

    assert list(adg.tsd.time_structure_graph(graph).edges()) == [(0, 1),
                                                                 (1, 2)]

    # Test the case where not all vertices are link to the bottom vertex
    diagram = np.array([[0, 0, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT')

    assert list(adg.tsd.time_structure_graph(graph).edges()) == [(0, 1),
                                                                 (1, 2)]

    # Test for a cycle topography
    diagram = np.array([[0, 2, 2, 0],
                        [0, 0, 0, 2],
                        [0, 0, 0, 2],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT')

    assert list(adg.tsd.time_structure_graph(graph).edges()) == [(0, 1),
                                                                 (0, 2),
                                                                 (1, 3),
                                                                 (2, 3)]


def test_find_cycle():
    """Test for the find_cycle routine."""
    # Test the simplest possible cycle
    diagram = np.array([[0, 1, 1, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)

    assert adg.tsd.find_cycle(graph) == (0, 3)

    # Test a case where false positive could occur
    diagram = np.array([[0, 1, 0, 0, 0, 1],
                        [0, 0, 1, 1, 0, 0],
                        [0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 1, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)

    assert adg.tsd.find_cycle(graph) == (1, 4)
