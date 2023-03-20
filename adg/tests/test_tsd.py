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
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert list(adg.tsd.time_structure_graph(diag).edges()) == [(0, 1),
                                                                (1, 2)]

    # Test the case where not all vertices are link to the bottom vertex
    diagram = np.array([[0, 0, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert list(adg.tsd.time_structure_graph(diag).edges()) == [(0, 1),
                                                                (1, 2)]

    # Test for a cycle topography
    diagram = np.array([[0, 2, 2, 0],
                        [0, 0, 0, 2],
                        [0, 0, 0, 2],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert list(adg.tsd.time_structure_graph(diag).edges()) == [(0, 1),
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
    graph = nx.from_numpy_array(diagram,
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
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)

    assert adg.tsd.find_cycle(graph) == (1, 4)


def test_disentangle_cycle():
    """Test the disentangle_cycle routine."""
    # Test for the simplest cycle TSD
    diagram = np.array([[0, 1, 1, 0],
                        [0, 0, 0, 1],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    result = adg.tsd.disentangle_cycle(graph, (0, 3))

    assert len(result) == 2
    assert list(result[0].edges()) == [(0, 1), (1, 2), (2, 3)]
    assert list(result[1].edges()) == [(0, 2), (1, 3), (2, 1)]


def test_treat_tsds():
    """Test for the treat_tsds routine."""
    # Setup for the test with one non-tree and two isomorphic trees
    diagrams = [
        np.array([[0, 1, 1, 0], [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 0]]),
        np.array([[0, 2, 2, 0], [0, 0, 0, 0], [0, 0, 0, 2], [0, 0, 0, 0]]),
        np.array([[0, 2, 2, 0], [0, 0, 0, 2], [0, 0, 0, 0], [0, 0, 0, 0]])
    ]
    graphs = [nx.from_numpy_array(diagram,
                                  create_using=nx.MultiDiGraph(),
                                  parallel_edges=True)
              for diagram in diagrams]
    adg.diag.label_vertices(graphs, "BMBPT", -1)
    diags = [adg.bmbpt.BmbptFeynmanDiagram(graph, ind)
             for ind, graph in enumerate(graphs)]
    for diag in diags:
        diag.attribute_qp_labels()
    tsds = [adg.tsd.TimeStructureDiagram(diag) for diag in diags]

    result_tsds, nb_trees = adg.tsd.treat_tsds(tsds)

    assert nb_trees == 1
    assert len(result_tsds) == 2
    assert result_tsds[0].is_tree
    assert not result_tsds[1].is_tree
    assert len(result_tsds[1].equivalent_trees) == 2


def test___init__():
    """Test the creation of TSD objects."""
    # Test for a tree TSD
    diagram = np.array([[0, 2, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    assert tsd.is_tree
    assert tsd.resum == 1

    # Test for the simplest non-tree TSD
    diagram = np.array([[0, 2, 2, 0],
                        [0, 0, 0, 2],
                        [0, 0, 0, 2],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    assert not tsd.is_tree
    assert tsd.resum == 0


def test_treat_cycles():
    """Test for the treat_cycles method."""
    # Test the simplest non-tree tsd
    diagram = np.array([[0, 2, 2, 0],
                        [0, 0, 0, 2],
                        [0, 0, 0, 2],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    equivalent_trees = tsd.treat_cycles()

    assert len(equivalent_trees) == 2
    assert list(equivalent_trees[1].edges()) == [(0, 1), (1, 2), (2, 3)]
    assert list(equivalent_trees[0].edges()) == [(0, 2), (1, 3), (2, 1)]

    # Test for a case with a longer branch
    diagram = np.array([[0, 2, 2, 0, 0],
                        [0, 0, 0, 2, 0],
                        [0, 0, 0, 0, 2],
                        [0, 0, 0, 0, 2],
                        [0, 0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    equivalent_trees = tsd.treat_cycles()

    assert len(equivalent_trees) == 3


def test_resummation_power():
    """Test for the resummation_power method."""
    # Test for a linear tree TSD
    diagram = np.array([[0, 2, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    assert tsd.resummation_power() == 1

    # Test for a non-linear tree TSD
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    assert tsd.resummation_power() == 2
