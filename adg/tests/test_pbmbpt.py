"""Tests for the adg.pbmbpt module."""

import numpy as np
import networkx as nx
import adg.pbmbpt
import adg.diag


def test_generate_anomalous_diags():
    """Test for the generate_anomalous_diags routine."""
    # Generate a simple BMBPT diagram
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT")
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()

    for edge in diag.graph.edges(keys=True, data='anomalous'):
        diag.graph[edge[0]][edge[1]][edge[2]]['anomalous'] = False

    two_body_diags = adg.pbmbpt.generate_anomalous_diags(diag.graph, 2)
    assert len(two_body_diags) == 4

    three_body_diags = adg.pbmbpt.generate_anomalous_diags(diag.graph, 3)
    assert len(three_body_diags) == 9

    diagram = np.array([[0, 2, 0], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_matrix(diagram,
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT")
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    for edge in diag.graph.edges(keys=True, data='anomalous'):
        diag.graph[edge[0]][edge[1]][edge[2]]['anomalous'] = False

    two_body_diags = adg.pbmbpt.generate_anomalous_diags(diag.graph, 2)
    assert len(two_body_diags) == 6

    three_body_diags = adg.pbmbpt.generate_anomalous_diags(diag.graph, 3)
    assert len(three_body_diags) == 18
