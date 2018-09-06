"""Routines and class for Projected Bogoliubov MBPT diagrams."""

import copy
import itertools
import numpy as np
import networkx as nx
import adg.bmbpt
import adg.tsd
import adg.diag


def generate_anomalous_diags(graph, nbody_max):
    """Generate PBMBPT graphs with anomalous lines.

    Args:
        graph (NetworkX MultiDiGraph): The graph to generate children from.
        nbody_max (int): The maximal n-body character of a graph vertex.

    Returns:
        (list): The anomalous graphs generated.

    """
    anom_graphs = [graph]
    vertices = [vert for vert in graph if not graph.node[vert]['operator']]

    # Loop to generate anomalous propagators
    for iter_graph in anom_graphs:
        for vertex1 in vertices:
            edges = iter_graph.edges(vertex1, keys=True, data='anomalous')
            for vertex2 in vertices[vertex1 + 1:]:
                for edge in edges:
                    if (edge[1], edge[3]) == (vertex2, False):
                        new_graph = copy.deepcopy(iter_graph)
                        new_graph[vertex1][vertex2][edge[2]]['anomalous'] = True
                        anom_graphs.append(new_graph)
                        break

    # Loop to generate self-contractions
    for iter_graph in anom_graphs:
        for vertex in vertices:
            if iter_graph.degree(vertex) / 2 < nbody_max:
                new_graph = copy.deepcopy(iter_graph)
                new_graph.add_edge(vertex, vertex, anomalous=True)
                anom_graphs.append(new_graph)

    return anom_graphs


class ProjectedBmbptDiagram(adg.bmbpt.BmbptFeynmanDiagram):
    """Describes a PBMBPT diagram with its related properties.

    Attributes:
        two_or_three_body (int): The 2 or 3-body characted of the vertices.
        time_tag (int): The tag number associated to the diagram's
            associated TSD.
        tsd_is_tree (bool): The tree or non-tree character of the
            associated TSD.
        feynman_exp (str): The Feynman expression associated to the diagram.
        diag_exp (str): The Goldstone expression associated to the diagram.
        vert_exp (list): The expression associated to the vertices.
        hf_type (str): The Hartree-Fock, non-Hartree-Fock or Hartree-Fock for
            the energy operator only character of the graph.

    """

    def __init__(self, bmbpt_diag):
        """Generate a PBMBPT diagram by copying a BMBPT one.

        Args:
            bmbpt_diag (adg.bmbpt.BmbptFeynmanDiagram)): The graph of interest.

        """
        adg.bmbpt.BmbptFeynmanDiagram.__init__(self,
                                               bmbpt_diag.graph,
                                               bmbpt_diag.tags[0])
        self.tags = [bmbpt_diag.tags[0]]
        self.time_tag = bmbpt_diag.time_tag
        self.tsd_is_tree = bmbpt_diag.tsd_is_tree
        self.vert_exp = bmbpt_diag.vert_exp
        for edge in self.graph.edges(keys=True, data=True):
            edge[3]['anomalous'] = False
