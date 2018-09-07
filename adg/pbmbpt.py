"""Routines and class for Projected Bogoliubov MBPT diagrams."""

import copy
import itertools
import adg.bmbpt
import adg.tsd
import adg.diag


def generate_anomalous_diags(graph, nbody_max):
    """Generate PBMBPT graphs with anomalous lines, with some redundancy.

    Args:
        graph (NetworkX MultiDiGraph): The graph to generate children from.
        nbody_max (int): The maximal n-body character of a graph vertex.

    Returns:
        (list): The anomalous graphs generated.

    """
    anom_graphs = [graph]
    vertices = [vert for vert in graph if not graph.node[vert]['operator']]

    # Turn normal propagators in anomalous ones
    tweakable_edges = []
    for vert1 in vertices:
        for vert2 in vertices[vert1:]:
            for _ in range(graph.number_of_edges(vert1, vert2)):
                tweakable_edges.append((vert1, vert2))
    for comb in generate_combinations(tweakable_edges):
        new_graph = copy.deepcopy(graph)
        for edge in comb:
            key = sum(1 for prop
                      in new_graph.edges(edge[0], edge[1], keys=True)
                      if new_graph[prop[0]][prop[1]][prop[2]]['anomalous'])
            new_graph[edge[0]][edge[1]][key]['anomalous'] = True
        anom_graphs.append(new_graph)

    # Loop to generate self-contractions
    for idx in xrange(len(anom_graphs)-1, -1, -1):
        iter_graph = anom_graphs[idx]
        test_vertices = []
        for vert in vertices:
            for _ in range(nbody_max - (iter_graph.degree(vert) / 2)):
                test_vertices.append(vert)
        for comb in generate_combinations(test_vertices):
            new_graph = copy.deepcopy(iter_graph)
            for vert in comb:
                new_graph.add_edge(vert, vert, anomalous=True)
            anom_graphs.append(new_graph)

    return anom_graphs


def generate_combinations(iter_list):
    """Generate all possible combinations of length 1 to total.

    Attributes:
        iter_list (list): A list of iterable objects.

    Returns:
        (list): A list with all the possible combinations of all lengths.

    """
    combinations = []
    for i in range(1, len(iter_list) + 1):
        combinations += [k for k in itertools.combinations(iter_list, i)]
    # Remove duplicates
    combinations = list(set(combinations))
    return combinations


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
