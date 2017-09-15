"""Module containg methods relative to time-stucture diagrams,
    to be called by ADG."""

import string
import networkx as nx
import methods as mth


def time_structure_graph(diagram):
    """Return the time-structure graph associated to the diagram."""
    time_diag = diagram.to_directed()
    if time_diag.node[0]['operator']:
        for vertex in xrange(1, len(time_diag)):
            time_diag.add_edge(0, vertex)
    for vertex_a in time_diag:
        for vertex_b in time_diag:
            while time_diag.number_of_edges(vertex_a, vertex_b) > 1:
                time_diag.remove_edge(vertex_a, vertex_b)
            if len(list(nx.all_simple_paths(time_diag,
                                            vertex_a, vertex_b))) > 1:
                while len(nx.shortest_path(time_diag,
                                           vertex_a, vertex_b)) == 2:
                    time_diag.remove_edge(vertex_a, vertex_b)
    return time_diag


def has_tree_time_structure(diagram):
    """Return True if the time structure of the diagram is a tree."""
    diag_copy = diagram.to_directed()
    for vertex in xrange(1, len(diag_copy)):
        if diag_copy.in_degree(vertex) == 0:
            diag_copy.add_edge(0, vertex)
    time_diag = nx.DiGraph()
    for vertex_i in diag_copy:
        for vertex_j in diag_copy:
            lgst_path = []
            for path in nx.all_simple_paths(diag_copy,
                                            source=vertex_i,
                                            target=vertex_j):
                if len(path) > len(lgst_path):
                    lgst_path = path
            time_diag.add_path(lgst_path)
    return nx.is_arborescence(time_diag)


def tree_time_structure_den(time_diagram):
    """Return the denominator associated to a tree time-structure diagram."""
    denominator = ""
    labels = list(string.ascii_lowercase)
    i = 0
    for vertex in time_diagram:
        if not time_diagram.node[vertex]['operator']:
            time_diagram.node[vertex]['label'] = labels[i]
            i += 1
    for vertex in time_diagram:
        if not time_diagram.node[vertex]['operator']:
            if time_diagram.out_degree(vertex) == 0:
                denominator += time_diagram.node[vertex]['label']
            else:
                denominator += "(" + time_diagram.node[vertex]['label']
                for descendant in nx.descendants(time_diagram, vertex):
                    denominator += "+" + time_diagram.node[descendant]['label']
                denominator += ")"
    return denominator


class TimeStructureDiagram(mth.Diagram):
    """Describes a time-structure diagram with its related properties."""

    def __init__(self, bmbpt_diag, tag_num):
        mth.Diagram.__init__(self, time_structure_graph(bmbpt_diag.graph))
        self.tags = [tag_num]
        if nx.is_arborescence(self.graph):
            self.expr = "".join("\\frac{1}{"
                                + tree_time_structure_den(self.graph) + "}")
        else:
            self.expr = ""
