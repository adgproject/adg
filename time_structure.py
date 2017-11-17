"""Module with functions relative to time-stucture diagrams, called by ADG."""

import os
import string
import networkx as nx
import methods as mth


def time_structure_graph(graph):
    """Return the time-structure graph associated to the graph."""
    time_graph = graph.to_directed()
    if time_graph.node[0]['operator']:
        for vertex in xrange(1, len(time_graph)):
            time_graph.add_edge(0, vertex)
    return mth.to_skeleton(time_graph)


def has_tree_time_structure(graph):
    """Return True if the time structure of the graph is a tree."""
    diag_copy = graph.to_directed()
    for vertex in xrange(1, len(diag_copy)):
        if diag_copy.in_degree(vertex) == 0:
            diag_copy.add_edge(0, vertex)
    time_diag = nx.DiGraph()
    for vertex_i in diag_copy:
        for vertex_j in diag_copy:
            lgst_path = []
            for path in nx.all_simple_paths(diag_copy,
                                            vertex_i,
                                            vertex_j):
                if len(path) > len(lgst_path):
                    lgst_path = path
            time_diag.add_path(lgst_path)
    return nx.is_arborescence(time_diag)


def tree_time_structure_den(time_graph):
    """Return the denominator associated to a tree time-structure graph."""
    denominator = ""
    labels = list(string.ascii_lowercase)
    i = 0
    for vertex in time_graph:
        if not time_graph.node[vertex]['operator']:
            time_graph.node[vertex]['label'] = labels[i]
            i += 1
    for vertex in time_graph:
        if not time_graph.node[vertex]['operator']:
            if time_graph.out_degree(vertex) == 0:
                denominator += time_graph.node[vertex]['label']
            else:
                denominator += "(%s" % time_graph.node[vertex]['label'] \
                    + "".join("+ %s" % time_graph.node[descendant]['label']
                              for descendant in nx.descendants(time_graph,
                                                               vertex)) + ")"
    return denominator


def write_time_diagrams_section(latex_file, directory, pdiag, pdraw,
                                time_diagrams):
    """Write the appropriate section for tst diagrams in the LaTeX file."""
    latex_file.write("\\section{Associated time-structure diagrams}\n\n")
    for tdiag in time_diagrams:
        latex_file.write("\\paragraph{Time-structure diagram T%i:}\n"
                         % (tdiag.tags[0]+1))
        if pdiag and pdraw:
            time_file = open(directory
                             + "/Diagrams/time_%i.tex" % tdiag.tags[0])
            latex_file.write('\n\\begin{center}\n%s\n\\end{center}\n\n'
                             % time_file.read())
        latex_file.write("Tree: Yes\n\n" if tdiag.is_tree else "Tree: No\n\n")
        latex_file.write("\\begin{equation}\n%s\\end{equation}\n" % tdiag.expr)
        if not tdiag.is_tree:
            latex_file.write("Equivalent tree diagrams:\n\n")
            latex_file.write('\n\\begin{center}\n')
            for index, graph in enumerate(tdiag.equivalent_trees):
                mth.feynmf_generator(graph,
                                     'MBPT',
                                     'equivalent%i_%i' % (tdiag.tags[0],
                                                          index))
                diag_file = open("equivalent%i_%i.tex" % (tdiag.tags[0],
                                                          index))
                latex_file.write(diag_file.read())
                diag_file.close()
                os.unlink("./equivalent%i_%i.tex" % (tdiag.tags[0], index))
            latex_file.write('\n\\end{center}\n\n')
        feynman_diags = ",".join(" %i" % (tag+1) for tag in tdiag.tags[1:]) \
            + ".\n\n"
        latex_file.write("Related Feynman diagrams:%s" % feynman_diags)


def treat_cycles(time_graph):
    """Find and treat cycles in a TSD diagram."""
    graphs = [time_graph]
    tree_graphs = []
    cycles_left = True
    while cycles_left:
        for gr_index in xrange(len(graphs)-1, -1, -1):
            graphs += disentangle_cycle(graphs[gr_index],
                                        find_cycle(graphs[gr_index]))
            del graphs[gr_index]
        cycles_left = False
        for graph_indx in xrange(len(graphs)-1, -1, -1):
            if nx.is_arborescence(graphs[graph_indx]):
                tree_graphs.append(graphs[graph_indx])
                del graphs[graph_indx]
            else:
                cycles_left = True
    tree_graphs_uniq = []
    for t_graph in tree_graphs:
        new_graph = True
        for t_graph_uniq in tree_graphs_uniq:
            if nx.edges(t_graph) == nx.edges(t_graph_uniq):
                new_graph = False
                break
        if new_graph:
            tree_graphs_uniq.append(t_graph)
    return tree_graphs_uniq


def disentangle_cycle(time_graph, cycle_nodes):
    """Separate a cycle in a sum of tree diagrams."""
    paths = list(nx.all_simple_paths(time_graph,
                                     cycle_nodes[0],
                                     cycle_nodes[1]))
    new_graphs = []
    for insert_node in paths[0][1:-1]:
        for daughter_node in paths[1][1:]:
            new_graph = time_graph.to_directed()
            new_graph.add_edge(insert_node, daughter_node)
            for test_node in paths[1]:
                if test_node in list(time_graph.predecessors(daughter_node)):
                    mother_node = test_node
            new_graph.add_edge(mother_node, insert_node)
            mth.to_skeleton(new_graph)
            new_graphs.append(new_graph)
    return new_graphs


def find_cycle(graph):
    """Return start and end nodes for an elementary cycle."""
    cycle_found = False
    for node_a in graph:
        if graph.out_degree(node_a) >= 2:
            for node_b in graph:
                if (graph.in_degree(node_b) == 2
                        and len(list(nx.all_simple_paths(graph,
                                                         node_a,
                                                         node_b))) == 2) \
                    or (graph.in_degree(node_b) > 2
                        and len(list(nx.all_simple_paths(graph,
                                                         node_a,
                                                         node_b))) > 2
                        and graph.out_degree(node_a) == graph.in_degree(node_b)):
                    cycle_nodes = (node_a, node_b)
                    cycle_found = True
                    break
        if cycle_found:
            break
    return cycle_nodes


class TimeStructureDiagram(mth.Diagram):
    """Describes a time-structure diagram with its related properties."""

    def __init__(self, bmbpt_diag, tag_num):
        """Generate a TST diagram out of a BMBPT one."""
        mth.Diagram.__init__(self, time_structure_graph(bmbpt_diag.graph))
        self.tags = [tag_num]
        self.equivalent_trees = []
        if nx.is_arborescence(self.graph):
            self.is_tree = True
            self.expr = "\\frac{1}{%s}" % tree_time_structure_den(self.graph)
        else:
            self.is_tree = False
            self.expr = ""
