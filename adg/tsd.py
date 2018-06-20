"""Module with functions relative to time-stucture diagrams, called by ADG."""

import os
import math
import networkx as nx
import adg.diag


def time_structure_graph(graph):
    """Return the time-structure graph associated to the graph.

    Args:
        graph (NetwrokX MultiDiGraph): The BMBPT graph of interest.

    Returns:
        (NetworkX MultiDiGraph): The time-structure diagram.

    """
    time_graph = graph.to_directed()
    if time_graph.node[0]['operator']:
        for vertex in xrange(1, len(time_graph)):
            time_graph.add_edge(0, vertex)
    return adg.diag.to_skeleton(time_graph)


def tree_time_structure_den(time_graph):
    """Return the denominator associated to a tree time-structure graph.

    Args:
        time_graph (NetworkX MultiDiGraph): The TSD of interest.

    Returns:
        (str): The denominator associated to the TSD.

    """
    denominator = ""
    i = 0
    for vertex in time_graph:
        if not time_graph.node[vertex]['operator']:
            time_graph.node[vertex]['label'] = 'a_%i' % (i + 1)
            i += 1
    for vertex in time_graph:
        if not time_graph.node[vertex]['operator']:
            if time_graph.out_degree(vertex) == 0:
                denominator += time_graph.node[vertex]['label']
            else:
                denominator += "(%s" % time_graph.node[vertex]['label'] \
                    + "".join("+ %s" % time_graph.node[descendant]['label']
                              for descendant
                              in nx.descendants(time_graph, vertex)) + ")"
    return denominator


def equivalent_labelled_tsds(equivalent_trees, labelled_tsds):
    """Return the list of labelled TSDs corresponding to equivalent TSDs.

    Args:
        equivalent_trees (list): The equivalent tree TSDs of a non-tree TSD.
        labelled_tsds (list): The labelled TSDs obtained from BMBPT diagrams.

    Returns:
        (str): The list of tag numbers of the equivalent TSDs.

    """
    op_nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
    eq_labelled_tsds = ""
    for eq_tree_graph in equivalent_trees:
        for comp_tdiag in labelled_tsds:
            if sorted(tuple((eq_tree_graph.in_degree(node),
                             eq_tree_graph.out_degree(node))
                            for node in eq_tree_graph)) \
                    == comp_tdiag.io_degrees \
                    and comp_tdiag.is_tree:
                if nx.is_isomorphic(eq_tree_graph, comp_tdiag.graph, op_nm):
                    eq_labelled_tsds += " T%s," % (comp_tdiag.tags[0]+1)
                    break
    return "".join("%s." % eq_labelled_tsds.strip(','))


def write_section(latex_file, directory, pdiag, time_diagrams, nb_tree_tsds):
    """Write the appropriate section for tsd diagrams in the LaTeX file.

    Args:
        latex_file (file): The LaTeX output file of the program.
        directory (str): Path to the output folder.
        pdiag (bool): ``True`` if diagrams are to be drawn.
        time_diagrams (list): The ensemble of TSDs.
        nb_tree_tsds (int): Number of tree TSDs.

    """
    latex_file.write("\\section{Time-structure diagrams}\n\n"
                     + "\\subsection{Tree diagrams}\n\n")
    for tdiag in time_diagrams:
        if tdiag.tags[0] == nb_tree_tsds:
            latex_file.write("\\subsection{Non-tree diagrams}\n\n")
        latex_file.write("\\paragraph{Time-structure diagram T%i:}\n"
                         % (tdiag.tags[0]+1))
        if pdiag:
            time_file = open(directory
                             + "/Diagrams/time_%i.tex" % tdiag.tags[0])
            latex_file.write('\n\\begin{center}\n%s\n\\end{center}\n\n'
                             % time_file.read())
        latex_file.write("\\begin{equation}\n\\text{T%i} = "
                         "%s\\end{equation}\n\n"
                         % (tdiag.tags[0]+1, tdiag.expr))
        if tdiag.is_tree:
            latex_file.write("Resummation power: %i\n\n" % tdiag.resum)
        else:
            latex_file.write("Equivalent tree diagrams: %s\n\n"
                             % equivalent_labelled_tsds(tdiag.equivalent_trees,
                                                        time_diagrams))
            latex_file.write('\n\\begin{center}\n')
            tdiag.draw_equivalent_tree_tsds(latex_file)
            latex_file.write('\n\\end{center}\n\n')
        latex_file.write("Number of related Feynman diagrams: %i.\n\n"
                         % (len(tdiag.tags)-1))
        feynman_diags = ",".join(" %i" % (tag+1) for tag in tdiag.tags[1:])
        latex_file.write("Related Feynman diagrams:%s.\n\n" % feynman_diags)


def disentangle_cycle(time_graph, cycle_nodes):
    """Separate a cycle in a sum of tree diagrams.

    Args:
        time_graph (NetworkXn MultiDiGraph): A time-structure diagram.
        cycle_nodes (tuple): Integers encoding the positions of the end nodes
            of the cycle.

    Returns:
        (list): New graphs produced from treating the cycles in the TSD.

    """
    paths = list(nx.all_simple_paths(time_graph,
                                     cycle_nodes[0],
                                     cycle_nodes[1]))
    new_graphs = []
    insert_node = paths[0][1]
    for daughter_node in paths[1][1:]:
        new_graph = time_graph.to_directed()
        new_graph.add_edge(insert_node, daughter_node)
        for test_node in paths[1]:
            if test_node in list(time_graph.predecessors(daughter_node)):
                mother_node = test_node
                break
        new_graph.add_edge(mother_node, insert_node)
        adg.diag.to_skeleton(new_graph)
        new_graphs.append(new_graph)
    return new_graphs


def find_cycle(graph):
    """Return start and end nodes for an elementary cycle.

    Args:
        graph (NetworkX MultiDiGraph): A TSD with cycle(s) to be treated.

    Returns:
        (tuple): Positions of the two end nodes of a cycle in the graph.

    """
    cycle_found = False
    for node_a in (node for node in graph if graph.out_degree(node) >= 2):
        for node_b in (node for node in graph if graph.in_degree(node) >= 2):
            paths = list(nx.all_simple_paths(graph, node_a, node_b))
            if len(paths) >= 2:
                cycle_nodes = (node_a, node_b)
                cycle_found = True
                # Avoid false positive when node_a has 2+ branches out
                # but only one goes to node_b
                for test_node in paths[0][1:-1]:
                    if test_node in paths[1][1:-1]:
                        cycle_found = False
                        break
                break
        if cycle_found:
            break
    return cycle_nodes


def treat_tsds(diagrams_time):
    """Order TSDs, produce their expressions, return also number of trees.

    Args:
        diagrams_time (list): All the associated TSDs.

    Returns:
        (tuple): List of TSDs, number of tree TSDs

    """
    tree_tsds = []
    for i_diag in xrange(len(diagrams_time)-1, -1, -1):
        if diagrams_time[i_diag].is_tree:
            tree_tsds.append(diagrams_time[i_diag])
            del diagrams_time[i_diag]

    adg.diag.topologically_distinct_diagrams(tree_tsds)
    adg.diag.topologically_distinct_diagrams(diagrams_time)

    diagrams_time = tree_tsds + diagrams_time

    for index, t_diag in enumerate(diagrams_time):
        t_diag.tags.insert(0, index)
        if not t_diag.is_tree:
            t_diag.equivalent_trees = t_diag.treat_cycles()
            t_diag.expr = " + ".join("\\frac{1}{%s}"
                                     % adg.tsd.tree_time_structure_den(graph)
                                     for graph
                                     in t_diag.equivalent_trees)
    return diagrams_time, len(tree_tsds)


class TimeStructureDiagram(adg.diag.Diagram):
    """Describes a time-structure diagram with its related properties.

    Attributes:
        perms (dict): The permutations on the vertices for all the BMBPT
            diagrams associated to this TSD.
        equivalent_trees (list): The tag numbers of the equivalent tree TSDs
            associated to a non-tree TSD.
        is_tree (bool): The tree or non-tree character of a TSD.
        expr (str): The Goldstone denominator associated to the TSD.

    """

    def __init__(self, bmbpt_diag, tag_num):
        """Generate a tsd diagram out of a BMBPT one.

        Args:
            bmbpt_diag (BmbptFeynmanDiagram): The BMBPT graph used to be
                turned into a TSD.
            tag_num (int): The number associated to the TSD.

        """
        adg.diag.Diagram.__init__(self, time_structure_graph(bmbpt_diag.graph))
        self.tags = [tag_num]
        self.perms = {tag_num: {i: i for i in xrange(len(self.graph))}}
        self.equivalent_trees = []
        if nx.is_arborescence(self.graph):
            self.is_tree = True
            self.expr = "\\frac{1}{%s}" % tree_time_structure_den(self.graph)
            self.resum = self.resummation_power()
        else:
            self.is_tree = False
            self.expr = ""
            self.resum = 0

    def treat_cycles(self):
        """Find and treat cycles in a TSD diagram.

        Returns:
            (list): The unique tree TSDs associated to a non-tree TSD.

        """
        graphs = [self.graph]
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

    def draw_equivalent_tree_tsds(self, latex_file):
        """Draw the equivalent tree TSDs for a given non-tree TSD.

        Args:
            latex_file (file): The output LaTeX file of the priogram.

        """
        for index, graph in enumerate(self.equivalent_trees):
            adg.diag.feynmf_generator(graph,
                                      'MBPT',
                                      'equivalent%i_%i' % (self.tags[0],
                                                           index))
            diag_file = open("equivalent%i_%i.tex" % (self.tags[0], index))
            latex_file.write(diag_file.read())
            diag_file.close()
            os.unlink("./equivalent%i_%i.tex" % (self.tags[0], index))

    def resummation_power(self):
        """Calculate the resummation power of the tree TSD.

        Returns:
            (int): The resummation power associated to the TSD.abs

        """
        power = math.factorial(len(self.graph) - 1)

        for node in range(1, len(self.graph)):
            power /= 1 + len(nx.descendants(self.graph, node))

        return power
