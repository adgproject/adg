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

    for edge in graph.edges(keys=True, data=True):
        edge[3]['anomalous'] = False

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
                      in new_graph.edges(keys=True, data=True)
                      if (edge[0], edge[1]) == (prop[0], prop[1])
                      and prop[3]['anomalous'])
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


def order_diagrams(diagrams):
    """Order the PBMBPT diagrams and return number of diags for each type.

    Args:
        diagrams (list): Possibly redundant BmbptFeynmanDiagrams.

    Returns:
        (tuple): First element is the list of topologically unique, ordered
            diagrams. Second element is a dict with the number of diagrams
            for each major type. Third element is a dict with the identifiers
            of diagrams starting each output file section.

    """
    diagrams_2_hf = []
    diagrams_2_ehf = []
    diagrams_2_not_hf = []
    diagrams_3_hf = []
    diagrams_3_ehf = []
    diagrams_3_not_hf = []

    for i_diag in xrange(len(diagrams)-1, -1, -1):
        if diagrams[i_diag].two_or_three_body == 2:
            if diagrams[i_diag].hf_type == "HF":
                diagrams_2_hf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "EHF":
                diagrams_2_ehf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "noHF":
                diagrams_2_not_hf.append(diagrams[i_diag])
        elif diagrams[i_diag].two_or_three_body == 3:
            if diagrams[i_diag].hf_type == "HF":
                diagrams_3_hf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "EHF":
                diagrams_3_ehf.append(diagrams[i_diag])
            elif diagrams[i_diag].hf_type == "noHF":
                diagrams_3_not_hf.append(diagrams[i_diag])
        del diagrams[i_diag]

    diagrams = diagrams_2_hf + diagrams_2_ehf + diagrams_2_not_hf \
        + diagrams_3_hf + diagrams_3_ehf + diagrams_3_not_hf

    diags_nb_per_type = {
        'nb_2_hf': len(diagrams_2_hf),
        'nb_2_ehf': len(diagrams_2_ehf),
        'nb_2_not_hf': len(diagrams_2_not_hf),
        'nb_3_hf': len(diagrams_3_hf),
        'nb_3_ehf': len(diagrams_3_ehf),
        'nb_3_not_hf': len(diagrams_3_not_hf),
        'nb_diags': len(diagrams),
        'nb_2': (len(diagrams_2_hf) + len(diagrams_2_ehf)
                 + len(diagrams_2_not_hf)),
        'nb_3': (len(diagrams_3_hf) + len(diagrams_3_ehf)
                 + len(diagrams_3_not_hf))
    }

    section_flags = {
        'two_body_hf': diagrams_2_hf[0].unique_id if diagrams_2_hf else -1,
        'two_body_ehf': diagrams_2_ehf[0].unique_id if diagrams_2_ehf else -1,
        'two_body_not_hf': diagrams_2_not_hf[0].unique_id
                           if diagrams_2_not_hf else -1,
        'three_body_hf': diagrams_3_hf[0].unique_id if diagrams_3_hf else -1,
        'three_body_ehf': diagrams_3_ehf[0].unique_id
                          if diagrams_3_ehf else -1,
        'three_body_not_hf': diagrams_3_not_hf[0].unique_id
                             if diagrams_3_not_hf else -1
    }

    return diagrams, diags_nb_per_type, section_flags


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

    __slots__ = ()

    def __init__(self, graph, unique_id, tag, child_tag):
        """Generate a PBMBPT diagram by copying a BMBPT one.

        Args:
            graph (NetworkX MultiDiGraph)): The graph of interest.

        """
        adg.bmbpt.BmbptFeynmanDiagram.__init__(self, graph, unique_id)
        self.tags = [tag, child_tag]

    def attribute_qp_labels(self):
        """Attribute the appropriate qp labels to the graph's propagators."""
        idx_counter = 1
        for prop in self.graph.edges(keys=True, data=True):
            if prop[3]['anomalous']:
                prop[3]['qp_state'] = "k_{%i}k_{%i}" % (idx_counter,
                                                        idx_counter + 1)
                idx_counter += 2
            else:
                prop[3]['qp_state'] = "k_{%i}" % idx_counter
                idx_counter += 1

    def vertex_expression(self, vertex):
        """Return the expression associated to a given vertex.

        Args:
            vertex (int): The vertex of interest in the graph.

        """
        expression = r"\epsilon^{" \
            + "".join("%s"
                      % prop[3]['qp_state']
                      for prop
                      in self.graph.out_edges(vertex, keys=True, data=True)
                      if not prop[3]['anomalous']) \
            + "}_{" \
            + "".join("%s"
                      % (prop[3]['qp_state'].split("}")[0] + "}")
                      for prop
                      in self.graph.in_edges(vertex, keys=True, data=True)) \
            + "".join("%s"
                      % (prop[3]['qp_state'].split("}")[1] + "}")
                      for prop
                      in self.graph.out_edges(vertex, keys=True, data=True)
                      if prop[3]['anomalous']) \
            + "}"
        return expression

    def write_graph(self, latex_file, directory, write_time):
        """Write the PBMBPT graph and its associated TSD to the LaTeX file.

        Args:
            latex_file (file): The LaTeX output file of the program.
            directory (str): The path to the result folder.
            write_time (bool): ``True`` if we want informations on the
                associated TSDs.

        """
        latex_file.write('\n\\begin{center}\n')
        adg.diag.draw_diagram(directory, latex_file,
                              "%i_%i" % (self.tags[0], self.tags[1]), 'diag')
        if write_time:
            latex_file.write(
                '\\hspace{10pt} $\\rightarrow$ \\hspace{10pt} T%i:'
                % (self.time_tag + 1))
            adg.diag.draw_diagram(directory, latex_file,
                                  str(self.time_tag), 'time')
        latex_file.write('\n\\end{center}\n\n')
