"""Routines and class for Projected Bogoliubov MBPT diagrams."""
from __future__ import division

from builtins import str
from builtins import range
from past.utils import old_div
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
    vertices = [vert for vert in graph if not graph.nodes[vert]['operator']]

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
                      and not prop[3]['anomalous']) - 1
            new_graph[edge[0]][edge[1]][key]['anomalous'] = True
        anom_graphs.append(new_graph)

    # Loop to generate self-contractions
    for idx in range(len(anom_graphs)-1, -1, -1):
        iter_graph = anom_graphs[idx]
        test_vertices = []
        for vert in vertices:
            for _ in range(nbody_max - (old_div(iter_graph.degree(vert), 2))):
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

    __slots__ = ()

    def __init__(self, graph, unique_id, tag, child_tag):
        """Generate a PBMBPT diagram by copying a BMBPT one.

        Args:
            graph (NetworkX MultiDiGraph)): The graph of interest.

        """
        adg.bmbpt.BmbptFeynmanDiagram.__init__(self, graph, unique_id)
        self.tags = [tag, child_tag]
        self.set_io_degrees()

    def set_io_degrees(self):
        """Attribute the correct in- and out-degrees to a PBMBPT diagram."""
        unsort_io_degrees = []
        for node in self.graph:
            # Edges going out that are anomalous are annihilators going in
            nb_anom_out_edges = sum(1 for edge
                                    in self.graph.out_edges(nbunch=node,
                                                            data=True,
                                                            keys=True)
                                    if edge[3]['anomalous'])
            # And thus should be decounted to NetworkX's out_degree...
            out_degree = self.graph.out_degree(node) - nb_anom_out_edges
            # ...and added to NetworkX's in_degree
            in_degree = self.graph.in_degree(node) + nb_anom_out_edges

            unsort_io_degrees.append((in_degree, out_degree))
        self.unsort_io_degrees = tuple(unsort_io_degrees)
        self.io_degrees = tuple(sorted(self.unsort_io_degrees))

    def attribute_qp_labels(self):
        """Attribute the appropriate qp labels to the graph's propagators."""
        idx_counter = 1
        for prop in self.graph.edges(keys=True, data=True):
            if prop[3]['anomalous']:
                prop[3]['qp_state'] = "k_{%i}k_{%i}" % (idx_counter + 1,
                                                        idx_counter)
                idx_counter += 2
            else:
                prop[3]['qp_state'] = "k_{%i}" % idx_counter
                idx_counter += 1

    def vertex_expression(self, vertex):
        """Return the expression associated to a given vertex.

        Args:
            vertex (int): The vertex of interest in the graph.

        Returns:
            (str): The LaTeX expression associated to the vertex.

        """
        expression = r"\epsilon^{" \
            + "".join("%s"
                      % prop[3]['qp_state']
                      for prop
                      in self.graph.out_edges(vertex, keys=True, data=True)
                      if not prop[3]['anomalous']) \
            + "}_{" \
            + "".join("%s"
                      % (prop[3]['qp_state'].split("}")[1] + "}"
                         if prop[3]['anomalous'] else prop[3]['qp_state'])
                      for prop
                      in self.graph.in_edges(vertex, keys=True, data=True)) \
            + "".join("%s"
                      % (prop[3]['qp_state'].split("}")[0] + "}")
                      for prop
                      in self.graph.out_edges(vertex, keys=True, data=True)
                      if prop[3]['anomalous']) \
            + "}"
        return expression

    def extract_numerator(self):
        """Return the numerator associated to a PBMBPT graph.

        Returns:
            (str): The numerator of the graph.

        """
        graph = self.graph
        numerator = ""
        for vertex in graph:
            # Attribute the correct operator to each vertex
            numerator += r"\tilde{O}" if graph.nodes[vertex]['operator'] \
                else "\\Omega"
            # Attribute the good "type number" to each vertex
            numerator += "^{%i%i}_{" % (self.unsort_io_degrees[vertex][1],
                                        self.unsort_io_degrees[vertex][0])
            # First add the qp states corresponding to propagators going out
            numerator += "".join(prop[3]['qp_state']
                                 for prop in graph.out_edges(vertex,
                                                             keys=True,
                                                             data=True)
                                 if not prop[3]['anomalous'])

            # Add the qp states corresponding to propagators coming in
            previous_vertex = vertex - 1
            # Iterate over the lower vertices
            while previous_vertex >= 0:
                # Start with normal propagators
                numerator += "".join(
                    prop[3]['qp_state']
                    for prop in graph.in_edges(vertex, keys=True, data=True)
                    if prop[0] == previous_vertex and prop[0] != prop[1]
                    and not prop[3]['anomalous'])
                # Read anomalous propagators backwards
                anom_props = [prop for prop
                              in graph.in_edges(vertex, keys=True, data=True)
                              if prop[0] == previous_vertex
                              and prop[0] != prop[1]
                              and prop[3]['anomalous']]
                numerator += "".join(
                    (anom_props[idx][3]['qp_state'].split("}")[0] + "}")
                    for idx in range(len(anom_props)-1, -1, -1)
                    )
                previous_vertex -= 1
            # Add self-contractions
            numerator += "".join(
                (prop[3].split("}")[1] + "}" + prop[3].split("}")[0] + "}")
                for prop in graph.edges(vertex, keys=True, data='qp_state')
                if prop[0] == prop[1])
            # Add anomalous propagators entering from above
            numerator += "".join(prop[3]['qp_state'].split("}")[1] + "}"
                                 for prop in graph.out_edges(vertex,
                                                             keys=True,
                                                             data=True)
                                 if prop[3]['anomalous'] and prop[1] > prop[0])
            numerator += r"} (\varphi) " if graph.nodes[vertex]['operator'] \
                else "} "
        # Add the terms correspoding to anomalous propagators
        numerator += " ".join("R^{--}_{%s}(\\varphi)" % prop[3]['qp_state']
                              for prop
                              in graph.out_edges(keys=True, data=True)
                              if prop[3]['anomalous'])
        return numerator

    def multiplicity_symmetry_factor(self):
        """Return the symmetry factor associated with propagators multiplicity.

        Returns:
            (str): The symmetry factor associated with equivalent lines.

        """
        factor = ""
        # Account for up to three-body operators
        prop_multiplicity = [0 for _ in range(6)]
        for vertex_i in self.graph:
            for vertex_j in self.graph:
                nb_normal_props = 0
                nb_anomalous_props = 0
                for prop in self.graph.edges(vertex_i,
                                             keys=True,
                                             data='anomalous'):
                    if prop[1] == vertex_j and prop[3]:
                        nb_anomalous_props += 1
                    if prop[1] == vertex_j and not prop[3]:
                        nb_normal_props += 1
                if nb_anomalous_props >= 2:
                    prop_multiplicity[nb_anomalous_props - 1] += 1
                if nb_normal_props >= 2:
                    prop_multiplicity[nb_normal_props - 1] += 1

        for prop_id, multiplicity in enumerate(prop_multiplicity):
            if multiplicity == 1:
                factor += "(%i!)" % (prop_id+1)
            elif multiplicity >= 2:
                factor += "(%i!)" % (prop_id+1) + "^%i" % multiplicity
        return factor

    def write_diag_exps(self, latex_file, norder):
        """Write the expressions associated to a diagram in the LaTeX file.

        Args:
            latex_file (file): The LaTeX outputfile of the program.
            norder (int): The order in BMBPT formalism.

        """
        latex_file.write(
            "\\begin{align}\n\\text{PO}%i.%i.%i\n" % (norder,
                                                      (self.tags[0] + 1),
                                                      (self.tags[1] + 1))
            + "&= %s" % self.feynman_exp
            + r" \nonumber \\" + "\n"
            + "&= %s\\end{align}\n" % self.diag_exp)

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
