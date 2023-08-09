"""Routines and class for Projected Bogoliubov MBPT diagrams."""
from __future__ import division

from builtins import str
from builtins import range
from past.utils import old_div
from operator import xor
import copy
import itertools
import networkx as nx
import adg.bmbpt
import adg.tsd
import adg.diag


def generate_anomalous_diags(diag, nbody_max):
    """Generate PBMBPT graphs with anomalous lines, with some redundancy.

    Args:
        diag (BmbptFeynmanDiagram): The diagram to generate children from.
        nbody_max (int): The maximal n-body character of a graph vertex.

    Returns:
        list: The anomalous graphs generated.

    """
    graph = diag.graph
    anom_graphs = [graph]
    vertices = [vert for vert in graph if not graph.nodes[vert]['operator']]

    nx.set_edge_attributes(graph, False, 'anomalous')

    # Turn normal propagators in anomalous ones
    tweakable_edges = []
    for vert1 in vertices:
        for vert2 in vertices[vert1:]:
            tweakable_edges += [(vert1, vert2) for _
                                in range(graph.number_of_edges(vert1, vert2))]

    for comb in unique_edge_combinations(tweakable_edges,
                                         diag.equivalent_permutations()):
        new_graph = copy.deepcopy(graph)
        for edge in comb:
            key = sum(1 for prop
                      in new_graph.out_edges(edge[0], keys=True, data=True)
                      if edge[1] == prop[1] and not prop[3]['anomalous']) - 1
            new_graph[edge[0]][edge[1]][key]['anomalous'] = True
        anom_graphs.append(new_graph)

    # Loop to generate self-contractions
    for iter_graph in reversed(anom_graphs):
        test_vertices = []
        for v in vertices:
            test_vertices += [v for _
                              in range(nbody_max
                                       - (old_div(iter_graph.degree(v), 2)))]
        if test_vertices:
            for comb in unique_vertex_combinations(
                    test_vertices, equiv_generating_permutations(iter_graph)):
                new_graph = copy.deepcopy(iter_graph)
                new_graph.add_edges_from(((vert, vert) for vert in comb),
                                         anomalous=True)
                anom_graphs.append(new_graph)
    return anom_graphs


def generate_combinations(iter_list):
    """Generate all possible combinations of length 1 to total.

    Args:
        iter_list (list): A list of iterable objects.

    Returns:
        list: A list with all the possible combinations of all lengths.

    >>> print(generate_combinations([1, 2, 3]))
    [(1,), (1, 2), (1, 2, 3), (1, 3), (2,), (2, 3), (3,)]

    """
    combinations = []
    for i in range(1, len(iter_list) + 1):
        combinations += [k for k in itertools.combinations(iter_list, i)]
    # Remove duplicates
    return sorted(set(combinations))


def unique_edge_combinations(edges, permutations):
    """Return all edge combinations not producing equivalent anomalous graphs.

    Args:
        edges (list): The edges that can be modified.
        permutations (list): The permutation generating equivalent diagrams.

    Returns:
        list: The list of edges producing unique anomalous diagrams.

    >>> edges = [(1, 3), (2, 3)]
    >>> permutations = [{1: 1, 2: 2}, {1: 2, 2: 1}]
    >>> print(unique_edge_combinations(edges, permutations))
    [((1, 3), (2, 3)), ((2, 3),)]

    """
    edge_combs = generate_combinations(edges)

    if len(permutations) <= 1:
        return edge_combs

    unique_edge_combs = copy.deepcopy(edge_combs)
    for idx, comb1 in enumerate(edge_combs):
        for comb2 in edge_combs[idx+1:]:
            is_same_perm = False
            if len(comb1) == len(comb2):
                for permutation in permutations[1:]:
                    new_comb = []
                    for edge in comb2:
                        new_0 = permutation[edge[0]] \
                            if edge[0] in permutation.keys() else edge[0]
                        new_1 = permutation[edge[1]] \
                            if edge[1] in permutation.keys() else edge[1]
                        new_comb.append(tuple((new_0, new_1)))
                    if sorted(comb1) == sorted(new_comb):
                        is_same_perm = True
                        break
            if is_same_perm:
                unique_edge_combs.remove(comb1)
                break
    return unique_edge_combs


def unique_vertex_combinations(vertices, permutations):
    """Return vertex combinations generating unique anomalous diagrams.

    Return combinations of vertices on which self-contractions can be added
    without producing topologically equivalent PBMBPT diagrams.

    Args:
        vertices (list): Vertices that can be self-contracted.
        permutations (list): The permutations that generate equivalent diags.

    Returns:
        list: Vertex combinations that do not produce equivalent diags.

    >>> vertices = [1, 3]
    >>> permutations = [{1: 1, 3: 3}, {1: 3, 3: 1}]
    >>> print(unique_vertex_combinations(vertices, permutations))
    [(1, 3), (3,)]

    """
    combinations = generate_combinations(vertices)
    if len(combinations) <= 1:
        return combinations
    unique_combs = copy.deepcopy(combinations)
    for idx, comb1 in enumerate(combinations):
        for comb2 in combinations[idx+1:]:
            is_same_perm = False
            if len(comb1) == len(comb2):
                # Skip the identity
                for permutation in permutations[1:]:
                    new_comb2 = [permutation[item]
                                 if item in permutation.keys() else item
                                 for item in comb2]
                    if sorted(comb1) == sorted(new_comb2):
                        is_same_perm = True
                        break
            if is_same_perm:
                unique_combs.remove(comb1)
                break
    return unique_combs


def equiv_generating_permutations(graph):
    """Return the list of permutations generating equivalent PBMBPT diags.

    Args:
        graph (Networkx MultiDiGraph): The graph to be checked.

    Returns:
        list: The mappings giving equivalent graphs, inc. identity.

    """
    op_nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
    anom_em = nx.algorithms.isomorphism.categorical_multiedge_match('anomalous',
                                                                    False)
    unsort_io_degrees = []
    for node in graph:
        # Edges going out that are anomalous are annihilators going in
        nb_anom_out_edges = sum(1 for edge
                                in graph.out_edges(node, data=True, keys=True)
                                if edge[3]['anomalous'])
        # And thus should be decounted to NetworkX's out_degree...
        out_degree = graph.out_degree(node) - nb_anom_out_edges
        # ...and added to NetworkX's in_degree
        in_degree = graph.in_degree(node) + nb_anom_out_edges

        unsort_io_degrees.append((in_degree, out_degree))
    unsort_io_degrees = tuple(unsort_io_degrees)

    perm_vertices = [vertex for vertex, degrees in enumerate(unsort_io_degrees)
                     if not graph.nodes[vertex]['operator']
                     and unsort_io_degrees.count(degrees) >= 2]
    permutations = []
    doubled_graph = adg.diag.create_checkable_diagram(graph)
    for permutation in itertools.permutations(perm_vertices):
        permuted_graph = nx.relabel_nodes(doubled_graph,
                                          dict(list(zip(perm_vertices,
                                                        permutation))),
                                          copy=True)
        intersection = copy.deepcopy(doubled_graph)
        intersection.remove_edges_from(e for e in doubled_graph.edges()
                                       if e not in permuted_graph.edges())
        matcher = nx.algorithms.isomorphism.DiGraphMatcher(doubled_graph,
                                                           intersection,
                                                           node_match=op_nm,
                                                           edge_match=anom_em)
        if matcher.is_isomorphic():
            permutations.append(dict(list(zip(perm_vertices, permutation))))
    return permutations


def filter_new_diagrams(new_diags, old_diags):
    """Eliminate diagrams having a topologically equivalent diag.

    Args:
        new_diags (list): The list of newly created PBMBPT diagrams.
        old_diags (list): The list of already checked PBMBPT diagrams.

    """
    iso = nx.algorithms.isomorphism
    op_nm = iso.categorical_node_match('operator', False)
    anom_em = iso.categorical_multiedge_match('anomalous', False)
    for ind, new_diag in adg.tools.reversed_enumerate(new_diags):
        # Only anomalous props that are not self-contractions can produce
        # topologically equivalent diagrams through our generation process
        if new_diag.has_anom_non_selfcontracted_props():
            for new_diag_2 in new_diags[:ind]:
                if new_diag_2.io_degrees == new_diag.io_degrees \
                        and new_diag_2.has_anom_non_selfcontracted_props():
                    matcher = iso.DiGraphMatcher(new_diag.check_graph,
                                                 new_diag_2.check_graph,
                                                 node_match=op_nm,
                                                 edge_match=anom_em)
                    if matcher.is_isomorphic():
                        del new_diags[ind]
                        break
            else:
                # Go through the list backwards: first items are BMBPT diags
                for old_diag in reversed(old_diags):
                    if not isinstance(old_diag, ProjectedBmbptDiagram):
                        # All PBMBPT diags have been parsed
                        break
                    if old_diag.io_degrees == new_diag.io_degrees \
                            and old_diag.has_anom_non_selfcontracted_props():
                        matcher = iso.DiGraphMatcher(old_diag.check_graph,
                                                     new_diag.check_graph,
                                                     node_match=op_nm,
                                                     edge_match=anom_em)
                        if matcher.is_isomorphic():
                            del new_diags[ind]
                            break


class ProjectedBmbptDiagram(adg.bmbpt.BmbptFeynmanDiagram):
    """Describes a PBMBPT diagram with its related properties.

    Args:
        graph (NetworkX MultiDiGraph): The graph of interest.
        unique_id (int): The unique number associated to the diagram.
        tag (int): Tag of the parent BMBPT diagram.
        child_tag (int): Identifier of the PBMBPT diagram within the children
            of the BMBPT diagram.

    """

    __slots__ = ('_check_graph',)

    def __init__(self, graph, unique_id, tag, child_tag):
        """Generate a PBMBPT diagram by copying a BMBPT one."""
        adg.bmbpt.BmbptFeynmanDiagram.__init__(self, graph, unique_id)
        self.tags = [tag, child_tag]
        self.set_io_degrees()
        self._check_graph = None
        """NetworkX MultiDiGraph: A copy of the graph that can be used for
        topological equivalence checks (lazy-initialized)."""

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
        self._check_graph = None

    def extract_integral(self):
        """Return the integral part of the Feynman expression of the diag.

        Returns:
            str: The integral part of its Feynman expression.

        """
        pert_vertex_indices = list(range(1, len(self.graph)))
        integral = "".join("\\mathrm{d}\\tau_%i" % vertex
                           for vertex in pert_vertex_indices)

        normal_props = [prop for prop in self.graph.edges(data='anomalous')
                        if not prop[2]]

        if len(pert_vertex_indices) > 1:
            for vertex_i in pert_vertex_indices:
                integral += "".join("\\theta(\\tau_%i-\\tau_%i) " % (vertex_j,
                                                                     vertex_i)
                                    for vertex_j in pert_vertex_indices
                                    if (vertex_i, vertex_j, False)
                                    in normal_props)
        integral += "".join("e^{-\\tau_%i %s}"
                            % (vertex, self.vert_exp[vertex])
                            for vertex in pert_vertex_indices)
        return integral

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
            str: The LaTeX expression associated to the vertex.

        """
        expression = r"\epsilon^{" \
            + "".join("%s"
                      % prop[3]['qp_state']
                      for prop
                      in self.graph.out_edges(vertex, keys=True, data=True)
                      if not prop[3]['anomalous']) \
            + "}_{" \
            + "".join("%s"
                      % (prop[3]['qp_state'].split("}")[0] + "}"
                         if prop[3]['anomalous'] else prop[3]['qp_state'])
                      for prop
                      in self.graph.in_edges(vertex, keys=True, data=True)) \
            + "".join("%s"
                      % (prop[3]['qp_state'].split("}")[1] + "}")
                      for prop
                      in self.graph.out_edges(vertex, keys=True, data=True)
                      if prop[3]['anomalous']) \
            + "}"
        return expression

    def has_sign_factor(self):
        """Return True if a sign factor is associated to the diagram.

        Wrapper allowing for easy refactoring of expression code.

        Returns:
            bool: The presence of a sign factor.

        """
        # Use exclusive or for the sign factor
        return xor(self.has_crossing_sign(), self.has_anom_props_linked_sign())

    def extract_numerator(self):
        """Return the numerator associated to a PBMBPT graph.

        Returns:
            str: The numerator of the graph.

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
                    if prop[0] == previous_vertex and not prop[3]['anomalous'])
                # Read anomalous propagators backwards
                anom_props = [prop for prop
                              in graph.in_edges(vertex, keys=True, data=True)
                              if prop[0] == previous_vertex
                              and prop[3]['anomalous']]
                numerator += "".join(
                    (anom_prop[3]['qp_state'].split("}")[0] + "}")
                    for anom_prop in reversed(anom_props)
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

    def symmetry_factor(self):
        """Return the overall symmetry factor of the diagram.

        Returns:
            str: The combination of all symmetry factors.

        """
        sym_factor = ""
        vertex_sym = 1
        anom_factor = self.anomalous_contractions_factor()
        for vertex_degrees in self.unsort_io_degrees:
            if self.unsort_io_degrees.count(vertex_degrees) >= 2:
                vertex_sym = self.vertex_exchange_sym_factor
                break
        sym_factor += "%i" % (vertex_sym*anom_factor) \
            if (vertex_sym*anom_factor) > 1 else ""
        sym_factor += self.multiplicity_symmetry_factor()
        return sym_factor

    def multiplicity_symmetry_factor(self):
        """Return the symmetry factor associated with propagators multiplicity.

        Returns:
            str: The symmetry factor associated with equivalent lines.

        """
        factor = ""
        # Account for up to three-body operators
        prop_multiplicity = [0 for _ in range(6)]
        for vertex_i in self.graph:
            for vertex_j in self.graph:
                nb_normal_props = sum(1 for prop
                                      in self.graph.edges(vertex_i,
                                                          keys=True,
                                                          data='anomalous')
                                      if prop[1] == vertex_j and not prop[3])
                nb_anomalous_props = sum(1 for prop
                                         in self.graph.edges(vertex_i,
                                                             keys=True,
                                                             data='anomalous')
                                         if prop[1] == vertex_j and prop[3])
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

    def anomalous_contractions_factor(self):
        """Return the factor associated with anomalous self-contractions.

        Returns:
            int: The anomalous self-contractions factor.

        """
        nb_self_contractions = len(list(nx.selfloop_edges(self.graph)))
        return pow(2, nb_self_contractions)

    def equivalent_permutations(self):
        """Return the permutations generating equivalent diagrams.

        Returns:
            list: Vertices permutations as dictionnaries.

        """
        op_nm = nx.algorithms.isomorphism.categorical_node_match('operator',
                                                                 False)
        anom_em = nx.algorithms.isomorphism.categorical_multiedge_match(
            'anomalous', False)

        perm_vertices = [vertex for vertex, degrees
                         in enumerate(self.unsort_io_degrees)
                         if not self.graph.nodes[vertex]['operator']
                         and self.unsort_io_degrees.count(degrees) >= 2]
        permutations = []
        for perm in itertools.permutations(perm_vertices):
            permuted_graph = nx.relabel_nodes(self.check_graph,
                                              dict(list(zip(perm_vertices,
                                                            perm))),
                                              copy=True)
            # Check for a permutation that leaves the graph unchanged
            # (only way to keep the edge list of the same length)
            intersection = copy.deepcopy(self.check_graph)
            intersection.remove_edges_from(e for e in self.check_graph.edges()
                                           if e not in permuted_graph.edges())
            check = nx.algorithms.isomorphism.DiGraphMatcher(self.check_graph,
                                                             intersection,
                                                             node_match=op_nm,
                                                             edge_match=anom_em)
            if check.is_isomorphic():
                permutations.append(dict(list(zip(perm_vertices, perm))))
        return permutations

    def has_anom_props_linked_sign(self):
        """Return True if there is a minus sign associated to anom props.

        Anomalous propagators departing to higher vertices introduce a sign
        factor if a normal propagator is going to an even higher vertex, as it
        departs from the canonical representation used for numerator
        extraction.

        Returns:
            bool: The presence of the sign factor.

        """
        counter = 0
        for vertex in self.graph:
            for prop_a in self.graph.out_edges(vertex,
                                               keys=True,
                                               data='anomalous'):
                # Check that prop_a is anomalous and not a self-contraction
                if prop_a[3] and prop_a[1] != prop_a[0]:
                    for prop_b in self.graph.out_edges(vertex,
                                                       keys=True,
                                                       data='anomalous'):
                        # Check for a normal prop going to a higher vertex
                        if prop_b[1] > prop_a[1] and not prop_b[3]:
                            counter += 1
        # True if overall factor is odd
        return counter % 2 == 1

    def has_anom_non_selfcontracted_props(self):
        """Return True if the diagram has anomalous propagators.

        Returns:
            bool: The presence of anomalous propagators.

        """
        for prop in self.graph.edges(keys=True, data='anomalous'):
            if prop[3] and prop[1] != prop[0]:
                return True
        return False

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

    @property
    def check_graph(self):
        """Return a graph that can be used for topological equivalence checks.

        Lazy-initialized to reduce memory and CPU costs as this operation
        requires a deep copy.

        Returns:
            NetworkX MultiDiGraph: The graph with doubled anomalous props.

        """
        if self._check_graph is None:
            self._check_graph = adg.diag.create_checkable_diagram(self.graph)
        return self._check_graph

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
