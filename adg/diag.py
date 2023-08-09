"""Routines and class for all types of diagrams, inherited by others."""

from builtins import range
from builtins import object, str
from adg.tools import reversed_enumerate

import copy
import numpy
import networkx as nx


def no_trace(matrices):
    """Select matrices with full 0 diagonal.

    Args:
        matrices (list): A list of adjacency matrices.

    Returns:
        list: The adjacency matrices without non-zero diagonal elements.

    >>> test_matrices = [[[0, 1, 2], [2, 0, 1], [5, 2, 0]], \
    [[2, 2, 2], [1, 2, 3], [0, 0, 0]], \
    [[0, 1, 3], [2, 0, 8], [2, 1, 0]]]
    >>> no_trace(test_matrices)
    [[[0, 1, 2], [2, 0, 1], [5, 2, 0]], [[0, 1, 3], [2, 0, 8], [2, 1, 0]]]

    """
    traceless_matrices = []
    for matrix in matrices:
        for ind_i, line in enumerate(matrix):
            if line[ind_i] != 0:
                break
        else:
            traceless_matrices.append(matrix)
    return traceless_matrices


def check_vertex_degree(matrices, three_body_use, nbody_max_observable,
                        canonical_only, vertex_id):
    """Check the degree of a specific vertex in a set of matrices.

    Args:
        matrices (list): Adjacency matrices.
        three_body_use (bool): ``True`` if one uses three-body forces.
        nbody_max_observable (int): Maximum body number for the observable.
        canonical_only (bool): ``True`` if one draws only canonical diagrams.
        vertex_id (int): The position of the studied vertex.

    >>> test_matrices = [numpy.array([[0, 1, 2], [1, 0, 1], [0, 2, 0]]), \
        numpy.array([[2, 0, 2], [1, 2, 3], [1, 0, 0]]), \
        numpy.array([[0, 1, 3], [2, 0, 8], [2, 1, 0]])]
    >>> check_vertex_degree(test_matrices, True, 3, False, 0)
    >>> test_matrices #doctest: +NORMALIZE_WHITESPACE
    [array([[0, 1, 2], [1, 0, 1], [0, 2, 0]]),
     array([[2, 0, 2], [1, 2, 3], [1, 0, 0]])]
    >>> check_vertex_degree(test_matrices, False, 2, False, 0)
    >>> test_matrices #doctest: +NORMALIZE_WHITESPACE
    [array([[0, 1, 2], [1, 0, 1], [0, 2, 0]])]

    """
    authorized_deg = [4]
    if three_body_use:
        authorized_deg.append(6)
    if not canonical_only or vertex_id == 0:
        authorized_deg.append(2)
    authorized_deg = tuple(authorized_deg)

    for i_mat, matrix in reversed_enumerate(matrices):
        vertex_degree = sum(matrix[index][vertex_id] + matrix[vertex_id][index]
                            for index in list(range(matrix.shape[0])))
        vertex_degree -= matrix[vertex_id][vertex_id]

        if (vertex_id != 0 and vertex_degree not in authorized_deg) \
                or (vertex_id == 0 and vertex_degree > 2*nbody_max_observable):
            del matrices[i_mat]


def topologically_distinct_diagrams(diagrams):
    """Return a list of diagrams all topologically distinct.

    Args:
        diagrams (list): The Diagrams of interest.

    Returns:
        list: Topologically unique diagrams.

    """
    import adg.tsd
    iso = nx.algorithms.isomorphism
    op_nm = iso.categorical_node_match('operator', False)
    anom_em = iso.categorical_multiedge_match('anomalous', False)
    for i_diag, diag in reversed_enumerate(diagrams):
        graph = diag.graph
        diag_io_degrees = diag.io_degrees
        for i_comp_diag, comp_diag in reversed_enumerate(diagrams[:i_diag]):
            if diag_io_degrees == comp_diag.io_degrees:
                # Check anomalous character of props for PBMBPT
                if isinstance(diag, adg.pbmbpt.ProjectedBmbptDiagram):
                    doubled_graph = create_checkable_diagram(graph)
                    doubled_comp_diag = create_checkable_diagram(comp_diag.graph)
                    matcher = iso.DiGraphMatcher(doubled_graph,
                                                 doubled_comp_diag,
                                                 node_match=op_nm,
                                                 edge_match=anom_em)
                # Check for topologically equivalent diags considering vertex
                # properties but not edge attributes, i.e. anomalous character
                else:
                    matcher = iso.DiGraphMatcher(graph,
                                                 comp_diag.graph,
                                                 node_match=op_nm)
                if matcher.is_isomorphic():
                    # Store the set of permutations to recover the original TSD
                    if isinstance(diag, adg.tsd.TimeStructureDiagram):
                        diag.perms.update(
                            update_permutations(comp_diag.perms,
                                                comp_diag.tags[0],
                                                matcher.mapping)
                            )
                    diag.tags += comp_diag.tags
                    del diagrams[i_comp_diag]
                    break
    return diagrams


def update_permutations(comp_graph_perms, comp_graph_tag, mapping):
    """Update permutations associated to the BMBPT diags for a shared TSD.

    Args:
        comp_graph_perms (dict): Permutations to be updated.
        comp_graph_tag (int): The tag associated to the TSD configuration.
        mapping (dict): permutations to go from previous ref TSD to new one.

    """
    identity = {key: key for key in comp_graph_perms[comp_graph_tag]}
    # Do permutations only when necessary
    if mapping != identity:
        for graph_id in comp_graph_perms:
            # Create a dummy dictionary to avoid overwriting some nodes
            dummy_nodes = copy.deepcopy(comp_graph_perms[graph_id])
            # Permute the nodes according to the new mapping
            for node in comp_graph_perms[graph_id]:
                comp_graph_perms[graph_id][node] = dummy_nodes[mapping[node]]

    return comp_graph_perms


def create_checkable_diagram(pbmbpt_graph):
    """Return a graph with anomalous props going both ways for topo check.

    Args:
        pbmbpt_graph (NetworkX MultiDiGraph): The graph to be copied.

    Returns:
        NetworkX MultiDiGraph: Graph with double the anomalous props.

    """
    doubled_graph = copy.deepcopy(pbmbpt_graph)
    props_to_add = [(prop[1], prop[0]) for prop
                    in doubled_graph.edges(keys=True, data='anomalous')
                    if prop[3] and not prop[0] == prop[1]]
    doubled_graph.add_edges_from(props_to_add, anomalous=True, weight=1)
    return doubled_graph


def label_vertices(graphs_list, theory_type, switch_flag):
    """Account for different status of vertices in operator diagrams.

    Args:
        graphs_list (list): The Diagrams of interest.
        theory_type (str): The name of the theory of interest.
        switch_flag (int): When to switch A and B operators for BIMSRG.

    """
    if theory_type != 'BIMSRG':
        for graph in graphs_list:
            nx.set_node_attributes(graph, False, 'operator')
            if theory_type in ("BMBPT", "PBMBPT"):
                graph.nodes[0]['operator'] = True
    else:
        for idx, graph in enumerate(graphs_list):
            nx.set_node_attributes(graph, '', 'operator')
            if idx < switch_flag:
                graph.nodes[1]['operator'] = 'B'
                graph.nodes[2]['operator'] = 'A'
            else:
                graph.nodes[1]['operator'] = 'A'
                graph.nodes[2]['operator'] = 'B'


def feynmf_generator(graph, theory_type, diagram_name):
    """Generate the feynmanmp instructions corresponding to the diagram.

    Args:
        graph (NetworkX MultiDiGraph): The graph of interest.
        theory_type (str): The name of the theory of interest.
        diagram_name (str): The name of the studied diagram.

    """
    p_order = graph.number_of_nodes()
    diag_size = 20*p_order

    theories = ["MBPT", "BMBPT", "PBMBPT", "BIMSRG"]
    prop_types = ["half_prop", "prop_pm", "prop_pm", "half_prop"]
    propa = prop_types[theories.index(theory_type)]

    fmf_file = open(diagram_name + ".tex", 'w')
    fmf_file.write("\\parbox{40pt}{\\begin{fmffile}{%s}\n" % diagram_name
                   + "\\begin{fmfgraph*}(40,%i)\n" % diag_size)

    # Define the appropriate line propagator_style
    fmf_file.write(propagator_style(propa))
    if theory_type == "PBMBPT":
        fmf_file.write(propagator_style("prop_mm"))

    # Set the position of the vertices
    if theory_type == "BIMSRG":
        bimsrg_diagram_internals(graph, fmf_file, propa)

    else:
        fmf_file.write(vertex_positions(graph, p_order))

        # Special config for self-contraction
        if theory_type == "PBMBPT":
            fmf_file.write(self_contractions(graph))

        # Loop over all elements of the graph to draw associated propagators
        for vert_i in graph:
            for vert_j in list(graph.nodes())[vert_i+1:]:
                props_to_draw = [edge for edge in graph.edges([vert_i, vert_j],
                                                              data=True,
                                                              keys=True)
                                 if edge[1] in (vert_i, vert_j)
                                 and edge[0] != edge[1]]
                # Set the list of propagators directions to use
                props_dir = prop_directions(vert_j - vert_i, len(props_to_draw))
                # Draw the diagrams
                key = 0
                # Start with props going down, used in MBPT only
                for _ in (p for p in props_to_draw
                          if p[1] < p[0]
                          and not ('anomalous' in p[3] and p[3]['anomalous'])):
                    fmf_file.write("\\fmf{%s%s}{v%i,v%i}\n"
                                   % (propa, props_dir[key], vert_j, vert_i))
                    key += 1
                # Reinitialise the drawing configuration as we change direction
                key = 0
                for _ in (p for p in props_to_draw
                          if ('anomalous' in p[3] and p[3]['anomalous'])):
                    fmf_file.write("\\fmf{prop_mm%s}{v%i,v%i}\n"
                                   % (props_dir[key], vert_i, vert_j))
                    key += 1
                for _ in (p for p in props_to_draw
                          if p[0] < p[1]
                          and not ('anomalous' in p[3] and p[3]['anomalous'])):
                    fmf_file.write("\\fmf{%s%s}{v%i,v%i}\n"
                                   % (propa, props_dir[key], vert_i, vert_j))
                    key += 1

    fmf_file.write("\\end{fmfgraph*}\n\\end{fmffile}}\n")
    fmf_file.close()


def prop_directions(vert_distance, nb_props):
    """Return a list of possible propagators directions.

    Args:
        vert_distance (int): Distance between the two connected vertices.
        nb_props (int): Number of propagators to be drawn.

    Returns:
        list: Propagators directions stored as strings.

    """
    if nb_props < 7:
        directions = [",right=0.9", ",right=0.75", ",right=0.6", ",right=0.5",
                      "", ",left=0.5", ",left=0.6", ",left=0.75", ",left=0.9"]

        if vert_distance != 1:
            props_dir = directions[:3] + directions[-3:]
        else:
            props_dir = directions[:2] + directions[3:6] + directions[-2:]
            if nb_props % 2 != 1:
                props_dir = props_dir[:3] + props_dir[-3:]
            else:
                props_dir = props_dir[1:]
        if nb_props < 5:
            props_dir = props_dir[1:-1]
            if nb_props < 3:
                props_dir = props_dir[1:-1]

    elif vert_distance == 1:
        directions = [",right=0.%i" % angle for angle in range(90, 0, -10)] \
            + [",left=0.%i" % angle for angle in range(10, 100, 10)]

        props_dir = directions[nb_props//4:-(nb_props//4)]

    return props_dir


def propagator_style(prop_type):
    """Return the FeynMF definition for the appropriate propagator type.

    Args:
        prop_type (str): The type of propagators used in the diagram.

    Returns:
        str: The FeynMF definition for the propagator style used.

    """
    line_styles = {}

    line_styles['prop_pm'] = "\\fmfcmd{style_def prop_pm expr p =\n" \
        + "draw_plain p;\nshrink(.7);\n" \
        + "\tcfill (marrow (p, .25));\n" \
        + "\tcfill (marrow (p, .75))\n" \
        + "endshrink;\nenddef;}\n"

    line_styles['prop_mm'] = "\\fmfcmd{style_def prop_mm expr p =\n" \
        + "draw_plain p;\nshrink(.7);\n" \
        + "\tcfill (marrow (p, .75));\n" \
        + "\tcfill (marrow (reverse p, .75))\n" \
        + "endshrink;\nenddef;}\n"

    line_styles['half_prop'] = "\\fmfcmd{style_def half_prop expr p =\n" \
        + "draw_plain p;\nshrink(.7);\n" \
        + "\tcfill (marrow (p, .5))\n" \
        + "endshrink;\nenddef;}\n"

    return line_styles[prop_type]


def vertex_positions(graph, order):
    """Return the positions of the graph's vertices.

    Args:
        graph (NetworkX MultiDiGraph): The graph of interest.
        order (int): The perturbative order of the graph.

    Returns:
        str: The FeynMP instructions for positioning the vertices.

    """
    positions = "\\fmftop{v%i}\\fmfbottom{v0}\n" % (order-1)
    for vert in range(order-1):
        positions += "\\fmf{phantom}{v%i,v%i}\n" % (vert, (vert+1)) \
            + ("\\fmfv{d.shape=square,d.filled=full,d.size=3thick"
               if 'operator' in graph.nodes[vert] and graph.nodes[vert]['operator']
               else "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick") \
            + "}{v%i}\n" % vert
    positions += "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" \
        % (order-1) + "\\fmffreeze\n"
    return positions


def self_contractions(graph):
    """Return the instructions for drawing the graph's self-contractions.

    Args:
        graph (NetworkX MultiDiGraph): The graph being drawn.

    Returns:
        str: FeynMF instructions for drawing the self-contractions.

    """
    instructions = ""
    # Check for self-contractions before going further
    if [nx.selfloop_edges(graph)]:
        instructions += propagator_style("half_prop")
        for vert in graph:
            props_to_draw = [edge for edge
                             in nx.selfloop_edges(graph, data=True, keys=True)
                             if edge[0] == vert]
            positions = ["15pt", "-15pt"]
            key = 0
            for prop in props_to_draw:
                if prop[3]['anomalous']:
                    a_name = "a%i%i" % (vert, key)
                    instructions += ("\\fmfv{}{%s}\n" % a_name
                                     + "\\fmffixed{(%s,0)}{v%i,%s}\n"
                                     % (positions[key], vert, a_name)
                                     + "\\fmf{half_prop,right}{%s,v%i}\n"
                                     % (a_name, vert)
                                     + "\\fmf{half_prop,left}{%s,v%i}\n"
                                     % (a_name, vert))
                    key += 1
        instructions += "\\fmffreeze\n"
    return instructions


def bimsrg_diagram_internals(graph, fmf_file, prop_type):
    """Write to file the vertices and propagators of BIMSRG diagrams.

    Args:
        graph (NetworkX MultiDiGraph): The graph to be drawn.
        fmf_file (file): The FeymanMF file to be written to.
        prop_type (str): The FeymanMF type for drawing the propagators.

    """
    nbs_out_edges = (sum(1 for p in graph.in_edges(3, keys=True) if p[0] == 1),
                     sum(1 for p in graph.in_edges(3, keys=True) if p[0] == 2))
    nbs_in_edges = (sum(1 for p in graph.out_edges(0, keys=True) if p[1] == 1),
                    sum(1 for p in graph.out_edges(0, keys=True) if p[1] == 2))

    nb_out_edges = len(graph.in_edges(3, keys=True))
    nb_in_edges = len(graph.out_edges(0, keys=True))

    nb_top_vertices = nb_out_edges if nb_out_edges == 1 \
        else (max(2*nbs_out_edges[0], 2*nbs_out_edges[1]) + 1)
    nb_bot_vertices = nb_in_edges if nb_in_edges == 1 \
        else (max(2*nbs_in_edges[0], 2*nbs_in_edges[1]) + 1)

    positions = "\\fmfstraight\n" \
        + "\\fmftopn{t}{%i}" % nb_top_vertices \
        + "\\fmfbottomn{b}{%i}\n" % nb_bot_vertices \
        + "\\fmf{phantom}{b%i,v1}\n" % (nb_bot_vertices//2 + 1) \
        + "\\fmf{phantom}{v1,v2}\n" \
        + "\\fmf{phantom}{v2,t%i}\n" % (nb_top_vertices//2 + 1) \
        + "\\fmfv{d.shape=circle,d.filled=%s,d.size=3thick}{v1}\n" \
        % ('full' if graph.nodes[1]['operator'] == 'A' else 'empty') \
        + "\\fmfv{d.shape=circle,d.filled=%s,d.size=3thick}{v2}\n" \
        % ('full' if graph.nodes[2]['operator'] == 'A' else 'empty') \
        + "\\fmffreeze\n"
    fmf_file.write(positions)

    # Internal lines
    nb_props = sum(1 for edge in graph.edges(1, keys=True) if edge[1] == 2)
    # Set the list of propagators directions to use
    props_dir = prop_directions(1, nb_props)
    # Draw the propagators
    for idx in range(nb_props):
        fmf_file.write("\\fmf{%s%s}{v1,v2}\n" % (prop_type, props_dir[idx]))

    # Incoming external line
    for vertex in (1, 2):
        nb_props = sum(1 for edge in graph.edges(0, keys=True)
                       if edge[1] == vertex)
        if (nb_bot_vertices == 1) and (vertex == 1):
            orientation = ""
        else:
            orientation = ",left=0.4" if vertex == 2 else ",right=0.3"
        # Draw the propagators
        for key in range(nb_props):
            fmf_file.write("\\fmf{%s%s}{b%i,v%i}\n"
                           % (prop_type,
                              orientation,
                              key+1 if vertex == 2 else nb_bot_vertices - key,
                              vertex))

    # Outgoing external lines
    for vertex in (1, 2):
        nb_props = sum(1 for edge in graph.in_edges(3, keys=True)
                       if edge[0] == vertex)
        if (nb_top_vertices == 1) and (vertex == 2):
            orientation = ""
        else:
            orientation = ",right=0.4" if vertex == 1 else ",left=0.3"
        # Draw the propagators
        for key in range(nb_props):
            fmf_file.write("\\fmf{%s%s}{v%i,t%i}\n"
                           % (prop_type,
                              orientation,
                              vertex,
                              key+1 if vertex == 2 else nb_top_vertices - key))


def draw_diagram(directory, result_file, diagram_index, diag_type):
    """Copy the diagram feynmanmp instructions in the result file.

    Args:
        directory (str): The path to the output folder.
        result_file (file): The LaTeX ouput file of the program.
        diagram_index (str): The number associated to the diagram.
        diag_type (str): The type of diagram used here.

    """
    with open(directory+"/Diagrams/%s_%s.tex"
              % (diag_type, diagram_index)) as diag_file:
        result_file.write(diag_file.read())


def to_skeleton(graph):
    """Return the bare skeleton of a graph, i.e. only non-redundant links.

    Args:
        graph (NetworkX MultiDiGraph): The graph to be turned into a skeleton.

    Returns:
        NetworkX MultiDiGraph: The skeleton of the initial graph.

    """
    for vertex_a in graph:
        for vertex_b in graph:
            while graph.number_of_edges(vertex_a, vertex_b) > 1:
                graph.remove_edge(vertex_a, vertex_b)
            if len(list(nx.all_simple_paths(graph, vertex_a, vertex_b))) > 1:
                while len(nx.shortest_path(graph, vertex_a, vertex_b)) == 2:
                    graph.remove_edge(vertex_a, vertex_b)
    return graph


def extract_denom(start_graph, subgraph):
    """Extract the appropriate denominator using the subgraph rule.

    Args:
        start_graph (NetworkX MultiDiGraph): The studied graph.
        subgraph (NetworkX MultiDiGraph): The subgraph used for this particular
            denominator factor.

    Returns:
        str: The denominator factor for this subgraph.

    """
    denomin = r"\epsilon^{" \
        + "".join("%s"
                  % prop[3]['qp_state']
                  for prop
                  in start_graph.out_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(prop[0], prop[1], prop[2])
                  and not ('anomalous' in prop[3] and prop[3]['anomalous'])) \
        + "}_{" \
        + "".join("%s"
                  % prop[3]['qp_state']
                  for prop
                  in start_graph.in_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(prop[0], prop[1], prop[2])
                  and not ('anomalous' in prop[3] and prop[3]['anomalous'])) \
        + "".join("%s"
                  % prop[3]['qp_state']
                  for prop
                  in start_graph.in_edges(subgraph, keys=True, data=True)
                  if subgraph.has_edge(prop[0], prop[1], prop[2])
                  and ('anomalous' in prop[3] and prop[3]['anomalous'])) \
        + "".join("%s"
                  % (prop[3]['qp_state'].split("}")[0] + "}")
                  for prop
                  in start_graph.in_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(prop[0], prop[1], prop[2])
                  and ('anomalous' in prop[3] and prop[3]['anomalous'])) \
        + "".join("%s"
                  % (prop[3]['qp_state'].split("}")[1] + "}")
                  for prop
                  in start_graph.out_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(prop[0], prop[1], prop[2])
                  and ('anomalous' in prop[3] and prop[3]['anomalous'])) \
        + "}"
    return denomin


def print_adj_matrices(directory, diagrams):
    """Print a computer-readable file with the diagrams' adjacency matrices.

    Args:
        directory (str): The path to the output directory.
        diagrams (list): All the diagrams.

    """
    with open(directory+"/adjacency_matrices.txt", "w") as mat_file:
        for idx, diagram in enumerate(diagrams):
            mat_file.write("Diagram n: %i\n" % (idx + 1))
            numpy.savetxt(mat_file,
                          nx.to_numpy_array(diagram.graph, dtype=int),
                          fmt='%d')
            mat_file.write("\n")


class Diagram(object):
    """Describes a diagram with its related properties.

    Args:
        nx_graph (NetworkX MultiDiGraph): The graph of interest.

    """

    __slots__ = ('graph', 'unsort_degrees', 'degrees', 'unsort_io_degrees',
                 'io_degrees', 'max_degree', 'tags')

    def __init__(self, nx_graph):
        """Generate a Diagram object starting from the NetworkX graph."""
        self.graph = nx_graph
        """NetworkX MultiDiGraph: The actual graph."""
        self.unsort_degrees = tuple(nx_graph.degree(node) for node in nx_graph)
        """tuple: The degrees of the graph vertices"""
        self.degrees = tuple(sorted(self.unsort_degrees))
        """tuple: The ascendingly sorted degrees of the graph vertices."""
        self.unsort_io_degrees = tuple((nx_graph.in_degree(node),
                                        nx_graph.out_degree(node))
                                       for node in nx_graph)
        """tuple: The list of in- and out-degrees for each vertex of the graph,
        stored in a (in, out) tuple."""
        self.io_degrees = tuple(sorted(self.unsort_io_degrees))
        """tuple: The sorted version of unsort_io_degrees."""
        self.max_degree = self.degrees[-1]
        """int: The maximal degree of a vertex in the graph."""
        self.tags = [0]
        """list: The tag numbers associated to a diagram."""

    def write_graph(self, latex_file, directory, write_time):
        """Write the graph of the diagram to the LaTeX file.

        Args:
            latex_file (file): The LaTeX ouput file of the program.
            directory (str): Path to the result folder.
            write_time (bool): (Here to emulate polymorphism).

        """
        latex_file.write('\n\\begin{center}\n')
        draw_diagram(directory, latex_file, str(self.tags[0]), 'diag')
        latex_file.write('\n\\end{center}\n\n')
