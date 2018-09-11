"""Routines and class for all types of diagrams, inherited by others."""

import numpy
import networkx as nx


def no_trace(matrices):
    """Select matrices with full 0 diagonal.

    Args:
        matrices (list): A list of adjacency matrices.

    Returns:
        (list): The adjacency matrices without non-zero diagonal elements.

    >>> test_matrices = [[[0, 1, 2], [2, 0, 1], [5, 2, 0]], \
    [[2, 2, 2], [1, 2, 3], [0, 0, 0]], \
    [[0, 1, 3], [2, 0, 8], [2, 1, 0]]]
    >>> no_trace(test_matrices)
    [[[0, 1, 2], [2, 0, 1], [5, 2, 0]], [[0, 1, 3], [2, 0, 8], [2, 1, 0]]]
    >>> no_trace()
    Traceback (most recent call last):
      File "/usr/lib/python2.7/doctest.py", line 1315, in __run
        compileflags, 1) in test.globs
      File "<doctest __main__.no_trace[4]>", line 1, in <module>
        no_trace()
    TypeError: no_trace() takes exactly 1 argument (0 given)

    """
    traceless_matrices = []
    for matrix in matrices:
        test_traceless = True
        for ind_i, line in enumerate(matrix):
            if line[ind_i] != 0:
                test_traceless = False
                break
        if test_traceless:
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
     array[[2, 0, 2], [1, 2, 3], [1, 0, 0]])]
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

    for i_mat in xrange(len(matrices)-1, -1, -1):
        matrix = matrices[i_mat]
        vertex_degree = sum(matrix[index][vertex_id] + matrix[vertex_id][index]
                            for index in range(matrix.shape[0]))
        vertex_degree -= matrix[vertex_id][vertex_id]

        if (vertex_id != 0 and vertex_degree not in authorized_deg) \
                or (vertex_id == 0 and vertex_degree > 2*nbody_max_observable):
            del matrices[i_mat]


def topologically_distinct_diagrams(diagrams):
    """Return a list of diagrams all topologically distinct.

    Args:
        diagrams (list): The Diagrams of interest.

    Returns:
        (list): Topologically unique diagrams.

    """
    import adg.tsd
    iso = nx.algorithms.isomorphism
    op_nm = iso.categorical_node_match('operator', False)
    for i_diag in xrange(len(diagrams)-1, -1, -1):
        graph = diagrams[i_diag].graph
        diag_io_degrees = diagrams[i_diag].io_degrees
        for i_comp_diag in xrange(i_diag+1, len(diagrams), 1):
            if diag_io_degrees == diagrams[i_comp_diag].io_degrees:
                matcher = iso.DiGraphMatcher(graph,
                                             diagrams[i_comp_diag].graph,
                                             node_match=op_nm)
                if matcher.is_isomorphic():
                    diagrams[i_diag].tags += diagrams[i_comp_diag].tags
                    if isinstance(diagrams[i_diag],
                                  adg.tsd.TimeStructureDiagram):
                        diagrams[i_diag].perms.update(
                            diagrams[i_comp_diag].perms)
                        diagrams[i_diag].perms[diagrams[i_comp_diag].tags[0]] \
                            = matcher.mapping
                    del diagrams[i_comp_diag]
                    break
    return diagrams


def label_vertices(graphs_list, theory_type):
    """Account for different status of vertices in operator diagrams.

    Args:
        graphs_list (list): The Diagrams of interest.
        theory_type (str): The name of the theory of interest.

    """
    for graph in graphs_list:
        for node in graph:
            graph.node[node]['operator'] = False
        if theory_type == "BMBPT" or "PBMBPT":
            graph.node[0]['operator'] = True


# def feynmf_generator(graph, theory_type, diagram_name):
#     """Generate the feynmanmp instructions corresponding to the diagram.
#
#     Args:
#         graph (NetworkX MultiDiGraph): The graph of interest.
#         theory_type (str): The name of the theory of interest.
#         diagram_name (str): The name of the studied diagram.
#
#     """
#     p_order = graph.number_of_nodes()
#     diag_size = 20*p_order
#
#     theories = ["MBPT", "BMBPT", "PBMBPT"]
#     prop_types = ["half_prop", "prop_pm", "prop_pm"]
#     propa = prop_types[theories.index(theory_type)]
#
#     fmf_file = open(diagram_name + ".tex", 'w')
#     fmf_file.write("\\parbox{%ipt}{\\begin{fmffile}{%s}\n" % (diag_size,
#                                                               diagram_name)
#                    + "\\begin{fmfgraph*}(%i,%i)\n" % (diag_size, diag_size))
#
#     # Define the appropriate line propagator_style
#     fmf_file.write(propagator_style(propa))
#     if theory_type == "PBMBPT":
#         fmf_file.write(propagator_style("prop_mm"))
#
#     # Set the position of the vertices
#     fmf_file.write(vertex_positions(graph, p_order))
#
#     # Loop over all elements of the graph to draw associated propagators
#     for vert_i in graph:
#         for vert_j in graph:
#             props_left_to_draw = graph.number_of_edges(vert_i, vert_j)
#             # Special config for consecutive vertices
#             if (props_left_to_draw % 2 == 1) and (abs(vert_i-vert_j) == 1):
#                 fmf_file.write("\\fmf{%s" % propa)
#                 # Check for specific MBPT configuration
#                 if graph.number_of_edges(vert_j, vert_i) == 1:
#                     fmf_file.write(",right=0.5")
#                 fmf_file.write("}{v%i,v%i}\n" % (vert_i, vert_j))
#                 props_left_to_draw -= 1
#             while props_left_to_draw > 0:
#                 fmf_file.write("\\fmf{%s," % propa)
#                 fmf_file.write("right=" if props_left_to_draw % 2 == 1
#                                else "left=")
#                 if props_left_to_draw in (5, 6):
#                     fmf_file.write("0.9")
#                 elif props_left_to_draw in (3, 4) \
#                         or (props_left_to_draw == 1
#                             and graph.number_of_edges(vert_j, vert_i) == 2):
#                     fmf_file.write("0.75")
#                 elif props_left_to_draw in (1, 2):
#                     fmf_file.write("0.5" if abs(vert_i-vert_j) == 1 else "0.6")
#                 fmf_file.write("}{v%i,v%i}\n" % (vert_i, vert_j))
#                 props_left_to_draw -= 1
#         # Special config for self-contraction
#         props_left_to_draw = len(list(edge for edge
#                                       in nx.selfloop_edges(graph,
#                                                            data=True,
#                                                            keys=True)
#                                       if edge[0] == vert_i))
#         while props_left_to_draw > 0:
#             if props_left_to_draw > 1:
#                 fmf_file.write("\\fmf{prop_mm,left=45}{v%i,v%i}\n"
#                                % (vert_i, vert_i))
#             else:
#                 fmf_file.write("\\fmf{prop_mm,right=45}{v%i,v%i}\n"
#                                % (vert_i, vert_i))
#             props_left_to_draw -= 1
#     fmf_file.write("\\end{fmfgraph*}\n\\end{fmffile}}\n")
#     fmf_file.close()


def feynmf_generator(graph, theory_type, diagram_name):
    """Generate the feynmanmp instructions corresponding to the diagram.

    Args:
        graph (NetworkX MultiDiGraph): The graph of interest.
        theory_type (str): The name of the theory of interest.
        diagram_name (str): The name of the studied diagram.

    """
    p_order = graph.number_of_nodes()
    diag_size = 20*p_order

    theories = ["MBPT", "BMBPT", "PBMBPT"]
    prop_types = ["half_prop", "prop_pm", "prop_pm"]
    propa = prop_types[theories.index(theory_type)]

    fmf_file = open(diagram_name + ".tex", 'w')
    fmf_file.write("\\parbox{%ipt}{\\begin{fmffile}{%s}\n" % (diag_size,
                                                              diagram_name)
                   + "\\begin{fmfgraph*}(%i,%i)\n" % (diag_size, diag_size))

    # Define the appropriate line propagator_style
    fmf_file.write(propagator_style(propa))
    if theory_type == "PBMBPT":
        fmf_file.write(propagator_style("prop_mm"))

    # Set the position of the vertices
    fmf_file.write(vertex_positions(graph, p_order))

    directions = [",left=0.9", ",left=0.75", ",left=0.6", ",left=0.5", "",
                  ",right=0.5", ",right=0.6", ",right=0.75", ",right=0.9"]

    # Loop over all elements of the graph to draw associated propagators
    for vert_i in graph:
        for vert_j in list(graph.nodes())[vert_i+1:]:
            props_to_draw = [edge for edge in graph.edges(data=True, keys=True)
                             if edge[0] in (vert_i, vert_j)
                             and edge[1] in (vert_i, vert_j)
                             and edge[0] != edge[1]]
            # Set the list of propagators directions to use
            if vert_j - vert_i != 1:
                props_dir = directions[:3] + directions[-3:]
            else:
                props_dir = directions[:2] + directions[3:6] + directions[-2:]
                if len(props_to_draw) % 2 != 1:
                    props_dir = props_dir[:3] + props_dir[-3:]
                else:
                    props_dir = props_dir[1:]
            if len(props_to_draw) < 5:
                props_dir = props_dir[1:-1]
                if len(props_to_draw) < 3:
                    props_dir = props_dir[1:-1]
            # Draw the diagrams
            key = 0
            for prop in props_to_draw:
                if prop[1] < prop[0] \
                        and not graph[prop[0]][prop[1]][prop[2]]['anomalous']:
                    fmf_file.write("\\fmf{%s%s}{v%i,v%i}\n"
                                   % (propa, props_dir[key], vert_j, vert_i))
                    key += 1
            for prop in props_to_draw:
                if prop[0] < prop[1] \
                        and not graph[prop[0]][prop[1]][prop[2]]['anomalous']:
                    fmf_file.write("\\fmf{%s%s}{v%i,v%i}\n"
                                   % (propa, props_dir[key], vert_i, vert_j))
                    key += 1
            for prop in props_to_draw:
                if graph[prop[0]][prop[1]][prop[2]]['anomalous']:
                    fmf_file.write("\\fmf{prop_mm%s}{v%i,v%i}\n"
                                   % (props_dir[key], vert_i, vert_j))
                    key += 1

        # Special config for self-contraction
        props_to_draw = [edge for edge
                         in nx.selfloop_edges(graph, data=True, keys=True)
                         if edge[0] == vert_i]
        angle = [",right=45", ",left=45"]
        key = 0
        for prop in props_to_draw:
            if graph[prop[0]][prop[1]][prop[2]]['anomalous']:
                fmf_file.write("\\fmf{prop_mm%s}{v%i,v%i}\n"
                               % (angle[key], vert_i, vert_i))
                key += 1
    fmf_file.write("\\end{fmfgraph*}\n\\end{fmffile}}\n")
    fmf_file.close()


def propagator_style(prop_type):
    """Return the FeynMF definition for the appropriate propagator type.

    Args:
        prop_type (str): The type of propagators used in the diagram.

    Returns:
        (str): The FeynMF definition for the propagator style used.

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
        (str): The FeynMP instructions for positioning the vertices.

    """
    positions = "\\fmftop{v%i}\\fmfbottom{v0}\n" % (order-1)
    for vert in xrange(order-1):
        positions += "\\fmf{phantom}{v%i,v%i}\n" % (vert, (vert+1)) \
            + ("\\fmfv{d.shape=square,d.filled=full,d.size=3thick"
               if graph.node[vert]['operator']
               else "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick") \
            + "}{v%i}\n" % vert
    positions += "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" \
        % (order-1) + "\\fmffreeze\n"
    return positions


def draw_diagram(directory, result_file, diagram_index, diag_type):
    """Copy the diagram feynmanmp instructions in the result file.

    Args:
        directory (str): The path to the output folder.
        result_file (file): The LaTeX ouput file of the program.
        diagram_index (str): The number associated to the diagram.
        diag_type (str): The type of diagram used here.

    """
    diag_file = open(directory+"/Diagrams/%s_%s.tex" % (diag_type,
                                                        diagram_index))
    result_file.write(diag_file.read())
    diag_file.close()


def to_skeleton(graph):
    """Return the bare skeleton of a graph, i.e. only non-redundant links.

    Args:
        graph (NetworkX MultiDiGraph): The graph to be turned into a skeleton.

    Returns:
        (NetworkX MultiDiGraph): The skeleton of the initial graph.

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
        (str): The denominator factor for this subgraph.

    """
    denomin = r"\epsilon^{" \
        + "".join("%s"
                  % propa[3]['qp_state']
                  for propa
                  in start_graph.out_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2])) \
        + "}_{" \
        + "".join("%s"
                  % propa[3]['qp_state']
                  for propa
                  in start_graph.in_edges(subgraph, keys=True, data=True)
                  if not subgraph.has_edge(propa[0], propa[1], propa[2])) \
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
            numpy.savetxt(mat_file, diagram.adjacency_mat, fmt='%d')
            mat_file.write("\n")


class Diagram(object):
    """Describes a diagram with its related properties.

    Attributes:
        graph (NetworkX MultiDiGraph): The actual graph.
        unsorted_degrees (tuple): The degrees of the graph vertices
        degrees (tuple): The ascendingly sorted degrees of the graph vertices.
        unsort_io_degrees (tuple): The list of in- and out-degrees for each
            vertex of the graph, stored in a (in, out) tuple.
        io_degrees (tuple): The sorted version of unsort_io_degrees.
        max_degree (int): The maximal degree of a vertex in the graph.
        tags (list): The tag numbers associated to a diagram.
        adjacency_mat (NumPy array): The adjacency matrix of the graph.

    """

    def __init__(self, nx_graph):
        """Generate a Diagram object starting from the NetworkX graph.

        Args:
            nx_graph (NetworkX MultiDiGraph): The graph of interest.

        """
        self.graph = nx_graph
        self.unsort_degrees = tuple(nx_graph.degree(node) for node in nx_graph)
        self.degrees = sorted(self.unsort_degrees)
        self.unsort_io_degrees = tuple((nx_graph.in_degree(node),
                                        nx_graph.out_degree(node))
                                       for node in nx_graph)
        self.io_degrees = sorted(self.unsort_io_degrees)
        self.max_degree = self.degrees[-1]
        self.tags = [0]
        self.adjacency_mat = nx.to_numpy_matrix(self.graph, dtype=int)

    def write_graph(self, latex_file, directory, write_time):
        """Write the graph of the diagram to the LaTeX file.

        Args:
            latex_file (file): The LaTeX ouput file of the program.
            directory (str): Path to the result folder.
            write_time (bool): (Here to emulate polymorphism).

        """
        latex_file.write('\n\\begin{center}\n')
        draw_diagram(directory, latex_file, self.tags[0], 'diag')
        latex_file.write('\n\\end{center}\n\n')
