"""Routines and class for Bogoliubov IMSRG diagrams."""

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
import math
import networkx as nx
import numpy as np
import adg.diag
from adg.tools import reversed_enumerate


def two_partitions(number):
    """Return 2-partitions of the given integer.

    Args:
        numbr (int): The integer to partition.

    Returns:
        (list): All the 2-partitions as tuples.
    """
    return [(elem, number - elem) for elem in range(number+1)]


def diagrams_generation(orders):
    """Generate diagrams for B-IMSRG.

    Args:
        orders (tuple): The B-IMSRG (N_A, N_B, N_C) order of the diagrams.

    Returns:
        (tuple): NumPy arrays encoding the adjacency matrices of the graphs,
            and number of diagrams with the A vertex on top.

    """
    deg_ext_max = 2*orders[-1]
    diagrams_topA = diagrams_subset(2*orders[0], 2*orders[1], deg_ext_max)
    diagrams_topB = diagrams_subset(2*orders[1], 2*orders[0], deg_ext_max)

    return diagrams_topA + diagrams_topB, len(diagrams_topA)


def diagrams_subset(deg_max_top, deg_max_bot, deg_max_ext):
    """Generate diagrams for B-IMSRG.

    Args:
        orders (tuple): The B-IMSRG (N_A, N_B, N_C) order of the diagrams.

    Returns:
        (list): NumPy arrays encoding the adjacency matrices of the graphs.

    """
    matrices = []
    mat = np.zeros((4, 4), dtype=int)

    # Pick a valid vertex degree or zero for the external lines
    for deg_ext in range(0, deg_max_ext + 1, 2):
        # Split the valid vertex degree between external vertices
        for deg_0, deg_3 in two_partitions(deg_ext):
            # Split the vertex degree between the two internal vertices
            for part_0 in two_partitions(deg_0):
                mat[0,1], mat[0,2] = part_0
                for part_3 in two_partitions(deg_3):
                    mat[1,3], mat[2,3] = part_3
                    temp_deg_1 = mat[0,1] + mat[1,3]
                    temp_deg_2 = mat[0,2] + mat[2,3]
                    # Check that internal vertices are both odd/even, as lines
                    # between them will affect their degrees similarly
                    if ((abs(temp_deg_1 - temp_deg_2) % 2) == 0):
                        # Determine how many lines can connect the vertices
                        max_addition = min(deg_max_bot - temp_deg_1,
                                           deg_max_top - temp_deg_2)
                        # Check the odd/even character of vertices, since they
                        # must eventually be even
                        min_addition = 2 if (temp_deg_1 % 2 == 0) else 1
                        for addition in range(min_addition, max_addition+1, 2):
                            temp_mat = mat.copy()
                            temp_mat[1,2] = addition
                            matrices.append(temp_mat)
    return matrices


def order_diagrams(diagrams, order):
    """Order the BIMSRG diagrams and return the number of diags for each type.

    Args:
        diagrams (list): The unordered BimsrgDiagrams.
        order (int): The order of the B-IMSRG truncation.

    Returns:
        (tuple): First element are the ordered BimsrgDiagrams. Second element
        is the number of diagrams for each type. Third element is flags for the
        output processing.

    """
    diags_per_order = {}
    for n in range(1, order + 1):
        diags_per_order[n] = []

    for i_diag, diag in reversed_enumerate(diagrams):
        diags_per_order[diag.max_degree/2].append(diag)
        del diagrams[i_diag]

    diags_nb_per_type = {}

    for n in range(1, order + 1):
        diags_nb_per_type[n] = len(diags_per_order[n])

    diagrams = []
    for n in range(1, order + 1):
        diagrams += sorted(diags_per_order[n],
                           key=lambda diag: diag.ext_io_degree)
    diags_nb_per_type['nb_diags'] = len(diagrams)

    section_flags = {1: 0}

    for ind, diagram in enumerate(diagrams):
        diagram.tags[0] = ind
        if ind == 0:
            section_flags['new_op_struct'] = [ind]
        elif diagram.ext_io_degree != diagrams[ind-1].ext_io_degree:
            section_flags['new_op_struct'].append(ind)

    for n in range(2, order + 1):
        index = sum(len(diags_per_order[i]) for i in range(1, n))
        section_flags[n] = diagrams[index].tags[0] if diags_per_order[n] else -1

    return diagrams, diags_nb_per_type, section_flags


def write_header(tex_file, commands, diags_nbs):
    """Write overall header for BIMSRG result file.

    Args:
        tex_file (file): The ouput LaTeX file of the program.
        commands (Namespace): Flags for the program run.
        diags_nbs (dict): The number of diagrams per type.

    """
    tex_file.write("Valid diagrams: %i\n\n" % diags_nbs['nb_diags'])

    for n in range(1, commands.order[-1] + 1):
        tex_file.write("B-IMSRG(%i) diagrams: %i\n\n" % (n, diags_nbs[n]))


class BimsrgDiagram(adg.diag.Diagram):
    """Describes a B-IMSRG Feynman diagram with its related properties.

    Attributes:
        adjacency_mat (Numpy array): The adjacency matrix of the diagram.
        unique_id (int): A unique number associated to the diagram.
        expr (str): The B-IMSRG expression associated to the diagram.
        ext_io_degree (tuple): The degree of the operator component the diagram
            corresponds to.

    """

    __slots__ = ('adjacency_mat', 'unique_id', 'ext_io_degree',
                 '_vert_exchange_sym_fact', 'expr')

    def __init__(self, nx_graph, tag_num):
        """Generate a BMBPT diagrams using a NetworkX graph.

        Args:
            nx_graph (NetworkX MultiDiGraph): The graph of interest.
            tag_num (int): The tag number associated to the diagram.

        """
        adg.diag.Diagram.__init__(self, nx_graph)
        self.tags = [tag_num]
        self.unique_id = tag_num
        self.adjacency_mat = nx.to_numpy_array(self.graph, dtype=int)
        self.expr = self.attribute_expression()
        self.max_degree = max(self.unsort_degrees[0] + self.unsort_degrees[3],
                              self.unsort_degrees[1],
                              self.unsort_degrees[2])
        self.ext_io_degree = (self.unsort_degrees[0], self.unsort_degrees[3])

    def attribute_expression(self):
        """Returns the LaTeX expression of the diagram.

        Returns:
            (str): The LaTeX expression for the diagram.

        """
        A_degrees = self.unsort_io_degrees[1] \
            if self.graph.nodes[1]['operator'] == 'A' else self.unsort_io_degrees[2]
        B_degrees = self.unsort_io_degrees[1] \
            if self.graph.nodes[1]['operator'] == 'B' else self.unsort_io_degrees[2]
        name = " C^{%i%i}(%i%i,%i%i) = " \
            % (self.unsort_degrees[3], self.unsort_degrees[0],
               A_degrees[1], A_degrees[0], B_degrees[1], B_degrees[0])
        return name + self.permutator() + self.symmetry_factor() \
            + self.vertices_expression()

    def permutator(self):
        """Return the permutator associated to the diagram.

        Returns:
            (str): The permutator associated to the diagram in LaTeX format.

        """
        # Number of external lines tied to each vertex
        nb_out_2 = self.adjacency_mat[2, 3]
        nb_in_2 = self.adjacency_mat[0, 2]
        nb_out_1 = self.adjacency_mat[1, 3]
        nb_in_1 = self.adjacency_mat[0, 1]

        perm = ""

        # Permutator for out-going lines
        if (nb_out_2 != 0) and (nb_out_1 != 0):
            perm += "P(%s/%s) " \
                % ("".join('k_{%i}' % (idx + 1) for idx in range(nb_out_2)),
                   "".join('k_{%i}' % (idx + 1) for idx
                           in range(nb_out_2, nb_out_2 + nb_out_1)))

        # Permutator for incoming lines
        if (nb_in_2 != 0) and (nb_in_1 != 0):
            perm += "P(%s/%s) " \
                % ("".join('k_{%i}' % (idx + 1) for idx
                           in range(nb_out_2 + nb_out_1,
                                    nb_out_2 + nb_out_1 + nb_in_2)),
                   "".join('k_{%i}' % (idx + 1) for idx
                           in range(nb_out_2 + nb_out_1 + nb_in_2,
                                    nb_out_2 + nb_out_1 + nb_in_2 + nb_in_1)))
        return perm

    def symmetry_factor(self):
        """Returns the symmetry factor of the diagram in LaTeX format.

        Returns:
            (str): The LaTeX-formatted symmetry factor.

        """
        # Count the number of internal lines
        factor = math.factorial(self.adjacency_mat[1, 2])
        return "\\frac{1}{%i}" % factor if factor != 1 else ""

    def vertices_expression(self):
        """Return the expression associated to the vertices in LaTeX format.

        Returns:
            (str): The LaTeX-formatted expression for the vertices.
        """
        internal_lines = "".join("p_{%i}" % label for label
                                 in range(1, self.adjacency_mat[1, 2] + 1))
        # Number of external lines tied to each vertex
        nb_out_2 = self.adjacency_mat[2, 3]
        nb_in_2 = self.adjacency_mat[0, 2]
        nb_out_1 = self.adjacency_mat[1, 3]
        nb_in_1 = self.adjacency_mat[0, 1]
        # Expression associated to the upper vertex
        expr_2 = "%s^{%i%i}_{%s %s}" \
            % (self.graph.nodes[2]['operator'],
               self.unsort_io_degrees[2][1], self.unsort_io_degrees[2][0],
               "".join('k_{%i}' % (idx + 1) for idx in range(nb_out_2))
               + "".join('k_{%i}' % (idx + 1) for idx
                         in range(nb_out_2 + nb_out_1,
                                  nb_out_2 + nb_out_1 + nb_in_2)),
               internal_lines)
        # Expression associated to the lower vertex
        expr_1 = "%s^{%i%i}_{%s %s}" \
            % (self.graph.nodes[1]['operator'],
               self.unsort_io_degrees[1][1], self.unsort_io_degrees[1][0],
               internal_lines,
               "".join('k_{%i}' % (idx + 1) for idx
                       in range(nb_out_2, nb_out_2 + nb_out_1))
               + "".join('k_{%i}' % (idx + 1) for idx
                         in range(nb_out_2 + nb_out_1 + nb_in_2,
                                  nb_out_2 + nb_out_1 + nb_in_2 + nb_in_1)))

        return "\\sum_{%s} %s %s" % (internal_lines, expr_2, expr_1)

    def write_graph(self, latex_file, directory, write_time):
        """Write the BMBPT graph and its associated TSD to the LaTeX file.

        Args:
            latex_file (file): The LaTeX output file of the program.
            directory (str): The path to the result folder.
            write_time (bool): ``True`` if we want informations on the
                associated TSDs.

        """
        latex_file.write('\n\\begin{center}\n')
        adg.diag.draw_diagram(directory, latex_file, str(self.tags[0]), 'diag')
        latex_file.write('\n\\end{center}\n\n')

    def write_section(self, result, commands, section_flags):
        """Write section and subsections for BMBPT result file.

        Args:
            result (file): The LaTeX output file of the program.
            commands (dict): The flags associated with run management.
            section_flags (dict): UniqueIDs of diags starting each section.

        """
        for n in range(1, commands.order[-1] + 1):
            if self.tags[0] == section_flags[n]:
                result.write("\\section{B-IMSRG(%i)}\n\n" % n)
        if self.tags[0] in section_flags['new_op_struct']:
            result.write("\\subsection{$C^{%i%i}$}\n\n"
                         % (self.ext_io_degree[1], self.ext_io_degree[0]))
        result.write("\\paragraph{Diagram %i:}\n" % (self.tags[0] + 1))
        result.write("\\begin{equation}\\n%s\\n\\end{equation}\\n" % self.expr)
