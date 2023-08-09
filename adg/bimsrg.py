"""Routines and class for Bogoliubov IMSRG diagrams."""

from future import standard_library
standard_library.install_aliases()
from builtins import range
import math
import networkx as nx
import numpy as np
from itertools import permutations, product
from more_itertools import powerset
import adg.diag
from adg.tools import reversed_enumerate


def two_partitions(number):
    """Return 2-partitions of the given integer.

    Args:
        numbr (int): The integer to partition.

    Returns:
        list: All the 2-partitions as tuples.

    >>> two_partitions(3)
    [(0, 3), (1, 2), (2, 1), (3, 0)]

    >>> two_partitions(0)
    [(0, 0)]

    >>> two_partitions(-1)
    []

    """
    return [(elem, number - elem) for elem in range(number+1)]


def diagrams_generation(orders):
    """Generate diagrams for B-IMSRG.

    Args:
        orders (tuple): The B-IMSRG (N_A, N_B, N_C) order of the diagrams.

    Returns:
        tuple: NumPy arrays encoding the adjacency matrices of the graphs,
            and number of diagrams with the A vertex on top.

    """
    deg_ext_max = 2*orders[-1]
    diagrams_topA = diagrams_subset(2*orders[0], 2*orders[1], deg_ext_max)
    diagrams_topB = diagrams_subset(2*orders[1], 2*orders[0], deg_ext_max)

    return diagrams_topA + diagrams_topB, len(diagrams_topA)


def diagrams_subset(deg_max_top, deg_max_bot, deg_max_ext):
    """Generate diagrams for B-IMSRG.

    Args:
        orders (tuple): The max ranks (2*N_A, 2*N_B, 2*N_C) of the vertices.

    Returns:
        list: NumPy arrays encoding the adjacency matrices of the graphs.

    >>> diagrams_subset(2, 2, 0) # doctest: +NORMALIZE_WHITESPACE
    [array([[0, 0, 0, 0],
           [0, 0, 2, 0],
           [0, 0, 0, 0],
           [0, 0, 0, 0]])]
    >>> len(diagrams_subset(2, 2, 2))
    5
    >>> len(diagrams_subset(4, 4, 4))
    41
    >>> len(diagrams_subset(0, 0, 0))
    0

    """
    matrices = []
    mat = np.zeros((4, 4), dtype=int)

    # Pick a valid vertex degree or zero for the external lines
    for deg_ext in range(0, deg_max_ext + 1, 2):
        # Split the valid vertex degree between external vertices
        for deg_0, deg_3 in two_partitions(deg_ext):
            # Split the vertex degree between the two internal vertices
            for part_0, part_3 in product(two_partitions(deg_0),
                                          two_partitions(deg_3)):
                mat[0,1], mat[0,2] = part_0
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
        tuple: First element are the ordered BimsrgDiagrams. Second element
        is the number of diagrams for each type. Third element is flags for the
        output processing.

    """
    diags_per_order = {n: [] for n in range(1, order + 1)}

    for i_diag, diag in reversed_enumerate(diagrams):
        diags_per_order[diag.max_degree/2].append(diag)
        del diagrams[i_diag]

    diags_nb_per_type = {n: len(diags_per_order[n])
                         for n in range(1, order + 1)}

    diagrams = []
    for n in range(1, order + 1):
        diagrams += sorted(diags_per_order[n],
                           key=lambda diag: (diag.ext_io_degree,
                                             not diag.is_AB))
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
    tex_file.write("$C=\\left[A,B\\right]$ "
                   + "with $N_A = %i$, $N_B = %i$ and $N_C = %i$\n\n"
                   % commands.order)
    tex_file.write("$d_\\mathrm{max} \\equiv \\mathrm{max}(d_A, d_B, d_C)\n\n")
    tex_file.write("\\par \\vspace{\\baselineskip}\n")
    tex_file.write("Valid diagrams: %i\n\n" % diags_nbs['nb_diags'])

    for n in range(1, commands.order[-1] + 1):
        tex_file.write("$d_\\mathrm{max} = %i$ diagrams: %i\n\n"
                       % (n, diags_nbs[n]))


def write_permutator_section(tex_file, commands):
    """Write the section defining permutation operators.

    Args:
        tex_file (file): The ouput LaTeX file of the program.
        commands (Namespace): Flags for the program run.

    """
    tex_file.write("\\section{Permutators definitions}\n\n"
                   + "\\begin{align*}\n")
    # Number of external legs is at most two times N_C
    for max_index in range(2, 2*commands.order[-1] + 1):
        indices = list(range(1, max_index + 1))
        # Split our qp states in two subgroups of varying length
        for partition in ((indices[:n], indices[n:])
                          for n in range(1, len(indices)//2 + 1)):
            tex_file.write(permutator(partition[0], partition[1]))
    tex_file.write("\\end{align*}\n")


def permutator(set_1, set_2):
    """Write the definition of the permutation operator given by the two sets.

    Args:
        set_1 (list): The list of the left-hand-side qp labels.
        set_2 (list): The list of the right-hand-side qp labels.

    Returns:
        str: The LaTeX expression for the permutation operator.

    >>> print(permutator([1, 2], [3])) # doctest: +NORMALIZE_WHITESPACE
    P(k_{1}k_{2}/k_{3}) &= 1 - P_{k_{1} k_{3}} - P_{k_{2} k_{3}} \\\\\n

    """
    perm_op = 'P(%s/%s) &= 1 ' % ("".join('k_{%i}' % label for label in set_1),
                                  "".join('k_{%i}' % label for label in set_2))
    perm_counter = 0
    # Combine all possible subsets of both sets
    for subset_1, subset_2 in product(powerset(set_1), powerset(set_2)):
        # Ensure each qp state has a permutation partner, exclude the empty set
        if (len(subset_1) == len(subset_2)) and (len(subset_1) > 0):
            # Run through all permutations of the first subset
            for permutation in permutations(subset_1):
                perm_op += '+ ' if len(subset_1) % 2 == 0 else '- '
                for elem in zip(permutation, subset_2):
                    perm_op += 'P_{k_{%i} k_{%i}} ' % elem
                    perm_counter += 1
                # Avoid permutators lines being too long
                if perm_counter >= 10:
                    perm_op += '\\\\\n &\\phantom{=} '
                    perm_counter = 0
    if perm_op.endswith('\\\\\n &\\phantom{=} '):
        perm_op = perm_op[:-len(' &\\phantom{=} ')]
    else:
        perm_op += "\\\\\n"
    return perm_op


class BimsrgDiagram(adg.diag.Diagram):
    """Describes a B-IMSRG Feynman diagram with its related properties.

    Args:
        nx_graph (NetworkX MultiDiGraph): The graph of interest.
        tag_num (int): The tag number associated to the diagram.

    """

    __slots__ = ('adjacency_mat', 'unique_id', 'ext_io_degree',
                 '_vert_exchange_sym_fact', 'expr', 'is_AB')

    def __init__(self, nx_graph, tag_num):
        """Generate a BMBPT diagrams using a NetworkX graph."""
        adg.diag.Diagram.__init__(self, nx_graph)
        self.tags = [tag_num]
        self.unique_id = tag_num
        """int: A unique number associated to the diagram."""
        self.adjacency_mat = nx.to_numpy_array(self.graph, dtype=int)
        """Numpy array: The adjacency matrix of the diagram."""
        self.expr = self.attribute_expression()
        """str: The B-IMSRG expression associated to the diagram."""
        self.max_degree = max(self.unsort_degrees[0] + self.unsort_degrees[3],
                              self.unsort_degrees[1],
                              self.unsort_degrees[2])
        """int: Maximal degree over vertices or external lines (A, B, C)."""
        self.ext_io_degree = (self.unsort_degrees[0], self.unsort_degrees[3])
        """tuple: The degree of the operator component the diagram
        corresponds to."""
        self.is_AB = True if self.graph.nodes[2]['operator'] == 'A' else False
        """bool: True if the diagram contributes to +AB, false if to -BA."""

    def attribute_expression(self):
        """Returns the LaTeX expression of the diagram.

        The expression is extracted in a way that assumes the canonical
        representation of the vertices as well as a labeling of the external
        lines that would correspond to a canonical representation of the C
        operator vertex. As such, there is no sign factor associated to
        departing from the canonical representation of the diagram. This
        additionally prevents any line crossing from appearing.

        Returns:
            str: The LaTeX expression for the diagram.

        """
        A_degrees = self.unsort_io_degrees[1] \
            if self.graph.nodes[1]['operator'] == 'A' \
            else self.unsort_io_degrees[2]
        B_degrees = self.unsort_io_degrees[1] \
            if self.graph.nodes[1]['operator'] == 'B' \
            else self.unsort_io_degrees[2]
        name = " C^{%i%i}_{%s}(%i%i,%i%i) = " \
            % (self.unsort_degrees[3], self.unsort_degrees[0],
               "".join('k_{%i}' % (idx + 1) for idx
                       in range(self.unsort_degrees[3]
                                + self.unsort_degrees[0])),
               A_degrees[1], A_degrees[0], B_degrees[1], B_degrees[0])
        return name + self.sign() + self.permutator() \
            + self.symmetry_factor() + self.vertices_expression()

    def sign(self):
        """Return the sign of the diagram.

        As the diagrams are made with vertices in canonical representation and
        avoiding crossings, the only sign left is associated with the
        commutator.

        Returns:
            str: The sign of the diagram.

        """
        return '- ' if self.graph.nodes[2]['operator'] == 'B' else ''

    def permutator(self):
        """Return the permutator associated to the diagram.

        The labelling of the external lines correspond to the canonical
        representation of the C operator vertex.

        Returns:
            str: The permutator associated to the diagram in LaTeX format.

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
            str: The LaTeX-formatted symmetry factor.

        """
        # Count the number of internal lines
        factor = math.factorial(self.adjacency_mat[1, 2])
        return "\\frac{1}{%i}" % factor if factor != 1 else ""

    def vertices_expression(self):
        """Return the expression associated to the vertices in LaTeX format.

        The expressions are extracted in a way that assumes the canonical
        representation of the vertices as well as a labeling of the external
        lines that would correspond to a canonical representation of the C
        operator vertex. This prevents any crossing as well as additional sign
        tied to departing from the canonical representation.

        Returns:
            str: The LaTeX-formatted expression for the vertices.
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

    def write_section(self, result, commands, section_flags):
        """Write section and subsections for BMBPT result file.

        Args:
            result (file): The LaTeX output file of the program.
            commands (dict): The flags associated with run management.
            section_flags (dict): UniqueIDs of diags starting each section.

        """
        for n in range(1, commands.order[-1] + 1):
            if self.tags[0] == section_flags[n]:
                result.write("\\section{$d_\\mathrm{max} = %i$}\n\n" % n)
        if self.tags[0] in section_flags['new_op_struct']:
            result.write("\\subsection{$C^{%i%i}$}\n\n"
                         % (self.ext_io_degree[1], self.ext_io_degree[0]))
        result.write("\\paragraph{Diagram %i (%s):}\n"
                     % (self.tags[0] + 1, '$+AB$' if self.is_AB else '$-BA$'))
        result.write("\\begin{equation}\n%s\n\\end{equation}\n" % self.expr)
