"""Routines and class for Bogoliubov IMSRG diagrams."""

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
import numpy as np
import adg.diag

def two_partitions(number):
    """Return 2-partitions of the given integer.

    Args:
        numbr (int): The integer to partition.

    Returns:
        (list): All the 2-partitions as tuples.
    """
    return [(elem, number - elem) for elem in range(number+1)]


def diagrams_generation(order):
    """Generate diagrams for B-IMSRG.

    Args:
        order (int): The B-IMSRG order of the studied diagrams.

    Returns:
        (list): NumPy arrays encoding the adjacency matrices of the graphs.

    """
    matrices = []
    deg_max = 2*order
    mat = np.zeros((4, 4), dtype=int)

    # Pick a valid vertex degree or zero for the external lines
    for deg_ext in range(deg_max + 1, 2):
        for deg_0 in range(deg_ext + 1):
            # Split the valid vertex degree between external vertices
            deg_3 = deg_ext - deg_0
            # Split the vertex degree between the two internal vertices
            for part_0 in two_partitions(deg_0):
                for part_3 in two_partitions(deg_3):
                    mat[0,1], mat[0,2] = part_0
                    mat[1,3], mat[2,3] = part_3
                    temp_deg_1 = mat[0,1] + mat[1,3]
                    temp_deg_2 = mat[0,2] + mat[2,3]
                    # Check that internal vertices are both odd / even, as lines
                    # between them will affect their degrees similarly
                    if ((abs(temp_deg_1 - temp_deg_2) % 2) != 0):
                        continue
                    # Determine how many lines can connect internal vertices
                    max_addition = min(deg_max - temp_deg_1, deg_max - temp_deg_2)
                    min_addition = 2 if (temp_deg_1 % 2 == 0) else 1
                    for addition in range(min_addition, max_addition+1, 2):
                        temp_mat = mat.copy()
                        temp_mat[1,2] = addition
                        matrices.append(temp_mat)
    return matrices


def write_header(tex_file, commands, diags_nbs):
    """Write overall header for BMBPT result file.

    Args:
        tex_file (file): The ouput LaTeX file of the program.
        commands (Namespace): Flags for the program run.
        diags_nbs (dict): The number of diagrams per type.

    """


class BimsrgDiagram(adg.diag.Diagram):
    """Describes a B-IMSRG Feynman diagram with its related properties.

    Attributes:
        two_or_three_body (int): The 2 or 3-body characted of the vertices.
        unique_id (int): A unique number associated to the diagram.

    """

    __slots__ = ('two_or_three_body', 'unique_id',
                 '_vert_exchange_sym_fact')

    def __init__(self, nx_graph, tag_num):
        """Generate a BMBPT diagrams using a NetworkX graph.

        Args:
            nx_graph (NetworkX MultiDiGraph): The graph of interest.
            tag_num (int): The tag number associated to the diagram.

        """
        adg.diag.Diagram.__init__(self, nx_graph)
        self.two_or_three_body = 3 if self.max_degree == 6 else 2
        self.tags = [tag_num]
        self.unique_id = tag_num

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
        result.write("\\paragraph{Diagram %i:}\n" % (self.tags[0] + 1))
