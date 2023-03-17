"""Tests for the adg.bmbpt module."""

import numpy as np
import networkx as nx
import adg.bmbpt
import adg.diag


def test_attribute_expressions():
    """Test the attribute_expressions method."""
    # Test for a simple diagram associated to a tree TSD
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT", -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    # Normally attributed in produce_expressions
    diag.tsd_is_tree = True

    diag.attribute_expressions(tsd)

    feyn_ref = '\\lim\\limits_{\\tau \\to \\infty}\\frac{(-1)^2 }{2(2!)^2}' \
        + '\\sum_{k_i}O^{40}_{k_{1}k_{2}k_{3}k_{4}} \\Omega^{02}_{k_{1}k_{2}}'\
        + ' \\Omega^{02}_{k_{3}k_{4}} \\int_{0}^{\\tau}\\mathrm{d}\\tau_1' \
        + '\\mathrm{d}\\tau_2e^{-\\tau_1 \\epsilon^{}_{k_{1}k_{2}}}' \
        + 'e^{-\\tau_2 \\epsilon^{}_{k_{3}k_{4}}}\n'

    diag_ref = '\\frac{(-1)^2 }{2(2!)^2}\\sum_{k_i}' \
        + '\\frac{O^{40}_{k_{1}k_{2}k_{3}k_{4}} \\Omega^{02}_{k_{1}k_{2}} ' \
        + '\\Omega^{02}_{k_{3}k_{4}} }{\\epsilon^{}_{k_{1}k_{2}}\\ ' \
        + '\\epsilon^{}_{k_{3}k_{4}}\\ } \n'

    assert diag.feynman_exp == feyn_ref
    assert diag.diag_exp == diag_ref

    # Test for a diagram associated to a non-tree TSD
    diagram = np.array([[0, 2, 2, 0],
                        [0, 0, 0, 2],
                        [0, 0, 0, 2],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT", -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)
    tsd.equivalent_trees = tsd.treat_cycles()

    diag.attribute_expressions(tsd)

    feyn_ref = '\\lim\\limits_{\\tau \\to \\infty}\\frac{(-1)^3 }{2(2!)^4}' \
        + '\\sum_{k_i}O^{40}_{k_{1}k_{2}k_{3}k_{4}} ' \
        + '\\Omega^{22}_{k_{5}k_{6}k_{1}k_{2}} ' \
        + '\\Omega^{22}_{k_{7}k_{8}k_{3}k_{4}} ' \
        + '\\Omega^{04}_{k_{7}k_{8}k_{5}k_{6}} \\int_{0}^{\\tau}' \
        + '\\mathrm{d}\\tau_1\\mathrm{d}\\tau_2\\mathrm{d}\\tau_3' \
        + '\\theta(\\tau_3-\\tau_1) \\theta(\\tau_3-\\tau_2) ' \
        + 'e^{-\\tau_1 \\epsilon^{k_{5}k_{6}}_{k_{1}k_{2}}}' \
        + 'e^{-\\tau_2 \\epsilon^{k_{7}k_{8}}_{k_{3}k_{4}}}' \
        + 'e^{-\\tau_3 \\epsilon^{}_{k_{5}k_{6}k_{7}k_{8}}}\n'

    diag_ref = '\\frac{(-1)^3 }{2(2!)^4}\\sum_{k_i}' \
        + 'O^{40}_{k_{1}k_{2}k_{3}k_{4}} \\Omega^{22}_{k_{5}k_{6}k_{1}k_{2}} '\
        + '\\Omega^{22}_{k_{7}k_{8}k_{3}k_{4}} ' \
        + '\\Omega^{04}_{k_{7}k_{8}k_{5}k_{6}} \\left[' \
        + '\\frac{1}{\\epsilon^{}_{k_{1}k_{2}k_{7}k_{8}}\\ ' \
        + '\\epsilon^{}_{k_{1}k_{2}k_{3}k_{4}}\\ ' \
        + '\\epsilon^{}_{k_{5}k_{6}k_{7}k_{8}}\\ } + ' \
        + '\\frac{1}{\\epsilon^{}_{k_{1}k_{2}k_{3}k_{4}}\\ ' \
        + '\\epsilon^{}_{k_{3}k_{4}k_{5}k_{6}}\\ ' \
        + '\\epsilon^{}_{k_{5}k_{6}k_{7}k_{8}}\\ } \\right]\n'

    assert diag.feynman_exp == feyn_ref
    assert diag.diag_exp == diag_ref


def test_vertex_expression():
    """Test the vertex_expression method."""
    # Set up the test case
    diagram = np.array([[0, 2, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    # attribute_qp_labels is used and thus tested in the process
    diag.attribute_qp_labels()

    assert diag.vertex_expression(0) == '\\epsilon^{k_{1}k_{2}k_{3}k_{4}}_{}'
    assert diag.vertex_expression(1) == '\\epsilon^{k_{5}k_{6}}_{k_{1}k_{2}}'
    assert diag.vertex_expression(2) == '\\epsilon^{}_{k_{3}k_{4}k_{5}k_{6}}'


def test_vertex_exchange_sym_factor():
    """Test the vertex_exchange_sym_factor method."""
    # Test case with symmetry factor
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.vertex_exchange_sym_factor == 2

    # Test case with higher symmetry factor
    diagram = np.array([[0, 2, 2, 2],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.vertex_exchange_sym_factor == 6

    # Test cases without symmetry factor
    diagram = np.array([[0, 2, 2], [0, 0, 2], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.vertex_exchange_sym_factor == 1

    diagram = np.array([[0, 0, 1, 3], [0, 0, 3, 1],
                        [0, 0, 0, 0], [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], 'BMBPT', -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.vertex_exchange_sym_factor == 1


def test_extract_integral():
    """Test for the extract_integral method."""
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT", -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()

    # Vertex expressions are needed prior to integral extraction
    diag.vert_exp = [diag.vertex_expression(vertex) for vertex in diag.graph]

    integral_ref = '\\mathrm{d}\\tau_1\\mathrm{d}\\tau_2' \
        + 'e^{-\\tau_1 \\epsilon^{}_{k_{1}k_{2}}}' \
        + 'e^{-\\tau_2 \\epsilon^{}_{k_{3}k_{4}}}'
    assert diag.extract_integral() == integral_ref


def test_extract_numerator():
    """Test for the extract_numerator method."""
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT", -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()

    numerator_ref = 'O^{40}_{k_{1}k_{2}k_{3}k_{4}} \\Omega^{02}_{k_{1}k_{2}}' \
        + ' \\Omega^{02}_{k_{3}k_{4}} '
    assert diag.extract_numerator() == numerator_ref


def test_has_crossing_sign():
    """Test for the has_crossing_sign method."""
    # Test case without crossing sign
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert not diag.has_crossing_sign()

    # Test case with a crossing sign
    diagram = np.array([[0, 2, 1, 1],
                        [0, 0, 1, 1],
                        [0, 0, 0, 1],
                        [0, 0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.has_crossing_sign()


def test_multiplicity_symmetry_factor():
    """Test the multiplicity_symmetry_factor method."""
    # Test a simple case with twice two lines
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.multiplicity_symmetry_factor() == '(2!)^2'

    # Test a case with once three lines
    diagram = np.array([[0, 3, 1], [0, 0, 1], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)

    assert diag.multiplicity_symmetry_factor() == '(3!)'


def test_time_tree_denominator():
    """Test for the time_tree_denominator method."""
    diagram = np.array([[0, 2, 2], [0, 0, 0], [0, 0, 0]])
    graph = nx.from_numpy_array(diagram,
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT", -1)
    diag = adg.bmbpt.BmbptFeynmanDiagram(graph, 0)
    diag.attribute_qp_labels()
    tsd = adg.tsd.TimeStructureDiagram(diag)

    denom_ref = "\\epsilon^{}_{k_{1}k_{2}}\\ \\epsilon^{}_{k_{3}k_{4}}\\ "
    assert diag.time_tree_denominator(tsd.graph) == denom_ref
