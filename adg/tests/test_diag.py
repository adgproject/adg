"""Tests for the adg.diag module."""

import os
import pytest
import networkx as nx
import numpy as np
import adg.diag


def test_feynmf_generator():
    """Test the feynmf_generator routine."""
    graph = nx.from_numpy_array(np.array([[0, 6], [0, 0]]),
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT", -1)
    adg.diag.feynmf_generator(graph, "BMBPT", "diag0")
    feynmf_ref = "\\parbox{40pt}{\\begin{fmffile}{diag0}\n" \
        + "\\begin{fmfgraph*}(40,40)\n" \
        + "\\fmfcmd{style_def prop_pm expr p =\n" \
        + "draw_plain p;\n" \
        + "shrink(.7);\n" \
        + "\tcfill (marrow (p, .25));\n" \
        + "\tcfill (marrow (p, .75))\n" \
        + "endshrink;\n" \
        + "enddef;}\n" \
        + "\\fmftop{v1}\\fmfbottom{v0}\n" \
        + "\\fmf{phantom}{v0,v1}\n" \
        + "\\fmfv{d.shape=square,d.filled=full,d.size=3thick}{v0}\n" \
        + "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v1}\n" \
        + "\\fmffreeze\n" \
        + "\\fmf{prop_pm,right=0.9}{v0,v1}\n" \
        + "\\fmf{prop_pm,right=0.75}{v0,v1}\n" \
        + "\\fmf{prop_pm,right=0.5}{v0,v1}\n" \
        + "\\fmf{prop_pm,left=0.5}{v0,v1}\n" \
        + "\\fmf{prop_pm,left=0.75}{v0,v1}\n" \
        + "\\fmf{prop_pm,left=0.9}{v0,v1}\n" \
        + "\\end{fmfgraph*}\n" \
        + "\\end{fmffile}}\n"

    with open("diag0.tex", "r") as result_file:
        assert result_file.read() == feynmf_ref
    os.unlink("diag0.tex")

    graph = nx.from_numpy_array(np.array([[0, 1], [1, 0]]),
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "MBPT", -1)
    adg.diag.feynmf_generator(graph, "MBPT", "diag1")
    feynmf_ref = "\\parbox{40pt}{\\begin{fmffile}{diag1}\n" \
        + "\\begin{fmfgraph*}(40,40)\n" \
        + "\\fmfcmd{style_def half_prop expr p =\n" \
        + "draw_plain p;\n" \
        + "shrink(.7);\n" \
        + "\tcfill (marrow (p, .5))\n" \
        + "endshrink;\n" \
        + "enddef;}\n" \
        + "\\fmftop{v1}\\fmfbottom{v0}\n" \
        + "\\fmf{phantom}{v0,v1}\n" \
        + "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v0}\n" \
        + "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v1}\n" \
        + "\\fmffreeze\n" \
        + "\\fmf{half_prop,right=0.5}{v1,v0}\n" \
        + "\\fmf{half_prop,right=0.5}{v0,v1}\n" \
        + "\\end{fmfgraph*}\n" \
        + "\\end{fmffile}}\n"

    with open("diag1.tex", "r") as result_file:
        assert result_file.read() == feynmf_ref
        os.unlink("diag1.tex")

    graph = nx.from_numpy_array(np.array([[0, 0, 0, 0],
                                          [0, 0, 1, 1],
                                          [0, 0, 0, 1],
                                          [0, 0, 0, 0]]),
                                create_using=nx.MultiDiGraph(),
                                parallel_edges=True)
    adg.diag.label_vertices([graph], "BIMSRG", 2)
    adg.diag.feynmf_generator(graph, "BIMSRG", "diag2")
    feynmf_ref = "\\parbox{40pt}{\\begin{fmffile}{diag2}\n" \
        + "\\begin{fmfgraph*}(40,80)\n" \
        + "\\fmfcmd{style_def half_prop expr p =\n" \
        + "draw_plain p;\n" \
        + "shrink(.7);\n" \
        + "\tcfill (marrow (p, .5))\n" \
        + "endshrink;\n" \
        + "enddef;}\n" \
        + "\\fmfstraight\n" \
        + "\\fmftopn{t}{3}\\fmfbottomn{b}{1}\n" \
        + "\\fmf{phantom}{b1,v1}\n" \
        + "\\fmf{phantom}{v1,v2}\n" \
        + "\\fmf{phantom}{v2,t2}\n" \
        + "\\fmfv{d.shape=circle,d.filled=empty,d.size=3thick}{v1}\n" \
        + "\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v2}\n" \
        + "\\fmffreeze\n" \
        + "\\fmf{half_prop}{v1,v2}\n" \
        + "\\fmf{half_prop,right=0.4}{v1,t3}\n" \
        + "\\fmf{half_prop,left=0.3}{v2,t1}\n" \
        + "\\end{fmfgraph*}\n" \
        + "\\end{fmffile}}\n"

    with open("diag2.tex", "r") as result_file:
        assert result_file.read() == feynmf_ref
        os.unlink("diag2.tex")


def test_propagator_style():
    """Test for expected and unexpected styles of propagator."""
    # Test for implemented cases
    assert adg.diag.propagator_style('prop_pm') == (
        "\\fmfcmd{style_def prop_pm expr p =\n"
        "draw_plain p;\n"
        "shrink(.7);\n"
        "\tcfill (marrow (p, .25));\n"
        "\tcfill (marrow (p, .75))\n"
        "endshrink;\n"
        "enddef;}\n"
    )

    assert adg.diag.propagator_style('half_prop') == (
        "\\fmfcmd{style_def half_prop expr p =\n"
        "draw_plain p;\n"
        "shrink(.7);\n"
        "\tcfill (marrow (p, .5))\n"
        "endshrink;\n"
        "enddef;}\n"
    )

    assert adg.diag.propagator_style('prop_mm') == (
        "\\fmfcmd{style_def prop_mm expr p =\n"
        "draw_plain p;\nshrink(.7);\n"
        "\tcfill (marrow (p, .75));\n"
        "\tcfill (marrow (reverse p, .75))\n"
        "endshrink;\nenddef;}\n"
    )

    # Test an unimplemented case
    with pytest.raises(KeyError):
        adg.diag.propagator_style('prop_pp')


def test_prop_directions():
    """Test the prop_directions function."""
    # Simplest case
    assert adg.diag.prop_directions(2, 5) == [',right=0.9', ',right=0.75',
                                              ',right=0.6', ',left=0.6',
                                              ',left=0.75', ',left=0.9']
    # Close vertices
    assert adg.diag.prop_directions(1, 5) == [',right=0.75', ',right=0.5', '',
                                              ',left=0.5', ',left=0.75',
                                              ',left=0.9']
    # Numerous props with close vertices
    assert adg.diag.prop_directions(1, 8) == [',right=0.70', ',right=0.60',
                                              ',right=0.50', ',right=0.40',
                                              ',right=0.30', ',right=0.20',
                                              ',right=0.10', ',left=0.10',
                                              ',left=0.20', ',left=0.30',
                                              ',left=0.40', ',left=0.50',
                                              ',left=0.60', ',left=0.70']


def test_to_skeleton():
    """Test the function returning the skeleton of a graph."""
    # Test for a graph with one redundant link
    graph = nx.MultiDiGraph()
    graph.add_edges_from([(0, 1), (1, 2), (0, 2)])
    assert list(adg.diag.to_skeleton(graph).edges()) == [(0, 1), (1, 2)]

    # Test for a graph without redundant link
    graph = nx.MultiDiGraph()
    graph.add_edges_from([(0, 1), (1, 2), (0, 3)])
    assert list(adg.diag.to_skeleton(graph).edges()) == [(0, 1), (0, 3), (1, 2)]

    # Test for a graph with twice the same edge
    graph = nx.MultiDiGraph()
    graph.add_edges_from([(0, 1), (1, 2), (1, 2)])
    assert list(adg.diag.to_skeleton(graph).edges()) == [(0, 1), (1, 2)]

    # Test for an empty graph
    graph = nx.MultiDiGraph()
    assert list(adg.diag.to_skeleton(graph).edges()) == []
