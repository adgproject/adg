"""Tests for the adg.diag module."""

import os
import pytest
import networkx as nx
import numpy as np
import adg.diag


def test_feynmf_generator():
    """Test the feynmf_generator routine."""
    graph = nx.from_numpy_matrix(np.array([[0, 6], [0, 0]]),
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], "BMBPT")
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
        + "\\fmf{prop_pm,left=0.9}{v0,v1}\n" \
        + "\\fmf{prop_pm,right=0.9}{v0,v1}\n" \
        + "\\fmf{prop_pm,left=0.75}{v0,v1}\n" \
        + "\\fmf{prop_pm,right=0.75}{v0,v1}\n" \
        + "\\fmf{prop_pm,left=0.5}{v0,v1}\n" \
        + "\\fmf{prop_pm,right=0.5}{v0,v1}\n" \
        + "\\end{fmfgraph*}\n" \
        + "\\end{fmffile}}\n"

    with open("diag0.tex", "r") as result_file:
        assert result_file.read() == feynmf_ref
    os.unlink("diag0.tex")

    graph = nx.from_numpy_matrix(np.array([[0, 1], [1, 0]]),
                                 create_using=nx.MultiDiGraph(),
                                 parallel_edges=True)
    adg.diag.label_vertices([graph], "MBPT")
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
        + "\\fmf{half_prop,right=0.5}{v0,v1}\n" \
        + "\\fmf{half_prop,right=0.5}{v1,v0}\n" \
        + "\\end{fmfgraph*}\n" \
        + "\\end{fmffile}}\n"

    with open("diag1.tex", "r") as result_file:
        assert result_file.read() == feynmf_ref
        os.unlink("diag1.tex")


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

    # Test an unimplemented case
    with pytest.raises(KeyError):
        adg.diag.propagator_style('prop_pp')


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
