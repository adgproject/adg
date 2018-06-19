"""Tests for the adg.diag module."""

import pytest
import networkx as nx
import adg.diag


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
