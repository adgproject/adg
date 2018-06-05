"""Tests for the adg.diag module."""

import pytest
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
