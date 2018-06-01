"""Some tests for adg.run."""

import argparse
import pytest
import adg.run


def test_generate_diagrams():
    """Unit and regression test for generate_diagrams."""
    com = argparse.Namespace()

    # Tests for the number of diagrams produced for simple known cases

    com.theory, com.order, com.with_three_body = 'BMBPT', 2, False
    assert len(adg.run.generate_diagrams(com)) == 2

    com.theory, com.order, com.with_three_body = 'MBPT', 2, False
    assert len(adg.run.generate_diagrams(com)) == 1

    com.theory, com.order, com.with_three_body = 'BMBPT', 2, True
    assert len(adg.run.generate_diagrams(com)) == 3

    # Test for anomalous cases

    with pytest.raises(TypeError):
        adg.run.generate_diagrams()

    com.theory, com.order, com.with_three_body = 'SCGF', 2, True
    with pytest.raises(UnboundLocalError):
        adg.run.generate_diagrams(com)


def test_order_diagrams():
    """Unit and regression test for order_diagrams."""
    com = argparse.Namespace()

    # Tests for the number of diagrams produced for simple known cases

    com.theory, com.order, com.with_three_body = 'BMBPT', 2, False

    # Use generate_diagrams as a seed
    diagrams = adg.run.generate_diagrams(com)
    assert len(diagrams) == 2

    # Test that diagrams are processed
    assert len(adg.run.order_diagrams(diagrams, com)[0]) == 2

    # Test for the removing of topologically identical diagrams
    diagrams = adg.run.generate_diagrams(com) + adg.run.generate_diagrams(com)
    assert len(diagrams) == 4
    assert len(adg.run.order_diagrams(diagrams, com)[0]) == 2

    # Test for the ordering with simple three-body case
    com.theory, com.order, com.with_three_body = 'BMBPT', 2, True
    diagrams = adg.run.generate_diagrams(com)
    assert len(diagrams) == 3
    diagrams, diag_nbs = adg.run.order_diagrams(diagrams, com)
    assert len(diagrams) == 3
    assert diagrams[-1].two_or_three_body == 3
    assert diag_nbs['nb_2'] == 2
    assert diag_nbs['nb_3'] == 1
    assert diag_nbs['nb_2_hf'] == 1
    assert diag_nbs['nb_2_ehf'] == 0
    assert diag_nbs['nb_2_not_hf'] == 1

    # Test that an empty list does not break the code
    assert not adg.run.order_diagrams([], com)[0]
