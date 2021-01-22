"""Some tests for adg.run."""

import argparse
import pytest
import adg.run
import adg.tools


def test_parse_command_line():
    """Test the parse_command_line function."""
    # Test with simple BMBPT case
    arguments = ['-t', 'BMBPT', '-o', '2']
    commands = adg.run.parse_command_line(arguments)
    assert commands.order == 2
    assert commands.theory == 'BMBPT'
    assert commands.nbody_observable == 2
    assert not commands.with_3NF
    assert not commands.interactive
    assert not commands.draw_diags
    assert not commands.canonical
    assert not commands.cd_output
    assert not commands.draw_tsds
    assert not commands.compile

    # Test with complex BMBPT case
    arguments = ['-t', 'BMBPT', '-o', '2', '-d', '--with_3NF', '-nobs', '3']
    commands = adg.run.parse_command_line(arguments)
    assert commands.order == 2
    assert commands.theory == 'BMBPT'
    assert commands.nbody_observable == 3
    assert commands.with_3NF
    assert not commands.interactive
    assert commands.draw_diags
    assert not commands.canonical
    assert not commands.cd_output
    assert not commands.draw_tsds
    assert not commands.compile

    # Test MBPT case
    arguments = ['-t', 'MBPT', '-o', '2', '-d', '-cd']
    commands = adg.run.parse_command_line(arguments)
    assert commands.order == 2
    assert commands.theory == 'MBPT'
    assert commands.nbody_observable == 2
    assert not commands.with_3NF
    assert not commands.interactive
    assert commands.draw_diags
    assert not commands.canonical
    assert commands.cd_output
    assert not commands.draw_tsds
    assert not commands.compile

    # Standard BIMSRG case
    arguments = ['-t', 'BIMSRG', '-o', '2', '-d']
    commands = adg.run.parse_command_line(arguments)
    assert commands.order == (2, 2, 2)
    assert commands.theory == 'BIMSRG'
    assert commands.nbody_observable == 2
    assert not commands.with_3NF
    assert not commands.interactive
    assert commands.draw_diags
    assert not commands.canonical
    assert not commands.cd_output
    assert not commands.draw_tsds
    assert not commands.compile

    # BIMSRG case with commutator truncation
    arguments = ['-t', 'BIMSRG', '-o', '3', '2', '-d']
    commands = adg.run.parse_command_line(arguments)
    assert commands.order == (3, 3, 2)
    assert commands.theory == 'BIMSRG'
    assert commands.nbody_observable == 2
    assert not commands.with_3NF
    assert not commands.interactive
    assert commands.draw_diags
    assert not commands.canonical
    assert not commands.cd_output
    assert not commands.draw_tsds
    assert not commands.compile

    # BIMSRG case with different truncation for each operator
    arguments = ['-t', 'BIMSRG', '-o', '3', '2', '1', '-d']
    commands = adg.run.parse_command_line(arguments)
    assert commands.order == (3, 2, 1)
    assert commands.theory == 'BIMSRG'
    assert commands.nbody_observable == 2
    assert not commands.with_3NF
    assert not commands.interactive
    assert commands.draw_diags
    assert not commands.canonical
    assert not commands.cd_output
    assert not commands.draw_tsds
    assert not commands.compile

    # Anomalous BIMSRG case with too many orders
    arguments = ['-t', 'BIMSRG', '-o', '3', '2', '1', '4']
    with pytest.raises(SystemExit):
        commands = adg.run.parse_command_line(arguments)

    # Anomalous non-BIMSRG case with too many orders
    arguments = ['-t', 'MBPT', '-o', '3', '2']
    with pytest.raises(SystemExit):
        commands = adg.run.parse_command_line(arguments)

    # Anomalous case with missing flag
    arguments = ['-t', 'BIMSRG']
    with pytest.raises(SystemExit):
        commands = adg.run.parse_command_line(arguments)


def test_generate_diagrams():
    """Unit and regression test for generate_diagrams."""
    com = argparse.Namespace()
    id_gen = adg.tools.UniqueID()

    # Tests for the number of diagrams produced for simple known cases

    com.theory, com.order, com.with_3NF = 'BMBPT', 1, False
    com.nbody_observable, com.canonical = 2, False
    assert len(adg.run.generate_diagrams(com, id_gen)) == 2

    com.theory, com.order, com.with_3NF = 'MBPT', 4, False
    com.nbody_observable, com.canonical = 2, False
    assert len(adg.run.generate_diagrams(com, id_gen)) == 39

    com.theory, com.order, com.with_3NF = 'BMBPT', 1, True
    com.nbody_observable, com.canonical = 3, False
    assert len(adg.run.generate_diagrams(com, id_gen)) == 3

    com.theory, com.order, com.with_3NF = 'PBMBPT', 1, False
    com.nbody_observable, com.canonical = 2, False
    assert len(adg.run.generate_diagrams(com, id_gen)) == 3

    com.theory, com.order, com.with_3NF = 'BIMSRG', (2, 2, 2), False
    assert len(adg.run.generate_diagrams(com, id_gen)) == 82

    # Test for anomalous cases

    with pytest.raises(TypeError):
        adg.run.generate_diagrams()

    com.theory, com.order, com.with_three_body = 'SCGF', 2, True
    com.nbody_observable, com.canonical = 2, False
    with pytest.raises(SystemExit):
        adg.run.generate_diagrams(com, id_gen)


def test_order_diagrams():
    """Unit and regression test for order_diagrams."""
    com = argparse.Namespace()
    id_gen = adg.tools.UniqueID()

    # Tests for the number of diagrams produced for simple known cases

    com.theory, com.order, com.with_3NF = 'BMBPT', 1, False
    com.nbody_observable, com.canonical = 2, False

    # Use generate_diagrams as a seed
    diagrams = adg.run.generate_diagrams(com, id_gen)
    assert len(diagrams) == 2

    # Test that diagrams are processed
    assert len(adg.run.order_diagrams(diagrams, com)[0]) == 2

    # Test for the ordering with simple three-body case
    com = argparse.Namespace()
    com.theory, com.order, com.with_3NF = 'BMBPT', 1, True
    com.nbody_observable, com.canonical = 3, False
    diagrams = adg.run.generate_diagrams(com, id_gen)
    assert len(diagrams) == 3
    diagrams, diag_nbs, _ = adg.run.order_diagrams(diagrams, com)
    assert len(diagrams) == 3
    assert diagrams[-1].two_or_three_body == 3
    assert diag_nbs['nb_2'] == 2
    assert diag_nbs['nb_3'] == 1
    assert diag_nbs['nb_2_hf'] == 1
    assert diag_nbs['nb_2_ehf'] == 0
    assert diag_nbs['nb_2_not_hf'] == 1

    # Similar test for PBMBPT
    com = argparse.Namespace()

    # Tests for the number of diagrams produced for simple known cases

    com.theory, com.order, com.with_3NF = 'PBMBPT', 1, False
    com.nbody_observable, com.canonical = 2, False

    # Use generate_diagrams as a seed
    diagrams = adg.run.generate_diagrams(com, id_gen)
    assert len(diagrams) == 3

    # Test that diagrams are processed
    assert len(adg.run.order_diagrams(diagrams, com)[0]) == 3

    # Test for the ordering with simple three-body case
    com = argparse.Namespace()
    com.theory, com.order, com.with_3NF = 'PBMBPT', 1, True
    com.nbody_observable, com.canonical = 3, False
    diagrams = adg.run.generate_diagrams(com, id_gen)
    assert len(diagrams) == 6
    diagrams, diag_nbs, _ = adg.run.order_diagrams(diagrams, com)
    assert len(diagrams) == 6
    assert diagrams[-1].two_or_three_body == 3
    assert diag_nbs['nb_2'] == 3
    assert diag_nbs['nb_3'] == 3
    assert diag_nbs['nb_2_hf'] == 1
    assert diag_nbs['nb_2_ehf'] == 1
    assert diag_nbs['nb_2_not_hf'] == 1

    # Test that an empty list does not break the code
    assert not adg.run.order_diagrams([], com)[0]

    # Simple tests calling for the ordering of MBPT diagrams
    com = argparse.Namespace()
    com.theory, com.order = 'MBPT', 3

    # Use generate_diagrams as a seed
    diagrams = adg.run.generate_diagrams(com, id_gen)
    assert len(diagrams) == 3

    # Test that diagrams are processed
    assert len(adg.run.order_diagrams(diagrams, com)[0]) == 3

    # Tests for the number of diagrams produced for simple BIMSRG cases
    com = argparse.Namespace()
    com.theory, com.order, com.with_3NF = 'BIMSRG', (2, 2, 2), False

    # Use generate_diagrams as a seed
    diagrams = adg.run.generate_diagrams(com, id_gen)
    assert len(diagrams) == 82

    # Test that diagrams are processed and ordered
    diagrams, diag_nbs, _ = adg.run.order_diagrams(diagrams, com)
    assert len(diagrams) == 82
    assert diag_nbs['nb_diags'] == 82
    assert diag_nbs[1] == 10
    assert diag_nbs[2] == 72


def test_print_diags_numbers(capsys):
    """Test the correct output od diagram numbers per type."""
    com = argparse.Namespace()

    # Tests for the number of diagrams produced for simple known cases

    com.theory = 'MBPT'

    diags_nb_per_type = {
        'nb_diags': 0,
        'singles': 1,
        'doubles': 2,
        'triples': 3,
        'quadruples': 4,
        'quintuples+': 5
    }

    adg.run.print_diags_numbers(com, diags_nb_per_type)
    output = capsys.readouterr()
    assert output.out == (
        "Number of connected diagrams:  0\n\n"
        "Valid diagrams: 0\n\n"
        "Singles: 1\n"
        "Doubles: 2\n"
        "Triples: 3\n"
        "Quadruples: 4\n"
        "Quintuples and higher: 5\n\n"
    )

    com.theory, com.with_3NF, com.canonical = 'BMBPT', False, False
    com.nbody_observable = 2

    diags_nb_per_type = {
        'nb_2_hf': 0,
        'nb_2_ehf': 1,
        'nb_2_not_hf': 2,
        'nb_3_hf': 3,
        'nb_3_ehf': 4,
        'nb_3_not_hf': 4,
        'nb_diags': 5,
        'nb_2': 6,
        'nb_3': 7
    }

    adg.run.print_diags_numbers(com, diags_nb_per_type)
    output = capsys.readouterr()
    assert output.out == (
        "Number of connected diagrams:  5\n\n"
        "2N valid diagrams: 6\n"
        "2N energy canonical diagrams: 0\n"
        "2N canonical diagrams for a generic operator only: 1\n"
        "2N non-canonical diagrams: 2\n\n"
    )

    com.theory, com.with_3NF, com.canonical = 'BMBPT', True, False
    com.nbody_observable = 3

    adg.run.print_diags_numbers(com, diags_nb_per_type)
    output = capsys.readouterr()
    assert output.out == (
        "Number of connected diagrams:  5\n\n"
        "2N valid diagrams: 6\n"
        "2N energy canonical diagrams: 0\n"
        "2N canonical diagrams for a generic operator only: 1\n"
        "2N non-canonical diagrams: 2\n\n"
        "3N valid diagrams: 7\n"
        "3N energy canonical diagrams: 3\n"
        "3N canonical diagrams for a generic operator only: 4\n"
        "3N non-canonical diagrams: 4\n\n"
    )

    com = argparse.Namespace()
    com.theory, com.with_3NF, com.canonical = 'BIMSRG', False, False
    com.order = (2, 2, 2)

    diags_nb_per_type = {'nb_diags': 82, 1: 10, 2: 72}

    adg.run.print_diags_numbers(com, diags_nb_per_type)
    output = capsys.readouterr()
    assert output.out == (
        "Number of connected diagrams:  82\n\n"
        "Valid diagrams: 82\n"
        "d_max = 1 diagrams: 10\n"
        "d_max = 2 diagrams: 72\n\n"
    )
