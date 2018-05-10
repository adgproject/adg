"""Main routine of the Automated Diagram Generator."""

from datetime import datetime
import adg.run
import adg.bmbpt
import adg.mbpt
import adg.tsd


def main():
    """Launch the ADG program."""
    run_commands = adg.run.parse_command_line()

    print "#####################"
    print "# Automatic Diagram #"
    print "#     Generator     #"
    print "#    RDL,JR,PA,MD   #"
    print "#####################"

    if run_commands.interactive:
        run_commands = adg.run.interactive_interface(run_commands)

    directory = adg.run.attribute_directory(run_commands)

    # Start computing everything
    print "Running"
    start_time = datetime.now()

    diagrams = adg.run.generate_diagrams(run_commands)

    # Ordering the diagrams in a convenient way and checking them for doubles
    diagrams, diags_per_type = adg.run.order_diagrams(diagrams,
                                                      run_commands)

    # Produce TSD for the expressions of BMBPT diagrams
    if run_commands.theory == "BMBPT" and not run_commands.norm:

        diagrams_time = [adg.tsd.TimeStructureDiagram(diagram, diagram.tags[0])
                         for diagram in diagrams]

        diagrams_time, nb_tree_tsds = adg.bmbpt.treat_tsds(diagrams_time)

        adg.bmbpt.produce_expressions(diagrams, diagrams_time)

    else:
        diagrams_time = []

    print "Time ellapsed: ", datetime.now() - start_time

    adg.run.print_diags_numbers(run_commands, diags_per_type)

    # Writing a feynmp file for each graph
    if run_commands.draw_diags:
        adg.run.prepare_drawing_instructions(directory, run_commands,
                                             diagrams, diagrams_time)

    # Write everything down in a nice LaTeX file
    latex_file = open(directory + '/result.tex', 'w')

    adg.run.write_file_header(latex_file, run_commands, diags_per_type)

    if run_commands.theory == "BMBPT" and run_commands.draw_tsds:
        adg.tsd.write_section(latex_file, directory, run_commands.draw_diags,
                              diagrams_time, nb_tree_tsds)
    for diag in diagrams:
        if run_commands.theory == "BMBPT":
            adg.bmbpt.write_section(latex_file, diag,
                                    run_commands, diags_per_type)
        elif run_commands.theory == "MBPT":
            adg.mbpt.write_section(latex_file, diag, diags_per_type)

        if run_commands.draw_diags:
            diag.write_graph(latex_file, directory, run_commands.draw_tsds)

        if run_commands.theory == 'BMBPT' and run_commands.draw_tsds:
            diag.write_tsd_info(diagrams_time, latex_file)

    latex_file.write("\\end{document}")
    latex_file.close()

    # Produce an output adapted to Christian Drischler's format
    if run_commands.cd_output:
        adg.mbpt.print_cd_output(directory, diagrams)

    if run_commands.compile:
        adg.run.compile_and_clean(directory, run_commands.draw_diags)


if __name__ == "__main__":
    main()