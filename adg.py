#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

from datetime import datetime
import cProfile
import pstats
import StringIO
import run_routines as run
import general_routines as gen
import bmbpt
import mbpt
import time_structure as tsd


run_commands = run.parse_command_line()

print "#####################"
print "# Automatic Diagram #"
print "#     Generator     #"
print "#    RDL,JR,PA,MD   #"
print "#####################"

if run_commands.interactive:
    run_commands = run.interactive_interface(run_commands)

directory = run.attribute_directory(run_commands)

pr = cProfile.Profile()
pr.enable()

# Start computing everything
print "Running"
START_TIME = datetime.now()

diagrams = run.generate_diagrams(run_commands)

# Ordering the diagrams in a convenient way and checking them for doubles
if run_commands.theory == "BMBPT":
    diagrams, diags_per_type = bmbpt.order_diagrams(diagrams)

elif run_commands.theory == "MBPT":
    diagrams, diags_per_type = mbpt.order_diagrams(diagrams)
    mbpt.attribute_conjugate(diagrams)

# Treatment of the algebraic expressions
if run_commands.theory == "BMBPT" and not run_commands.norm:

    diagrams_time = [tsd.TimeStructureDiagram(diagram, diagram.tags[0])
                     for diagram in diagrams]

    diagrams_time, nb_tree_TSDs = bmbpt.treat_TSDs(diagrams_time)

    bmbpt.produce_expressions(diagrams, diagrams_time)

else:
    diagrams_time = []

print "Time ellapsed: ", datetime.now() - START_TIME

run.print_diags_numbers(run_commands, diags_per_type)

pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
ps.dump_stats("stats.dat")

# Writing a feynmp file for each graph
if run_commands.draw_diags:
    run.prepare_drawing_instructions(directory, run_commands,
                                     diagrams, diagrams_time)

# Write everything down in a nice LaTeX file
latex_file = open(directory + '/result.tex', 'w')

run.write_file_header(latex_file, run_commands, diags_per_type)

if run_commands.theory == "BMBPT" and run_commands.draw_tsds:
    tsd.write_section(latex_file, directory, run_commands.draw_diags,
                      diagrams_time, nb_tree_TSDs)
for diag in diagrams:
    if run_commands.theory == "BMBPT":
        bmbpt.write_section(latex_file, diag.tags[0],
                            run_commands.with_three_body, run_commands.norm,
                            diags_per_type)
        latex_file.write("\\paragraph{Diagram %i:}\n" % (diag.tags[0] + 1))
        if not run_commands.norm:
            bmbpt.write_diag_exps(latex_file, diag, run_commands.order)
    elif run_commands.theory == "MBPT":
        mbpt.write_section(latex_file, diag.tags[0], diags_per_type)
        latex_file.write("\\paragraph{Diagram %i:}\n" % (diag.tags[0] + 1))
        if diag.complex_conjugate >= 0:
            latex_file.write("Complex conjugate diagram: %i\n"
                             % (diag.complex_conjugate + 1))
        mbpt.write_diag_exp(latex_file, diag)

    if run_commands.draw_diags:
        diag.write_graph(latex_file, directory, run_commands.draw_tsds)

    if run_commands.theory == 'BMBPT' and run_commands.draw_tsds:
        diag.write_tsd_info(diagrams_time, latex_file)

latex_file.write("\\end{document}")
latex_file.close()

# Produce an output adapted to Christian Drischler's format
if run_commands.cd_output:
    mbpt.print_cd_output(directory, diagrams)

if run_commands.compile:
    run.compile_and_clean(directory, run_commands.draw_diags)
