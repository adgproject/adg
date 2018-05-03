#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

from datetime import datetime
import shutil
import cProfile
import pstats
import StringIO
import general_routines as gen
import bmbpt
import mbpt
import time_structure as tsd


run_commands = gen.parse_command_line()

print "#####################"
print "# Automatic Diagram #"
print "#     Generator     #"
print "#    RDL,JR,PA,MD   #"
print "#####################"

if run_commands.interactive:
    run_commands = gen.interactive_interface(run_commands)

directory = gen.attribute_directory(run_commands)

pr = cProfile.Profile()
pr.enable()

# Start computing everything
print "Running"
START_TIME = datetime.now()

diagrams = gen.generate_diagrams(run_commands)

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
print "Number of connected diagrams, ", diags_per_type['nb_diags']

if run_commands.theory == "BMBPT":
    print "\n2N valid diagrams: %i" % diags_per_type['nb_2']
    print "2N energy canonical diagrams: %i" % diags_per_type['nb_2_hf']
    if not run_commands.norm:
        print "2N canonical diagrams for a generic operator only: %i" \
            % diags_per_type['nb_2_ehf']
    print "2N non-canonical diagrams: %i\n" % diags_per_type['nb_2_not_hf']
    if run_commands.with_three_body:
        print "3N valid diagrams: %i" % diags_per_type['nb_3']
        print "3N energy canonical diagrams: %i" % diags_per_type['nb_3_hf']
        if not run_commands.norm:
            print "3N canonical diagrams for a generic operator only: %i" \
                % diags_per_type['nb_3_ehf']
        print "3N non-canonical diagrams: %i" % diags_per_type['nb_3_not_hf']
elif run_commands.theory == "MBPT":
    print "\nValid diagrams: %i\n" % diags_per_type['nb_diags']
    print "Singles: %i" % diags_per_type['singles']
    print "Doubles: %i" % diags_per_type['doubles']
    print "Triples: %i" % diags_per_type['triples']
    print "Quadruples: %i" % diags_per_type['quadruples']
    print "Quintuples and higher: %i" % diags_per_type['quintuples+']


pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
ps.dump_stats("stats.dat")

# Writing a feynmp file for each graph
if run_commands.draw_diags:
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    gen.create_feynmanmp_files(diagrams, run_commands.theory,
                               directory, 'diag')
    if run_commands.draw_tsds:
        gen.create_feynmanmp_files(diagrams_time, run_commands.theory,
                                   directory, 'time')

# Write everything down in a nice LaTeX file
latex_file = open(directory + '/result.tex', 'w')
gen.write_file_header(latex_file, run_commands.draw_diags, run_commands.order,
                      run_commands.theory)

if run_commands.theory == "BMBPT":
    bmbpt.write_header(latex_file, run_commands.with_three_body,
                       run_commands.norm, diags_per_type)
elif run_commands.theory == "MBPT":
    mbpt.write_header(latex_file, diags_per_type)

latex_file.write("\\tableofcontents\n\n")

if run_commands.theory == "BMBPT" and run_commands.draw_tsds:
    tsd.write_section(latex_file, directory,
                      run_commands.draw_diags, run_commands.draw_diags,
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
    gen.compile_and_clean(directory, run_commands.draw_diags)
