#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

import os
from datetime import datetime
import shutil
import cProfile
import pstats
import StringIO
import numpy as np
import networkx as nx
import general_routines as gen
import bmbpt
import mbpt
import time_structure as tsd


print "#####################"
print "# Automatic Diagram #"
print "#     Generator     #"
print "#    RDL,JR,PA,MD   #"
print "#####################"


norder = int(raw_input('Order of the diagrams?\n'))
while norder < 2:
    print "Perturbative order too small!"
    norder = int(raw_input('Order of the diagrams?\n'))
theory = raw_input('MBPT or BMBPT?\n').upper()

three_N = False
norm = False
write_time = False
if theory == "BMBPT":
    three_N = raw_input("Include three-body forces? (y/N)").lower() == 'y'
    norm = raw_input(
        "Compute norm kernel instead of operator kernel? (y/N)").lower() == 'y'
    write_time = raw_input(
        "Draw time-structure diagrams? (y/N)").lower() == 'y'
directory = '%s/Order-%i' % (theory, norder)
if three_N:
    directory += 'with3N'
if norm:
    directory += '_Norm'
if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(directory+"/Diagrams"):
    os.makedirs(directory+"/Diagrams")

pr = cProfile.Profile()
pr.enable()

# Start computing everything
print "Running"
start_time = datetime.now()
if theory == "MBPT":
    diagrams = mbpt.diagram_generation(norder)
elif theory == "BMBPT":
    diagrams = bmbpt.BMBPT_generation(norder, three_N)
else:
    print "Invalid theory"
numdiag = len(diagrams)
print "Number of possible diagrams, ", numdiag

i = 0
with open(directory+"/Diagrams.list", "w") as f:
    for diagram in diagrams:
        f.write("Diagram n: %i" % i)
        np.savetxt(f, diagram)
        f.write("\n")
        i += 1

# Graph part (computing, writing, drawing)
G = [nx.from_numpy_matrix(diagram, create_using=nx.MultiDiGraph(),
                          parallel_edges=True) for diagram in diagrams]

for i_diag in xrange(len(G)-1, -1, -1):
    if (nx.number_weakly_connected_components(G[i_diag])) != 1:
        del G[i_diag]

# Specific check for loop diagrams in BMBPT
if theory == "BMBPT":
    for i_diag in xrange(len(G)-1, -1, -1):
        if not nx.is_directed_acyclic_graph(G[i_diag]):
            del G[i_diag]

gen.label_vertices(G, theory, norm)

if theory == 'BMBPT':
    diagrams = [bmbpt.BmbptFeynmanDiagram(graph, norm, ind)
                for ind, graph in enumerate(G)]
elif theory == 'MBPT':
    diagrams = [mbpt.MbptDiagram(graph, ind) for ind, graph in enumerate(G)]

# Ordering the diagrams in a convenient way and checking them for doubles
if theory == "BMBPT":
    diagrams2HF = []
    diagrams2EHF = []
    diagrams2noHF = []
    diagrams3HF = []
    diagrams3EHF = []
    diagrams3noHF = []

    for i_diag in xrange(len(diagrams)-1, -1, -1):
        if diagrams[i_diag].two_or_three_body == 2:
            if diagrams[i_diag].HF_type == "HF":
                diagrams2HF.append(diagrams[i_diag])
            elif diagrams[i_diag].HF_type == "EHF":
                diagrams2EHF.append(diagrams[i_diag])
            elif diagrams[i_diag].HF_type == "noHF":
                diagrams2noHF.append(diagrams[i_diag])
        elif diagrams[i_diag].two_or_three_body == 3:
            if diagrams[i_diag].HF_type == "HF":
                diagrams3HF.append(diagrams[i_diag])
            elif diagrams[i_diag].HF_type == "EHF":
                diagrams3EHF.append(diagrams[i_diag])
            elif diagrams[i_diag].HF_type == "noHF":
                diagrams3noHF.append(diagrams[i_diag])
        del diagrams[i_diag]

    gen.topologically_distinct_diagrams(diagrams2HF)
    gen.topologically_distinct_diagrams(diagrams2EHF)
    gen.topologically_distinct_diagrams(diagrams2noHF)
    gen.topologically_distinct_diagrams(diagrams3HF)
    gen.topologically_distinct_diagrams(diagrams3EHF)
    gen.topologically_distinct_diagrams(diagrams3noHF)

    diagrams = diagrams2HF + diagrams2EHF + diagrams2noHF \
        + diagrams3HF + diagrams3EHF + diagrams3noHF
    for ind, diagram in enumerate(diagrams):
        diagram.tags[0] = ind
    nb_2_HF = len(diagrams2HF)
    nb_2_EHF = len(diagrams2EHF)
    nb_2_noHF = len(diagrams2noHF)
    nb_2 = nb_2_HF + nb_2_EHF + nb_2_noHF
    nb_3_HF = len(diagrams3HF)
    nb_3_EHF = len(diagrams3EHF)
    nb_3_noHF = len(diagrams3noHF)
    nb_3 = nb_3_HF + nb_3_EHF + nb_3_noHF

elif theory == "MBPT":
    singles = []
    doubles = []
    triples = []
    quadruples = []
    quintuples_and_higher = []

    for i_diag in xrange(len(diagrams)-1, -1, -1):
        if diagrams[i_diag].excitation_level == 1:
            singles.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level == 2:
            doubles.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level == 3:
            triples.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level == 4:
            quadruples.append(diagrams[i_diag])
        elif diagrams[i_diag].excitation_level >= 5:
            quintuples_and_higher.append(diagrams[i_diag])
        else:
            print "Zero or negative excitation level!\n"
            exit()
        del diagrams[i_diag]

    diagrams = singles + doubles + triples + quadruples \
        + quintuples_and_higher
    nb_singles = len(singles)
    nb_doubles = len(doubles)
    nb_triples = len(triples)
    nb_quadruples = len(quadruples)
    nb_quintuples_and_higher = len(quintuples_and_higher)

numdiag = len(diagrams)

# Treatment of the algebraic expressions
if theory == "BMBPT" and not norm:

    # Production of the time-structure diagrams

    diagrams_time = [tsd.TimeStructureDiagram(diagram, diagram.tags[0])
                     for diagram in diagrams]

    tree_TSDs = []
    for i_diag in xrange(len(diagrams_time)-1, -1, -1):
        if diagrams_time[i_diag].is_tree:
            tree_TSDs.append(diagrams_time[i_diag])
            del diagrams_time[i_diag]

    gen.topologically_distinct_diagrams(tree_TSDs)
    nb_tree_TSDs = len(tree_TSDs)

    gen.topologically_distinct_diagrams(diagrams_time)
    diagrams_time = tree_TSDs + diagrams_time

    for index, t_diag in enumerate(diagrams_time):
        t_diag.tags.insert(0, index)
        if not t_diag.is_tree:
            t_diag.equivalent_trees = tsd.treat_cycles(t_diag.graph)
            t_diag.expr = " + ".join("\\frac{1}{%s}"
                                     % tsd.tree_time_structure_den(graph)
                                     for graph
                                     in t_diag.equivalent_trees)

    bmbpt.produce_expressions(diagrams, diagrams_time)

else:
    diagrams_time = []

print "Time ellapsed: ", datetime.now() - start_time
print "Number of connected diagrams, ", numdiag

if theory == "BMBPT":
    print "\n2N valid diagrams: %i" % nb_2
    print "2N energy canonical diagrams: %i" % nb_2_HF
    if not norm:
        print "2N canonical diagrams for a generic operator only: %i" \
            % nb_2_EHF
    print "2N non-canonical diagrams: %i\n" % nb_2_noHF
    if three_N:
        print "3N valid diagrams: %i" % nb_3
        print "3N energy canonical diagrams: %i" % nb_3_HF
        if not norm:
            print "3N canonical diagrams for a generic operator only: %i" \
                % nb_3_EHF
        print "3N non-canonical diagrams: %i" % nb_3_noHF
elif theory == "MBPT":
    print "\nValid diagrams: %i\n" % numdiag
    print "Singles: %i" % nb_singles
    print "Doubles: %i" % nb_doubles
    print "Triples: %i" % nb_triples
    print "Quadruples: %i" % nb_quadruples
    print "Quintuples and higher: %i" % nb_quintuples_and_higher


pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
ps.dump_stats("stats.dat")

# Writing a feynmp file for each graph
msg = 'Generate diagrams feymanmf instructions?'
pdraw = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdraw:
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    gen.create_feynmanmp_files(diagrams, theory, directory, 'diag')
    if write_time:
        gen.create_feynmanmp_files(diagrams_time, theory, directory, 'time')

pdiag = raw_input("Include diagrams in tex? (y/N) ").lower() == 'y'

# Write everything down in a nice LaTeX file
latex_file = open(directory + '/result.tex', 'w')
gen.write_file_header(directory, latex_file, pdiag, norder, theory)

if theory == "BMBPT":
    bmbpt.write_BMBPT_header(latex_file, numdiag, three_N, norm, nb_2_HF,
                             nb_2_EHF, nb_2_noHF, nb_3_HF, nb_3_EHF, nb_3_noHF)
elif theory == "MBPT":
    mbpt.write_MBPT_header(latex_file, numdiag, nb_singles, nb_doubles,
                           nb_triples, nb_quadruples,
                           nb_quintuples_and_higher)

latex_file.write("\\tableofcontents\n\n")

if theory == "BMBPT" and write_time:
    tsd.write_time_diagrams_section(latex_file, directory, pdiag, pdraw,
                                    diagrams_time, nb_tree_TSDs)
for diag_idx, diag in enumerate(diagrams):
    if theory == "BMBPT":
        bmbpt.write_BMBPT_section(latex_file, diag.tags[0], three_N, norm,
                                  nb_2, nb_2_HF, nb_2_EHF, nb_3_HF, nb_3_EHF)
        latex_file.write("\\paragraph{Diagram %i:}\n" % (diag.tags[0] + 1))
        if not norm:
            bmbpt.write_diag_exps(latex_file, diag, norder)
    elif theory == "MBPT":
        mbpt.write_MBPT_section(latex_file, diag_idx, nb_singles, nb_doubles,
                                nb_triples, nb_quadruples,
                                nb_quintuples_and_higher)
        mbpt.write_diag_exp(latex_file, diag)
    if pdiag and pdraw:
        latex_file.write('\n\\begin{center}\n')
        gen.draw_diagram(directory, latex_file, diag.tags[0], 'diag')
        if write_time:
            latex_file.write('\\hspace{10pt} $\\rightarrow$ \\hspace{10pt} T%i:'
                             % (diag.time_tag + 1))
            gen.draw_diagram(directory, latex_file, diag.time_tag, 'time')
        latex_file.write('\n\\end{center}\n\n')
    if theory == 'BMBPT' and write_time:
        for tdiag in diagrams_time:
            if diag.time_tag == tdiag.tags[0]:
                time_diag = tdiag
                break
        latex_file.write("\\begin{equation}\n\\text{T}%i = " % (diag.time_tag
                                                                + 1)
                         + "%s\\end{equation}\n" % time_diag.expr)
        bmbpt.write_vertices_values(latex_file,
                                    diag,
                                    time_diag.perms[diag.tags[0]])

latex_file.write("\\end{document}")
latex_file.close()

# Produce an output adapted to Christian Drischler's format
if theory == "MBPT":
    if raw_input("Produce a CD output file? (y/N) ").lower() == 'y':
        mbpt.print_CD_output(directory, diagrams)

if raw_input("Compile pdf? (y/N) ").lower() == 'y':
    gen.compile_and_clean(directory, pdiag, diagrams, write_time,
                          diagrams_time)
