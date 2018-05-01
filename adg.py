#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

import os
from datetime import datetime
import shutil
import cProfile
import pstats
import StringIO
import argparse
import numpy as np
import networkx as nx
import general_routines as gen
import bmbpt
import mbpt
import time_structure as tsd

parser = argparse.ArgumentParser(description="Automatic Diagram Generator\n\n"
                                 + "Generates diagrams at a given order for a "
                                 + "range of many-body formalisms")
parser.add_argument("-o", "--order", type=int, choices=range(2, 10),
                    help="order of the diagrams (>=2)")
parser.add_argument("-t", "--theory", type=str, choices=['MBPT', 'BMBPT'],
                    help="theory of interest: MBPT or BMBPT")
parser.add_argument("-i", "--interactive", action="store_true",
                    help="execute ADG in interactive mode")
parser.add_argument("-n", "--norm", action="store_true",
                    help="study norm BMBPT diagrams instead of operator ones")
parser.add_argument("-3N", "--with_three_body", action="store_true",
                    help="use two and three-body forces for BMBPT diagrams")
parser.add_argument("-dt", "--draw_tsds", action="store_true",
                    help="draw Time-Structure Diagrams (BMBPT)")
parser.add_argument("-d", "--draw_diags", action="store_true",
                    help="draw the diagrms using FeynMF")
parser.add_argument("-c", "--compile", action="store_true",
                    help="compile the LaTeX output file with PDFLaTeX")
parser.add_argument("-cd", "--cd_output", action="store_true",
                    help="produce output for C. Drischler's framework (MBPT)")
args = parser.parse_args()

if (not args.interactive) and ((args.order is None) or (args.theory is None)):
    print "\nPlease either run the interactive mode, or the batch mode by"
    print "providing the theory and the order for the desired diagrams.\n"
    print "Use 'python2 adg.py -h' for help.\n"
    exit()

print "#####################"
print "# Automatic Diagram #"
print "#     Generator     #"
print "#    RDL,JR,PA,MD   #"
print "#####################"

if args.interactive:
    args.order = int(raw_input('Order of the diagrams?\n'))
    while args.order < 2:
        print "Perturbative order too small!"
        args.order = int(raw_input('Order of the diagrams?\n'))
    args.theory = raw_input('MBPT or BMBPT?\n').upper()

    if args.theory == "BMBPT":
        args.with_three_body = raw_input(
            "Include three-body forces? (y/N)").lower() == 'y'
        args.norm = raw_input(
            "Compute norm diagrams, not operator ones? (y/N)").lower() == 'y'
        args.draw_tsds = raw_input(
            "Draw time-structure diagrams? (y/N)").lower() == 'y'

    args.draw_diags = raw_input(
        "Generate diagrams FeynMF instructions in TeX file? (y/N) ").lower() == 'y'

    PDIAG = raw_input("Include diagrams in tex? (y/N) ").lower() == 'y'
    if args.theory == "MBPT":
        args.cd_output = raw_input(
            "Produce a CD output file? (y/N) ").lower() == 'y'
    args.compile = raw_input("Compile pdf? (y/N) ").lower() == 'y'

THEORY = args.theory
THREE_N = args.with_three_body
NORM = args.norm
norder = args.order
WRITE_TIME = args.draw_tsds
DO_COMPILE = args.compile
CD_OUTPUT = args.cd_output
PDRAW = args.draw_diags
PDIAG = PDRAW

directory = '%s/Order-%i' % (THEORY, norder)
if THREE_N:
    directory += 'with3N'
if NORM:
    directory += '_NORM'
if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(directory+"/Diagrams"):
    os.makedirs(directory+"/Diagrams")

pr = cProfile.Profile()
pr.enable()

# Start computing everything
print "Running"
START_TIME = datetime.now()
if THEORY == "MBPT":
    diagrams = mbpt.diagrams_generation(norder)
elif THEORY == "BMBPT":
    diagrams = bmbpt.diagrams_generation(norder, THREE_N)
else:
    print "Invalid THEORY"
print "Number of possible diagrams, ", len(diagrams)

with open(directory+"/Diagrams.list", "w") as f:
    for idx, diagram in enumerate(diagrams):
        f.write("Diagram n: %i\n" % (idx + 1))
        np.savetxt(f, diagram, fmt='%d')
        f.write("\n")

# Graph part (computing, writing, drawing)
G = [nx.from_numpy_matrix(diagram, create_using=nx.MultiDiGraph(),
                          parallel_edges=True) for diagram in diagrams]

for i_diag in xrange(len(G)-1, -1, -1):
    if (nx.number_weakly_connected_components(G[i_diag])) != 1:
        del G[i_diag]

# Specific check for loop diagrams in BMBPT
if THEORY == "BMBPT":
    for i_diag in xrange(len(G)-1, -1, -1):
        if not nx.is_directed_acyclic_graph(G[i_diag]):
            del G[i_diag]

gen.label_vertices(G, THEORY, NORM)

if THEORY == 'BMBPT':
    diagrams = [bmbpt.BmbptFeynmanDiagram(graph, NORM, ind)
                for ind, graph in enumerate(G)]
elif THEORY == 'MBPT':
    diagrams = [mbpt.MbptDiagram(graph, ind) for ind, graph in enumerate(G)]

# Ordering the diagrams in a convenient way and checking them for doubles
if THEORY == "BMBPT":
    diagrams, NB_2_HF, NB_2_EHF, NB_2_NOT_HF, \
        NB_3_HF, NB_3_EHF, NB_3_NOT_HF = bmbpt.order_diagrams(diagrams)
    NB_2 = NB_2_HF + NB_2_EHF + NB_2_NOT_HF
    NB_3 = NB_3_HF + NB_3_EHF + NB_3_NOT_HF

elif THEORY == "MBPT":
    diagrams, NB_SINGLES, NB_DOUBLES, NB_TRIPLES, NB_QUADRUPLES, \
        NB_QUINTUPLES_AND_HIGHER = mbpt.order_diagrams(diagrams)
    mbpt.attribute_conjugate(diagrams)

NUMDIAG = len(diagrams)

# Treatment of the algebraic expressions
if THEORY == "BMBPT" and not NORM:

    diagrams_time = [tsd.TimeStructureDiagram(diagram, diagram.tags[0])
                     for diagram in diagrams]

    diagrams_time, nb_tree_TSDs = bmbpt.treat_TSDs(diagrams_time)

    bmbpt.produce_expressions(diagrams, diagrams_time)

else:
    diagrams_time = []

print "Time ellapsed: ", datetime.now() - START_TIME
print "Number of connected diagrams, ", NUMDIAG

if THEORY == "BMBPT":
    print "\n2N valid diagrams: %i" % NB_2
    print "2N energy canonical diagrams: %i" % NB_2_HF
    if not NORM:
        print "2N canonical diagrams for a generic operator only: %i" \
            % NB_2_EHF
    print "2N non-canonical diagrams: %i\n" % NB_2_NOT_HF
    if THREE_N:
        print "3N valid diagrams: %i" % NB_3
        print "3N energy canonical diagrams: %i" % NB_3_HF
        if not NORM:
            print "3N canonical diagrams for a generic operator only: %i" \
                % NB_3_EHF
        print "3N non-canonical diagrams: %i" % NB_3_NOT_HF
elif THEORY == "MBPT":
    print "\nValid diagrams: %i\n" % NUMDIAG
    print "Singles: %i" % NB_SINGLES
    print "Doubles: %i" % NB_DOUBLES
    print "Triples: %i" % NB_TRIPLES
    print "Quadruples: %i" % NB_QUADRUPLES
    print "Quintuples and higher: %i" % NB_QUINTUPLES_AND_HIGHER


pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
ps.dump_stats("stats.dat")

# Writing a feynmp file for each graph
if PDRAW:
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    gen.create_feynmanmp_files(diagrams, THEORY, directory, 'diag')
    if WRITE_TIME:
        gen.create_feynmanmp_files(diagrams_time, THEORY, directory, 'time')

# Write everything down in a nice LaTeX file
latex_file = open(directory + '/result.tex', 'w')
gen.write_file_header(latex_file, PDIAG, norder, THEORY)

if THEORY == "BMBPT":
    bmbpt.write_header(latex_file, NUMDIAG, THREE_N, NORM, NB_2_HF, NB_2_EHF,
                       NB_2_NOT_HF, NB_3_HF, NB_3_EHF, NB_3_NOT_HF)
elif THEORY == "MBPT":
    mbpt.write_header(latex_file, NUMDIAG, NB_SINGLES, NB_DOUBLES, NB_TRIPLES,
                      NB_QUADRUPLES, NB_QUINTUPLES_AND_HIGHER)

latex_file.write("\\tableofcontents\n\n")

if THEORY == "BMBPT" and WRITE_TIME:
    tsd.write_section(latex_file, directory, PDIAG, PDRAW, diagrams_time,
                      nb_tree_TSDs)
for diag in diagrams:
    if THEORY == "BMBPT":
        bmbpt.write_section(latex_file, diag.tags[0], THREE_N, NORM,
                            NB_2, NB_2_HF, NB_2_EHF, NB_3_HF, NB_3_EHF)
        latex_file.write("\\paragraph{Diagram %i:}\n" % (diag.tags[0] + 1))
        if not NORM:
            bmbpt.write_diag_exps(latex_file, diag, norder)
    elif THEORY == "MBPT":
        mbpt.write_section(latex_file, diag.tags[0], NB_SINGLES, NB_DOUBLES,
                           NB_TRIPLES, NB_QUADRUPLES)
        latex_file.write("\\paragraph{Diagram %i:}\n" % (diag.tags[0] + 1))
        if diag.complex_conjugate >= 0:
            latex_file.write("Complex conjugate diagram: %i\n"
                             % (diag.complex_conjugate + 1))
        mbpt.write_diag_exp(latex_file, diag)

    if PDIAG and PDRAW:
        diag.write_graph(latex_file, directory, WRITE_TIME)

    if THEORY == 'BMBPT' and WRITE_TIME:
        diag.write_tsd_info(diagrams_time, latex_file)

latex_file.write("\\end{document}")
latex_file.close()

# Produce an output adapted to Christian Drischler's format
if CD_OUTPUT:
    mbpt.print_cd_output(directory, diagrams)

if DO_COMPILE:
    gen.compile_and_clean(directory, PDIAG)
