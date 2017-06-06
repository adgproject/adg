#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

import os
import multiprocessing
from datetime import datetime
import shutil
import numpy as np
import networkx as nx
import methods as mth


print "#####################"
print "# Automatic Diagram #"
print "#     Generator     #"
print "#    RDL,JR,PA,MD   #"
print "#####################"

print "Parallel Mode"
num_cores = multiprocessing.cpu_count()
print "There is %i" % num_cores + " core(s) available"
norder = int(raw_input('Order of the diagrams?\n'))
while norder < 2:
    print "Perturbative order too small!"
    norder = int(raw_input('Order of the diagrams?\n'))
theory = raw_input('MBPT or BMBPT?\n').upper()

three_N = False
norm = False
if theory == "BMBPT":
    three_N = raw_input("Include three-body forces? (y/N)").lower() == 'y'
    norm = raw_input(
        "Compute norm kernel instead of operator kernel? (y/N)").lower() == 'y'
if three_N:
    directory = theory + '/Order-%i' % norder + 'with3N'
else:
    directory = theory + '/Order-%i' % norder
if norm:
    directory = directory + '_Norm'
if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(directory+"/Diagrams"):
    os.makedirs(directory+"/Diagrams")


# Start computing everything
print "Running"
start_time = datetime.now()
if theory == "MBPT":
    diagrams = mth.diagram_generation(norder)
elif theory == "BMBPT":
    diagrams = mth.BMBPT_generation(norder, three_N, norm)
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
G = []
for diagram in diagrams:
    G.append(nx.from_numpy_matrix(
        diagram, create_using=nx.MultiDiGraph(), parallel_edges=True))

G1 = []
for diag in G:
    if (nx.number_weakly_connected_components(diag)) == 1:
        G1.append(diag)
G = G1
# Specific check for loop diagrams in BMBPT
if theory == "BMBPT":
    G1 = []
    for diag in G:
        if nx.is_directed_acyclic_graph(diag):
            G1.append(diag)
    G = G1

for diag in G:
    # Account for different status of vertices in operator diagrams
    for node in diag:
        diag.node[node]['operator'] = False
    if (theory == "BMBPT") and not norm:
        diag.node[0]['operator'] = True

# Specific check for topologically identical diagrams in BMBPT
if theory == "BMBPT":
    G = mth.topologically_distinct_diags(G)

numdiag = len(G)
print "Time ellapsed: ", datetime.now() - start_time
print "Number of connected diagrams, ", numdiag

# Ordering the diagrams in a convenient way
if theory == "BMBPT":
    G2 = []
    G3 = []
    G2_HF = []
    G2_EHF = []
    G2_noHF = []
    G3_HF = []
    G3_EHF = []
    G3_noHF = []

    mth.order_2B_or_3B(G, G2, G3)
    mth.order_HF_or_not(G2, G2_HF, G2_EHF, G2_noHF, norm)
    mth.order_HF_or_not(G3, G3_HF, G3_EHF, G3_noHF, norm)

    G = G2_HF + G2_EHF + G2_noHF + G3_HF + G3_EHF + G3_noHF
    nb_2 = len(G2)
    nb_2_HF = len(G2_HF)
    nb_2_EHF = len(G2_EHF)
    nb_2_noHF = len(G2_noHF)
    nb_3 = len(G3)
    nb_3_HF = len(G3_HF)
    nb_3_EHF = len(G3_EHF)
    nb_3_noHF = len(G3_noHF)
    print "\n2N valid diagrams: %i" % nb_2
    print "2N energy canonical diagrams: %i" % nb_2_HF
    if not norm:
        print "2N canonical diagrams for a generic operator only: %i" % nb_2_EHF
    print "2N non-canonical diagrams: %i\n" % nb_2_noHF
    if three_N:
        print "3N valid diagrams: %i" % nb_3
        print "3N energy canonical diagrams: %i" % nb_3_HF
        if not norm:
            print "3N canonical diagrams for a generic operator only: %i" % nb_3_EHF
        print "3N non-canonical diagrams: %i" % nb_3_noHF

# Algebraic expressions:

# Treatment of the algebraic expressions
# To be extended to BMBPT in the future
if theory == "MBPT":
    mat_els = []
    denoms = []
    phases = []
    nedges_eq = []
    for diag in G:
        type_edg = []
        braket = ''
        # Beware of the sign convention !!!
        incidence = - nx.incidence_matrix(diag, oriented=True).todense()
        nrow = diag.number_of_nodes()
        ncol = diag.number_of_edges()
        n_holes = 0
        diffcols = set()
        for col in range(ncol):
            flat = list(incidence[:, col].A1)
            if flat.index(1) < flat.index(-1):
                n_holes += 1
                type_edg.append('h')
            else:
                type_edg.append('p')
            diffcols.add(repr(flat))

        for row in range(nrow):
            ket = ''
            bra = ''
            for col in range(ncol):
                ######### Mtrx Elements ###########
                if (incidence[row, col] == 1):
                    if type_edg[col] == 'h':
                        bra = bra + mth.line_label_h(col)
                    else:
                        bra = bra + mth.line_label_p(col)
                if (incidence[row, col] == -1):
                    if type_edg[col] == 'h':
                        ket = ket + mth.line_label_h(col)
                    else:
                        ket = ket + mth.line_label_p(col)
                ###################################
            braket = braket + '\\braket{'+bra+'|H|'+ket+'}'
        mat_els.append(braket)
        denom = ''
        for row in range(1, nrow):
            denom = denom + '('
            for col in range(ncol):
                val_test = incidence[0:row, col].sum()
                if val_test == 1:
                    if type_edg[col] == 'h':
                        denom = denom + ' +E_' + mth.line_label_h(col)
                    else:
                        denom = denom + ' +E_' + mth.line_label_p(col)
                if val_test == -1:
                    if type_edg[col] == 'h':
                        denom = denom + '-E_' + mth.line_label_h(col)
                    else:
                        denom = denom + '-E_' + mth.line_label_p(col)
            denom = denom + ')'
        if '( +' in denom:
            denom = denom.replace('( +', '(')
        denom = denom.strip(' ')
        denoms.append(denom)
        phases.append('(-1)^{%i' % n_holes + '+l}')
        # print incidence
        eq_lines = np.array(incidence.transpose())
        neq_lines = np.asarray(list(i for i in set(map(tuple, eq_lines))))
        n_sym = len(eq_lines)-len(neq_lines)
        # CAVEAT !!! Valid only for *MBPT*
        nedges_eq.append(2**n_sym)
        # print "After neqlines"
        # Loops


# Treatment of the algebraic expressions
if theory == "BMBPT" and not norm:
    feynman_expressions = []
    diag_expressions = []
    for diag in G:
        # Attribute a qp label to all propagators
        mth.attribute_qp_labels(diag)
        # Determine the numerator corresponding to the diagram
        numerator = mth.extract_numerator(diag)
        # Determine the type of structure we have
        # Create a subgraph without the operator vertex
        testdiag = mth.omega_subgraph(diag)
        # Use the subgraph to determine the structure of the overall graph
        test_adg_subgraphs = mth.has_only_adg_operator_subgraphs(testdiag)
        sink_number = mth.number_of_sinks(diag)
        # Determine the denominator depending on the graph structure
        denominator = ""
        if test_adg_subgraphs:
            for connected_subgraph in nx.weakly_connected_component_subgraphs(testdiag):
                for i in range(len(connected_subgraph)):
                    subgraph_stack = []
                    if len(connected_subgraph) > 1:
                        for j in range(i, len(connected_subgraph)):
                            vertex = nx.dag_longest_path(connected_subgraph)[j]
                            subgraph_stack.append(vertex)
                    elif len(connected_subgraph) == 1:
                        subgraph_stack.append(connected_subgraph.nodes()[0])
                    subdiag = connected_subgraph.subgraph(subgraph_stack)
                    denominator += "(" + mth.extract_denom(diag, subdiag) + ")"
        elif (norder == 4) and (sink_number == 1):
            for i in range(2):
                subgraph_stack = []
                subgraph_stack.append(nx.dag_longest_path(testdiag)[1])
                if i == 0:
                    subgraph_stack.append(nx.dag_longest_path(testdiag)[0])
                else:
                    for vertex_1 in nx.nodes(testdiag):
                        test_vertex = True
                        for vertex_2 in nx.dag_longest_path(testdiag):
                            if vertex_1 == vertex_2:
                                test_vertex = False
                        if test_vertex:
                            subgraph_stack.append(vertex_1)
                subdiag = testdiag.subgraph(subgraph_stack)
                denominator += "(" + mth.extract_denom(diag, subdiag) + ")"
        elif (norder == 4) and (sink_number == 2):
            denominator += "(" + mth.extract_denom(diag, testdiag) + ")"
            for vertex in diag:
                if diag.out_degree(vertex) == 0:
                    denominator += "(" \
                        + mth.extract_denom(diag, diag.subgraph(vertex)) + ")"
        # Determine the integral component in the Feynman expression
        integral = mth.extract_integral(diag)
        # Determine the pre-factor
        prefactor = "(-1)^%i " % (norder - 1)
        if mth.extract_BMBPT_crossing_sign(diag):
            prefactor = "-" + prefactor
        sym_fact = mth.vertex_exchange_sym_factor(diag) \
            + mth.multiplicity_symmetry_factor(diag)
        if sym_fact != "":
            prefactor = "\\frac{" + prefactor + "}{" \
                + sym_fact + "}\\sum_{k_i}"
        else:
            prefactor = prefactor + "\\sum_{k_i}"
        feynman_exp = prefactor + numerator + "\\int_{0}^{\\tau}" \
            + integral + "\n"
        if denominator != "":
            diag_exp = prefactor + "\\frac{ " + numerator \
                + " }{ " + denominator + " }\n"
        else:
            diag_exp = prefactor + numerator + "\n"
        if (norder == 4) and (sink_number == 1) and (not test_adg_subgraphs):
            for vertex in range(norder):
                if diag.out_degree(vertex) == 0:
                    subdiag = diag.subgraph(vertex)
            denominator_abc = mth.extract_denom(diag, subdiag)
            subdiag = mth.omega_subgraph(diag)
            denominator_a = mth.extract_denom(diag, subdiag)
            diag_exp += "\\left[ \\frac{1}{" + denominator_abc \
                + "} + \\frac{1}{" + denominator_a + "} \\right]"
        feynman_expressions.append(feynman_exp)
        diag_expressions.append(diag_exp)


# Writing a feynmp file for each graph
msg = 'Generate diagrams feymanmf instructions?'
pdraw = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdraw:
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    for i in range(0, numdiag):
        diag_name = 'diag_%i' % i
        mth.feynmf_generator(G[i], theory, diag_name)
        shutil.move(diag_name + '.tex',
                    directory + "/Diagrams/" + diag_name + '.tex')


msg = 'Include diagrams in tex?'
pdiag = raw_input("%s (y/N) " % msg).lower() == 'y'

# Write everything down in a nice LaTeX file
header = "\\documentclass[10pt,a4paper]{article}\n" \
    + "\\usepackage[utf8]{inputenc}\n" \
    + "\\usepackage{braket}\n\\usepackage{graphicx}\n" \
    + "\\usepackage[english]{babel}\n\\usepackage{amsmath}\n" \
    + "\\usepackage{amsfonts}\n\\usepackage{amssymb}\n"
if pdiag:
    header = header + "\\usepackage[force]{feynmp-auto}\n"
land = False
if norder > 3:
    msg = 'Expressions may be long, rotate pdf?'
    land = raw_input("%s (y/N) " % msg).lower() == 'y'
if land:
    header = header + "\\usepackage[landscape]{geometry}\n"

header = header \
    + "\\title{Diagrams and algebraic expressions at order %i" % norder \
    + " in " + theory + "}\n" \
    + "\\author{RDL, JR, PA, MD, AT}\n"
latex_file = open(directory + '/result.tex', 'w')
latex_file.write(header)
begdoc = "\\begin{document}\n"
enddoc = "\\end{document}"
begeq = "\\begin{equation}\n"
endeq = "\\end{equation}\n"
latex_file.write(begdoc)
latex_file.write("\\maketitle\n")
latex_file.write("\\graphicspath{{Diagrams/}}")

if theory == "BMBPT":
    latex_file.write("Valid diagrams: %i\n\n" % numdiag)
    latex_file.write("2N valid diagrams: %i\n\n" % nb_2)
    latex_file.write("2N canonical diagrams for the energy: %i\n\n" % nb_2_HF)
    if not norm:
        latex_file.write("2N canonical diagrams for a generic operator only: %i\n\n" % nb_2_EHF)
    latex_file.write("2N non-canonical diagrams: %i\n\n" % nb_2_noHF)
    if three_N:
        latex_file.write("3N valid diagrams: %i\n\n" % nb_3)
        latex_file.write("3N canonical diagrams for the energy: %i\n\n" % nb_3_HF)
        if not norm:
            latex_file.write("3N canonical diagrams for a generic operator only: %i\n\n" % nb_3_EHF)
        latex_file.write("3N non-canonical diagrams: %i\n\n" % nb_3_noHF)
    latex_file.write("\\section{Two-body diagrams}\n")
    latex_file.write("\\subsection{Two-body energy canonical diagrams}\n")

for i_diag in range(0, numdiag):
    if theory == "BMBPT":
        if (i_diag == nb_2_HF) and (not norm):
            latex_file.write("\\subsection{Two-body canonical diagrams for a generic operator only}\n")
        elif i_diag == nb_2_HF + nb_2_EHF:
            latex_file.write("\\subsection{Two-body non-canonical diagrams}\n")
        if three_N:
            if i_diag == nb_2:
                latex_file.write("\\section{Three-body diagrams}\n")
                latex_file.write("\\subsection{Three-body energy canonical diagrams}\n")
            elif (i_diag == nb_2 + nb_3_HF) and (not norm):
                latex_file.write("\\subsection{Three-body canonical diagrams for a generic operator only}\n")
            elif i_diag == nb_2 + nb_3_HF + nb_3_EHF:
                latex_file.write("\\subsection{Three-body non-canonical diagrams}\n")
        latex_file.write("Diagram %i:\n" % (i_diag+1))
        if not norm:
            diag_exp = diag_expressions[i_diag]
            feynman_exp = feynman_expressions[i_diag]
            latex_file.write(begeq)
            latex_file.write(feynman_exp)
            latex_file.write(endeq)
    elif theory == "MBPT":
        diag_exp = "\\dfrac{1}{%i}" % nedges_eq[i_diag] + phases[i_diag] \
            + "\\sum{\\dfrac{" + mat_els[i_diag] + "}{" \
            + denoms[i_diag] + "}}\n"
    if not norm:
        latex_file.write(begeq)
        latex_file.write(diag_exp)
        latex_file.write(endeq)
    if pdiag and pdraw:
        latex_file.write('\n\\begin{center}\n')
        diag_file = open(directory+"/Diagrams/diag_%i.tex" % i_diag)
        latex_file.write(diag_file.read())
        latex_file.write('\\end{center}\n\n')
latex_file.write(enddoc)
latex_file.close()

msg = 'Compile pdf?'
pdfcompile = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdfcompile:
    os.chdir(directory)
    os.system("pdflatex -shell-escape result.tex")
    if pdiag:
        # Second compilation needed
        os.system("pdflatex -shell-escape result.tex")
        # Get rid of undesired feynmp files to keep a clean directory
        for i_diag in range(0, numdiag):
            os.unlink("diag_%i.1" % i_diag)
            os.unlink("diag_%i.mp" % i_diag)
            os.unlink("diag_%i.log" % i_diag)
    print "Result saved in "+directory + '/result.pdf'
