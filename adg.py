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
use_parallel = raw_input("Use parallel processing? (y/N)").lower() == 'y'
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
    write_time = raw_input("Draw time-structure diagrams? (y/N)").lower() == 'y'
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


# Ordering the diagrams in a convenient way and checking them for doubles
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

    if use_parallel:
        if three_N:
            nb_procs_max = 6
        else:
            nb_procs_max = 3
        nb_processes = min(num_cores-1, nb_procs_max)
        pool = multiprocessing.Pool(nb_processes)
        r1 = pool.apply_async(mth.topologically_distinct_diags, (G2_HF, ))
        r2 = pool.apply_async(mth.topologically_distinct_diags, (G2_EHF, ))
        r3 = pool.apply_async(mth.topologically_distinct_diags, (G2_noHF, ))
        if three_N:
            r4 = pool.apply_async(mth.topologically_distinct_diags, (G3_HF, ))
            r5 = pool.apply_async(mth.topologically_distinct_diags, (G3_EHF, ))
            r6 = pool.apply_async(mth.topologically_distinct_diags, (G3_noHF, ))
        G2_HF = r1.get()
        G2_EHF = r2.get()
        G2_noHF = r3.get()
        if three_N:
            G3_HF = r4.get()
            G3_EHF = r5.get()
            G3_noHF = r6.get()
        pool.close()
        pool.join()

    else:
        G2_HF = mth.topologically_distinct_diags(G2_HF)
        G2_EHF = mth.topologically_distinct_diags(G2_EHF)
        G2_noHF = mth.topologically_distinct_diags(G2_noHF)
        G3_HF = mth.topologically_distinct_diags(G3_HF)
        G3_EHF = mth.topologically_distinct_diags(G3_EHF)
        G3_noHF = mth.topologically_distinct_diags(G3_noHF)

    G = G2_HF + G2_EHF + G2_noHF + G3_HF + G3_EHF + G3_noHF
    G2 = G2_HF + G2_EHF + G2_noHF
    G3 = G3_HF + G3_EHF + G3_noHF
    nb_2 = len(G2)
    nb_2_HF = len(G2_HF)
    nb_2_EHF = len(G2_EHF)
    nb_2_noHF = len(G2_noHF)
    nb_3 = len(G3)
    nb_3_HF = len(G3_HF)
    nb_3_EHF = len(G3_EHF)
    nb_3_noHF = len(G3_noHF)

numdiag = len(G)
print "Time ellapsed: ", datetime.now() - start_time
print "Number of connected diagrams, ", numdiag

if theory == "BMBPT":
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
nb_time_diags = 0
if theory == "BMBPT" and not norm:
    feynman_expressions = []
    diag_expressions = []
    if write_time:
        G_time = []
        time_indexes = []
    for diag in G:
        # Attribute a qp label to all propagators
        mth.attribute_qp_labels(diag)
        # Check for the time-structure diagram
        time_diag = mth.time_structure_graph(diag)
        if write_time:
            if G_time == []:
                G_time.append(time_diag)
                time_diag_num = G_time.index(time_diag)
            else:
                test = True
                nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
                for good_tdiag in G_time:
                    if nx.is_isomorphic(time_diag, good_tdiag, node_match=nm):
                        test = False
                        time_diag_num = G_time.index(good_tdiag)
                        break
                if test:
                    G_time.append(time_diag)
                    time_diag_num = G_time.index(time_diag)
            time_indexes.append(time_diag_num)
        # Determine the expression depending on the graph structure
        numerator = mth.extract_numerator(diag)
        denominator = ""
        extra_factor = ""
        if nx.is_arborescence(time_diag):
            denominator = mth.time_tree_denominator(
                diag, time_diag, denominator)

        elif (norder == 4) and (mth.number_of_sinks(diag) == 1):
            testdiag = mth.omega_subgraph(diag)
            for i in range(2):
                subgraph_stack = []
                subgraph_stack.append(nx.dag_longest_path(testdiag)[1])
                if i == 0:
                    subgraph_stack.append(nx.dag_longest_path(testdiag)[0])
                else:
                    for vertex_1 in testdiag:
                        test_vertex = True
                        for vertex_2 in nx.dag_longest_path(testdiag):
                            if vertex_1 == vertex_2:
                                test_vertex = False
                        if test_vertex:
                            subgraph_stack.append(vertex_1)
                subdiag = testdiag.subgraph(subgraph_stack)
                denominator += "(" + mth.extract_denom(diag, subdiag) + ")"
            for vertex in diag:
                if diag.out_degree(vertex) == 0:
                    subdiag = diag.subgraph(vertex)
            denominator_a = mth.extract_denom(diag, subdiag)
            denominator_abc = mth.extract_denom(diag, testdiag)
            extra_factor += "\\left[ \\frac{1}{" + denominator_a \
                + "} + \\frac{1}{" + denominator_abc + "} \\right]"
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
                + " }{ " + denominator + " }" + extra_factor + "\n"
        else:
            diag_exp = prefactor + numerator + extra_factor + "\n"
        feynman_expressions.append(feynman_exp)
        diag_expressions.append(diag_exp)
    nb_time_diags = len(G_time)


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
    if write_time:
        for i, time_diag in enumerate(G_time):
            time_name = 'time_%i' % i
            mth.feynmf_generator(time_diag, 'MBPT', time_name)
            shutil.move(time_name + '.tex',
                        directory + "/Diagrams/" + time_name + '.tex')


msg = 'Include diagrams in tex?'
pdiag = raw_input("%s (y/N) " % msg).lower() == 'y'

# Write everything down in a nice LaTeX file
latex_file = open(directory + '/result.tex', 'w')
mth.write_file_header(directory, latex_file, pdiag, norder, theory)

if theory == "BMBPT":
    mth.write_BMBPT_header(latex_file, numdiag, three_N, norm, nb_2_HF,
                           nb_2_EHF, nb_2_noHF, nb_3_HF, nb_3_EHF, nb_3_noHF)
    if write_time and pdiag and pdraw:
        latex_file.write("\\section{Associated time-structure diagrams}\n\n")
        for i in range(nb_time_diags):
            latex_file.write("Time-structure diagram %i:\n" % (i+1))
            latex_file.write('\n\\begin{center}\n')
            time_file = open(directory+"/Diagrams/time_%i.tex" % i)
            latex_file.write(time_file.read())
            latex_file.write('\n\\end{center}\n\n')
            latex_file.write("Related diagrams:")
            for i_diag in range(0, numdiag):
                if time_indexes[i_diag] == i:
                    latex_file.write(" %i," % (i_diag+1))
            latex_file.write("\n\n")
    latex_file.write("\\section{Two-body diagrams}\n\n")
    latex_file.write("\\subsection{Two-body energy canonical diagrams}\n\n")

begeq = "\\begin{equation}\n"
endeq = "\\end{equation}\n"
for i_diag in range(0, numdiag):
    if theory == "BMBPT":
        mth.write_BMBPT_section(latex_file, i_diag, three_N, norm,
                                nb_2, nb_2_HF, nb_2_EHF, nb_3_HF, nb_3_EHF)
        latex_file.write("Diagram %i:\n" % (i_diag + 1))
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
        mth.draw_diagram(directory, latex_file, i_diag, 'diag')
        if write_time:
            i_tdiag = time_indexes[i_diag]
            latex_file.write('\\hspace{10pt} $\\rightarrow$')
            mth.draw_diagram(directory, latex_file, i_tdiag, 'time')
        latex_file.write('\n\\end{center}\n\n')
enddoc = "\\end{document}"
latex_file.write(enddoc)
latex_file.close()

msg = 'Compile pdf?'
pdfcompile = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdfcompile:
    mth.compile_and_clean(directory, pdiag, numdiag, write_time, nb_time_diags)
