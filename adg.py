#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

import os
import multiprocessing
from datetime import datetime
import shutil
import cProfile
import pstats
import StringIO
import numpy as np
import networkx as nx
from meliae import scanner
import methods as mth
import bmbpt
import mbpt
import time_structure as tst


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
    print "Parallel Mode available"
    num_cores = multiprocessing.cpu_count()
    print "There is %i" % num_cores + " core(s) available"
    use_parallel = raw_input("Use parallel processing? (y/N)").lower() == 'y'
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

pr = cProfile.Profile()
pr.enable()

# Start computing everything
print "Running"
start_time = datetime.now()
if theory == "MBPT":
    diagrams = mbpt.diagram_generation(norder)
elif theory == "BMBPT":
    diagrams = bmbpt.BMBPT_generation(norder, three_N, norm)
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

for i_diag in xrange(len(G)-1, -1, -1):
    if (nx.number_weakly_connected_components(G[i_diag])) != 1:
        del G[i_diag]

# Specific check for loop diagrams in BMBPT
if theory == "BMBPT":
    for i_diag in xrange(len(G)-1, -1, -1):
        if not nx.is_directed_acyclic_graph(G[i_diag]):
            del G[i_diag]

mth.label_vertices(G, theory, norm)

diagrams = [bmbpt.BmbptFeynmanDiagram(diag, norm, i) for i, diag in enumerate(G)]
print len(diagrams)

# Ordering the diagrams in a convenient way and checking them for doubles
if theory == "BMBPT":
    G2_HF = []
    G2_EHF = []
    G2_noHF = []
    G3_HF = []
    G3_EHF = []
    G3_noHF = []

    # bmbpt.order_2B_or_3B(G, G2, G3)
    G2 = [diag for diag in G if sorted(diag.degree().values())[-1] < 5]
    G3 = [diag for diag in G if sorted(diag.degree().values())[-1] == 6]
    bmbpt.order_HF_or_not(G2, G2_HF, G2_EHF, G2_noHF, norm)
    bmbpt.order_HF_or_not(G3, G3_HF, G3_EHF, G3_noHF, norm)

    diagrams2 = [diag for diag in diagrams if diag.two_or_three_body == 2]
    diagrams3 = [diag for diag in diagrams if diag.two_or_three_body == 3]
    print len(diagrams2)
    print len(diagrams3)
    diagrams2HF = [diag for diag in diagrams2 if diag.HF_type == "HF"]
    diagrams2EHF = [diag for diag in diagrams2 if diag.HF_type == "EHF"]
    diagrams2noHF = [diag for diag in diagrams2 if diag.HF_type == "noHF"]
    diagrams3HF = [diag for diag in diagrams3 if diag.HF_type == "HF"]
    diagrams3EHF = [diag for diag in diagrams3 if diag.HF_type == "EHF"]
    diagrams3noHF = [diag for diag in diagrams3 if diag.HF_type == "noHF"]

    if use_parallel:
        nb_procs_max = 6 if three_N else 3
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
        diagrams2HF = mth.topologically_distinct_diagrams(diagrams2HF)
        diagrams2EHF = mth.topologically_distinct_diagrams(diagrams2EHF)
        diagrams2noHF = mth.topologically_distinct_diagrams(diagrams2noHF)
        diagrams3HF = mth.topologically_distinct_diagrams(diagrams3HF)
        diagrams3EHF = mth.topologically_distinct_diagrams(diagrams3EHF)
        diagrams3noHF = mth.topologically_distinct_diagrams(diagrams3noHF)

    G = G2_HF + G2_EHF + G2_noHF + G3_HF + G3_EHF + G3_noHF
    diagrams = diagrams2HF + diagrams2EHF + diagrams2noHF + diagrams3HF \
        + diagrams3EHF + diagrams3noHF
    nb_2_HF = len(G2_HF)
    nb_2_EHF = len(G2_EHF)
    nb_2_noHF = len(G2_noHF)
    nb_2 = nb_2_HF + nb_2_EHF + nb_2_noHF
    nb_3_HF = len(G3_HF)
    nb_3_EHF = len(G3_EHF)
    nb_3_noHF = len(G3_noHF)
    nb_3 = nb_3_HF + nb_3_EHF + nb_3_noHF
    print
    print len(diagrams)
    print

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
                        bra = bra + mbpt.line_label_h(col)
                    else:
                        bra = bra + mbpt.line_label_p(col)
                if (incidence[row, col] == -1):
                    if type_edg[col] == 'h':
                        ket = ket + mbpt.line_label_h(col)
                    else:
                        ket = ket + mbpt.line_label_p(col)
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
                        denom = denom + ' +E_' + mbpt.line_label_h(col)
                    else:
                        denom = denom + ' +E_' + mbpt.line_label_p(col)
                if val_test == -1:
                    if type_edg[col] == 'h':
                        denom = denom + '-E_' + mbpt.line_label_h(col)
                    else:
                        denom = denom + '-E_' + mbpt.line_label_p(col)
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
        diagrams_time = []
        time_indexes = []
        nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
    for diag in G:
        # Attribute a qp label to all propagators
        bmbpt.attribute_qp_labels(diag)
        # Check for the time-structure diagram
        time_diag = tst.time_structure_graph(diag)
        if write_time:
            if G_time == []:
                G_time.append(time_diag)
                time_diag_num = G_time.index(time_diag)
            else:
                test = True
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
        numerator = bmbpt.extract_numerator(diag)
        denominator = ""
        extra_factor = ""
        if nx.is_arborescence(time_diag):
            denominator = bmbpt.time_tree_denominator(
                diag, time_diag, denominator)

        elif (norder == 4) and (mth.number_of_sinks(diag) == 1):
            testdiag = bmbpt.omega_subgraph(diag)
            for i in range(2):
                subgraph_stack = [nx.dag_longest_path(testdiag)[1]]
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
                denominator += "(" + bmbpt.extract_denom(diag, subdiag) + ")"
            for vertex in diag:
                if diag.out_degree(vertex) == 0:
                    subdiag = diag.subgraph(vertex)
            denominator_a = bmbpt.extract_denom(diag, subdiag)
            denominator_abc = bmbpt.extract_denom(diag, testdiag)
            extra_factor += "\\left[ \\frac{1}{" + denominator_a \
                + "} + \\frac{1}{" + denominator_abc + "} \\right]"
        # Determine the integral component in the Feynman expression
        integral = bmbpt.extract_integral(diag)
        # Determine the pre-factor
        prefactor = "(-1)^%i " % (norder - 1)
        if bmbpt.extract_BMBPT_crossing_sign(diag):
            prefactor = "-" + prefactor
        sym_fact = bmbpt.vertex_exchange_sym_factor(diag) \
            + bmbpt.multiplicity_symmetry_factor(diag)
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
    # scanner.dump_all_objects("memory.dat")
    nb_time_diags = len(G_time)
    diagrams_time = [tst.TimeStructureDiagram(diagram, diagram.tags[0])
                     for diagram in diagrams]
    diagrams_time = mth.topologically_distinct_diagrams(diagrams_time)
    for diag in diagrams:
        for t_diag in diagrams_time:
            if diag.tags[0] in t_diag.tags:
                diag.time_tag = t_diag.tags[0]
                diag.tst_is_tree = True
                break
        diag.attribute_expressions(diagrams_time)
    print len(diagrams_time)
    # print [tst_diag.expr for tst_diag in diagrams_time]

pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
ps.dump_stats("stats.dat")
scanner.dump_all_objects("memory.dat")

# Writing a feynmp file for each graph
msg = 'Generate diagrams feymanmf instructions?'
pdraw = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdraw:
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    mth.create_feynmanmp_files(G, theory, directory, 'diag')
    if write_time:
        mth.create_feynmanmp_files(G_time, theory, directory, 'time')

msg = 'Include diagrams in tex?'
pdiag = raw_input("%s (y/N) " % msg).lower() == 'y'

# Write everything down in a nice LaTeX file
latex_file = open(directory + '/result.tex', 'w')
mth.write_file_header(directory, latex_file, pdiag, norder, theory)

if theory == "BMBPT":
    bmbpt.write_BMBPT_header(latex_file, numdiag, three_N, norm, nb_2_HF,
                             nb_2_EHF, nb_2_noHF, nb_3_HF, nb_3_EHF, nb_3_noHF)
    if write_time:
        latex_file.write("\\section{Associated time-structure diagrams}\n\n")
        time_diag_exps = {}
        for i in xrange(nb_time_diags):
            latex_file.write("\\paragraph{Time-structure diagram T%i:}\n"
                             % (i+1))
            if pdiag and pdraw:
                latex_file.write('\n\\begin{center}\n')
                time_file = open(directory + "/Diagrams/time_%i.tex" % i)
                latex_file.write(time_file.read())
                latex_file.write('\n\\end{center}\n\n')
            if nx.is_arborescence(G_time[i]):
                latex_file.write("Tree: Yes\n\n")
                latex_file.write("\\begin{equation}\n")
                time_diag_exps.setdefault(i, "\\frac{1}{"
                                          + tst.tree_time_structure_den(G_time[i])
                                          + "}\n")
                latex_file.write(time_diag_exps[i])
                latex_file.write("\\end{equation}\n")
            else:
                latex_file.write("Tree: No\n\n")
            latex_file.write("Related Feynman diagrams:")
            for i_diag in range(0, numdiag):
                if time_indexes[i_diag] == i:
                    latex_file.write(" %i," % (i_diag+1))
            latex_file.write("\n\n")
    latex_file.write("\\section{Two-body diagrams}\n\n")
    latex_file.write("\\subsection{Two-body energy canonical diagrams}\n\n")

for i_diag in range(0, numdiag):
    if theory == "BMBPT":
        bmbpt.write_BMBPT_section(latex_file, i_diag, three_N, norm,
                                  nb_2, nb_2_HF, nb_2_EHF, nb_3_HF, nb_3_EHF)
        latex_file.write("\\paragraph{Diagram %i:}\n" % (i_diag + 1))
        if not norm:
            bmbpt.write_diag_exps(latex_file, i_diag, norder,
                                  feynman_expressions[i_diag],
                                  diag_expressions[i_diag])
    elif theory == "MBPT":
        mbpt.write_diag_exp(latex_file, nedges_eq[i_diag], phases[i_diag],
                            mat_els[i_diag], denoms[i_diag])
    if pdiag and pdraw:
        latex_file.write('\n\\begin{center}\n')
        mth.draw_diagram(directory, latex_file, i_diag, 'diag')
        if write_time:
            i_tdiag = time_indexes[i_diag]
            latex_file.write('\\hspace{10pt} $\\rightarrow$ \\hspace{10pt} T%i:'
                             % (i_tdiag + 1))
            mth.draw_diagram(directory, latex_file, i_tdiag, 'time')
        latex_file.write('\n\\end{center}\n\n')
    if theory == 'BMBPT' and write_time:
        i_tdiag = time_indexes[i_diag]
        latex_file.write("\\begin{equation}\n\\text{T}%i = " % (i_tdiag + 1)
                         + time_diag_exps.get(i_tdiag, '')
                         + "\\end{equation}\n")
        bmbpt.write_vertices_values(latex_file, G[i_diag])
enddoc = "\\end{document}"
latex_file.write(enddoc)
latex_file.close()

msg = 'Compile pdf?'
pdfcompile = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdfcompile:
    mth.compile_and_clean(directory, pdiag, numdiag, write_time, nb_time_diags)
