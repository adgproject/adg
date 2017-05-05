#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

import copy
import os
import itertools
import multiprocessing
from datetime import datetime
import string
import shutil
import numpy as np
import networkx as nx


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


def seed(n):
    """Generate all 1-magic square of dimension n."""
    return [k for k in itertools.permutations(range(n), n)]


def no_trace(matrices):
    """Select matrices with full 0 diagonal."""
    traceless_matrices = []
    for matrix in matrices:
        test_traceless = True
        for i, n in enumerate(matrix):
            if n[i] == 1:
                test_traceless = False
                break
        if test_traceless:
            traceless_matrices.append(matrix)
    return traceless_matrices


def no_loop(matrices):
    """Select out matrices with loops between two vertices."""
    no_loop_matrices = []
    for matrix in matrices:
        test_no_loop = True
        for i in range(len(matrix[0])):
            for j in range(i+1):
                if (matrix[i][j] != 0) and (matrix[j][i] != 0):
                    test_no_loop = False
                    break
        if test_no_loop:
            no_loop_matrices.append(matrix)
    return no_loop_matrices


def check_degree(matrices, three_N_use):
    """Check the degrees of the vertices
    (i.e. its effective one-, two- or three-body structure).
    """
    deg_ok = []
    for matrix in matrices:
        test_degree = True
        for i in range(len(matrix[0])):
            degree = 0
            for j in range(len(matrix[0])):
                degree += matrix[i][j] + matrix[j][i]
            if (degree != 2) and (degree != 4):
                if (not three_N_use) or (degree != 6):
                    test_degree = False
                    break
        if test_degree:
            deg_ok.append(matrix)
    return deg_ok


def check_vertex_degree(matrices, three_N_use, vertex_id):
    """Check the degree of a specific vertex in a set of matrices"""
    good_matrices = []
    for matrix in matrices:
        vertex_degree_OK = True
        vertex_degree = 0
        for index in range(len(matrix[0])):
            vertex_degree += matrix[index][vertex_id] \
                + matrix[vertex_id][index]
        if (vertex_degree != 2) and (vertex_degree != 4):
            if (not three_N_use) or (vertex_degree != 6):
                vertex_degree_OK = False
        if vertex_degree_OK:
            good_matrices.append(matrix)
    return good_matrices


def empty_matrix_generation(size):
    """Generate an empty matrix of size (size,size)"""
    empty_matrix = []
    for line in range(size):
        empty_matrix.append([])
        for element in range(size):
            empty_matrix[line].append(0)
    return empty_matrix


def diagram_generation(n):
    """Generate the diagrams for the MBPT case."""
    seeds = seed(n)
    all = [[[0 if i != j else 1 for i in range(n)] for j in k] for k in seeds]
    traceless = no_trace(all)
    coeffs = [i for i in itertools.combinations_with_replacement(
        range(len(traceless)), 2)]
    double = []

    for coef in coeffs:
        matrix = copy.deepcopy(traceless[coef[0]])
        for i, line in enumerate(traceless[coef[1]]):
            for j, elem in enumerate(line):
                matrix[i][j] += elem
        double.append(matrix)
    doubleUniq = []
    for matrix in double:
        if matrix not in doubleUniq:
            doubleUniq.append(matrix)
    doubleUniq.sort(reverse=True)
    mbpt_diagrams = []
    for matrix in doubleUniq:
        mbpt_diagrams.append(np.array(matrix))
    return mbpt_diagrams


def BMBPT_generation(p_order, three_N_use, norm_diagrams):
    """Generate diagrams for BMBPT from bottom up."""

    deg_max = 4
    if three_N_use:
        deg_max = 6

    # Create a null oriented adjacency matrix of dimension (p_order,p_order)
    temp_matrices = []
    temp_matrices.append(empty_matrix_generation(p_order))

    # Generate oriented adjacency matrices going vertex-wise
    for vertex in range(p_order):
        for sum_index in range(vertex+1, p_order):
            matrices = []
            for mat in temp_matrices:
                matrices.append(mat)
                if mat[sum_index][vertex] == 0:
                    vert_degree = 0
                    for k in range(0, p_order):
                        vert_degree += mat[k][vertex] + mat[vertex][k]
                    elem = 1
                    while (elem + vert_degree) <= deg_max:
                        temp_mat = copy.deepcopy(mat)
                        temp_mat[vertex][sum_index] = elem
                        matrices.append(temp_mat)
                        elem += 1
            temp_matrices = copy.deepcopy(matrices)
            # Column not iterated upon for first vertex in operator diagrams
            if norm_diagrams or (vertex != 0):
                matrices = []
                for mat in temp_matrices:
                    matrices.append(mat)
                    if mat[vertex][sum_index] == 0:
                        vert_degree = 0
                        for k in range(0, p_order):
                            vert_degree += mat[vertex][k] + mat[k][vertex]
                        elem = 1
                        while (elem + vert_degree) <= deg_max:
                            temp_mat = copy.deepcopy(mat)
                            temp_mat[sum_index][vertex] = elem
                            matrices.append(temp_mat)
                            elem += 1
                temp_matrices = copy.deepcopy(matrices)
        temp_matrices = check_vertex_degree(matrices, three_N_use, vertex)
        matrices = copy.deepcopy(temp_matrices)

    # Checks to exclude non-conform matrices
    matrices = check_degree(matrices, three_N_use)
    matrices = no_loop(matrices)
    matricesUniq = []
    for mat in matrices:
        if mat not in matricesUniq:
            matricesUniq.append(mat)
    matricesUniq.sort(reverse=True)
    bmbpt_diagrams = []
    for mat in matricesUniq:
        bmbpt_diagrams.append(np.array(mat))
    return bmbpt_diagrams


def line_label_h(n):
    """Select appropriate label for hole line."""
    labels = list(string.ascii_lowercase)
    labels = labels[0:15]
    return labels[n]


def line_label_p(n):
    """Select appropriate label for particle line."""
    labels = list(string.ascii_lowercase)
    labels = labels[15:-1]
    return labels[n]


# Start computing everything
print "Running"
start_time = datetime.now()
if theory == "MBPT":
    diagrams = diagram_generation(norder)
elif theory == "BMBPT":
    diagrams = BMBPT_generation(norder, three_N, norm)
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
    G1 = []
    nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
    for diag in G:
        if G1 == []:
            G1.append(diag)
        else:
            test = True
            for good_diag in G1:
                if nx.is_isomorphic(diag, good_diag, node_match=nm):
                    test = False
                    break
            if test:
                G1.append(diag)
    G = G1
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
    for diag in G:
        max_deg = 0
        for node in diag:
            max_deg = max(max_deg, diag.degree(node))
        if max_deg == 6:
            G3.append(diag)
        else:
            G2.append(diag)
    for diag in G2:
        test_HF = True
        test_EHF = True
        for node in diag:
            if diag.degree(node) == 2:
                test_HF = False
                if node != 0:
                    test_EHF = False
        if test_HF:
            G2_HF.append(diag)
        elif (not test_EHF) or norm:
            G2_noHF.append(diag)
        else:
            G2_EHF.append(diag)
    for diag in G3:
        test_HF = True
        test_EHF = True
        for node in diag:
            if diag.degree(node) == 2:
                test_HF = False
                if node != 0:
                    test_EHF = False
        if test_HF:
            G3_HF.append(diag)
        elif (not test_EHF) or norm:
            G3_noHF.append(diag)
        else:
            G3_EHF.append(diag)
    G = G2_HF + G2_EHF + G2_noHF + G3_HF + G3_EHF + G3_noHF
    nb_2 = len(G2)
    nb_2_HF = len(G2_HF)
    nb_2_EHF = len(G2_EHF)
    nb_2_noHF = len(G2_noHF)
    nb_3 = len(G3)
    nb_3_HF = len(G3_HF)
    nb_3_EHF = len(G3_EHF)
    nb_3_noHF = len(G3_noHF)
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
                        bra = bra + line_label_h(col)
                    else:
                        bra = bra + line_label_p(col)
                if (incidence[row, col] == -1):
                    if type_edg[col] == 'h':
                        ket = ket + line_label_h(col)
                    else:
                        ket = ket + line_label_p(col)
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
                        denom = denom + ' +E_' + line_label_h(col)
                    else:
                        denom = denom + ' +E_' + line_label_p(col)
                if val_test == -1:
                    if type_edg[col] == 'h':
                        denom = denom + '-E_' + line_label_h(col)
                    else:
                        denom = denom + '-E_' + line_label_p(col)
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
if theory == "BMBPT":
    feynman_expressions = []
    diag_expressions = []
    for diag in G:
        # Attribute a qp label to all propagators
        i = 1
        nb_down_props = 0
        for prop in diag.edges_iter(keys=True):
            diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] = "k_{%i}" % i
            if prop[1] < prop[0]:
                nb_down_props += 1
            i += 1
        # Determine the numerator corresponding to the diagram
        numerator = ""
        for vertex in nx.nodes(diag):
            # Attribute the correct operator to each vertex
            if diag.node[vertex]['operator']:
                numerator += "O"
            else:
                numerator += "\\Omega"
            # Attribute the good "type number" to each vertex
            numerator = numerator + "^{%i" % diag.out_degree(vertex) \
                + "%i}_{" % diag.in_degree(vertex)
            # First add the qp states corresponding to propagators going out
            for prop in diag.out_edges_iter(vertex, keys=True):
                numerator = numerator \
                    + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state']
            # Add the qp states corresponding to propagators coming in
            for prop in diag.in_edges_iter(vertex, keys=True):
                numerator = numerator  \
                    + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state']
            numerator = numerator + "} "
        # Determine the denominator corresponding to the diagram
        # First determine the type of structure we have
        denominator = ""
        testgraph_stack = []
        # Create a subgraph without the operator vertex
        for vertex in range(1, norder):
            testgraph_stack.append(vertex)
        testdiag = diag.subgraph(testgraph_stack)
        test_adg_subgraphs = True
        # Use the subgraph to determine the structure of the overall graph
        for connected_subgraph in nx.weakly_connected_component_subgraphs(testdiag):
            if len(connected_subgraph) > 1:
                if nx.dag_longest_path_length(connected_subgraph) != (len(connected_subgraph)-1):
                    test_adg_subgraphs = False
        # If the graph has the appropriate structure, determine the denominator
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
                    denominator = denominator + "("
                    for prop in diag.in_edges_iter(subdiag, keys=True):
                        if subdiag.has_edge(prop[0], prop[1], prop[2]) is False:
                            denominator = denominator + " + E_{" + \
                                diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
                    for prop in diag.out_edges_iter(subdiag, keys=True):
                        if subdiag.has_edge(prop[0], prop[1], prop[2]) is False:
                            denominator = denominator + " - E_{" + \
                                diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
                    denominator = denominator + ")"
        # Determine the time structure of the graph
        diag_copy = diag.to_directed()
        for vertex in range(1, len(diag_copy)):
            if diag_copy.in_degree(vertex) == 0:
                diag_copy.add_edge(0, vertex)
        time_diag = nx.DiGraph()
        for vertex_i in diag_copy.nodes_iter():
            for vertex_j in diag_copy.nodes_iter():
                lgst_path = []
                for path in nx.all_simple_paths(diag_copy,
                                                source=vertex_i,
                                                target=vertex_j):
                    if len(path) > len(lgst_path):
                        lgst_path = path
                time_diag.add_path(lgst_path)
        TimeStructureIsTree = nx.is_arborescence(time_diag)
        # Determine the integral component in the Feynman expression
        integral = ""
        for vertex in range(1, norder):
            integral = integral + "\\mathrm{d}\\tau_%i" % vertex
        if norder > 2:
            for vertex_i in range(1, norder):
                for vertex_j in range(1, norder):
                    if diag.has_edge(vertex_i, vertex_j):
                        integral = integral + "\\theta(\\tau_%i" % vertex_j \
                            + "-\\tau_%i) " % vertex_i
        for vertex in range(1, norder):
            integral = integral + "e^{-\\tau_%i (" % vertex
            for prop in diag.in_edges_iter(vertex, keys=True):
                integral = integral + " + E_{" \
                    + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
            for prop in diag.out_edges_iter(vertex, keys=True):
                integral = integral + " - E_{" \
                    + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
            integral = integral + ")}"
        # Determine the pre-factor
        prefactor = "(-1)^%i " % (norder - 1)
        if nb_down_props % 2 != 0:
            prefactor = "-" + prefactor
        sym_fact = ""
        prop_multiplicity = []
        for i in range(6):
            prop_multiplicity.append(0)
        for vertex_i in diag:
            for vertex_j in diag:
                if diag.number_of_edges(vertex_i, vertex_j) >= 2:
                    prop_multiplicity[diag.number_of_edges(
                        vertex_i, vertex_j) - 1] += 1

        for prop_id, multiplicity in enumerate(prop_multiplicity):
            if multiplicity == 1:
                sym_fact = sym_fact + "(%i!)" % (prop_id+1)
            elif multiplicity >= 2:
                sym_fact = sym_fact + "(%i!)" % (prop_id+1) \
                    + "^%i" % multiplicity
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
        feynman_expressions.append(feynman_exp)
        diag_expressions.append(diag_exp)


def feynmf_generator(diag, theory_type, diagram_name):
    """Generate the feynmanmp instructions corresponding to the diagram."""
    p_order = diag.number_of_nodes()
    diag_size = 20*p_order

    theories = ["MBPT", "BMBPT", "SCGF"]
    prop_types = ["half_prop", "prop_pm", "double_arrow"]
    prop = prop_types[theories.index(theory_type)]

    fmf_file = open(diagram_name + ".tex", 'w')
    begin_file = "\\parbox{%i" % diag_size + "pt}{\\begin{fmffile}{" \
        + diagram_name + "}\n\\begin{fmfgraph*}(%i" % diag_size \
        + ",%i)\n" % diag_size
    fmf_file.write(begin_file)

    # Set the position of the vertices
    fmf_file.write("\\fmftop{v%i}\\fmfbottom{v0}\n" % (p_order-1))
    for vert in range(p_order-1):
        fmf_file.write("\\fmf{phantom}{v%i" % vert + ",v%i}\n" % (vert+1))
        if diag.node[vert]['operator']:
            fmf_file.write("\\fmfv{d.shape=square,d.filled=full,d.size=3thick")
        else:
            fmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick")
        fmf_file.write("}{v%i}\n" % vert)
    fmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}")
    fmf_file.write("{v%i}\n" % (p_order-1))
    fmf_file.write("\\fmffreeze\n")

    # Loop over all elements of the matrix to draw associated propagators
    for vert_i in range(0, p_order):
        for vert_j in range(0, p_order):
            props_left_to_draw = diag.number_of_edges(vert_i, vert_j)
            # Special config for consecutive vertices
            if (props_left_to_draw % 2 == 1) and (abs(vert_i-vert_j) == 1):
                fmf_file.write("\\fmf{" + prop)
                # Check for specific MBPT configuration
                if diag.number_of_edges(vert_j, vert_i) == 1:
                    fmf_file.write(",right=0.5")
                fmf_file.write("}{v%i," % vert_i + "v%i}\n" % vert_j)
                props_left_to_draw -= 1
            while props_left_to_draw > 0:
                fmf_file.write("\\fmf{" + prop + ",")
                if props_left_to_draw % 2 == 1:
                    fmf_file.write("right=")
                else:
                    fmf_file.write("left=")
                if (props_left_to_draw == 6) or (props_left_to_draw == 5):
                    fmf_file.write("0.9")
                elif (props_left_to_draw == 4) or (props_left_to_draw == 3):
                    fmf_file.write("0.75")
                elif (props_left_to_draw == 2) or (props_left_to_draw == 1):
                    if abs(vert_i-vert_j) == 1:
                        fmf_file.write("0.5")
                    else:
                        fmf_file.write("0.6")
                fmf_file.write("}{v%i," % vert_i + "v%i}\n" % vert_j)
                props_left_to_draw -= 1
    fmf_file.write("\\end{fmfgraph*}\n\\end{fmffile}}\n\n")
    fmf_file.close()


# Writing a feynmp file for each graph
msg = 'Generate diagrams feymanmf instructions?'
pdraw = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdraw:
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    for i in range(0, numdiag):
        diag_name = 'diag_%i' % i
        feynmf_generator(G[i], theory, diag_name)
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
        diag_exp = diag_expressions[i_diag]
        feynman_exp = feynman_expressions[i_diag]
        latex_file.write("Diagram %i:\n" % (i_diag+1))
        latex_file.write(begeq)
        latex_file.write(feynman_exp)
        latex_file.write(endeq)
    elif theory == "MBPT":
        diag_exp = "\\dfrac{1}{%i}" % nedges_eq[i_diag] + phases[i_diag] \
            + "\\sum{\\dfrac{" + mat_els[i_diag] + "}{" \
            + denoms[i_diag] + "}}\n"
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
