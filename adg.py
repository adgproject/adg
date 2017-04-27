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
norder = int(raw_input('Order of the diagrams ?\n'))
theory = raw_input('MBPT or BMBPT ?\n').upper()

three_N = False
norm = False
if theory == "BMBPT":
    three_N = raw_input("Include three-body forces ? (y/N) ").lower() == 'y'
    norm = raw_input("Compute norm kernel instead of operator kernels ? (y/N) ").lower() == 'y'

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
    traceless = []
    for matrix in matrices:
        test = True
        for i, n in enumerate(matrix):
            if n[i] == 1:
                test = False
                break
        if test:
            traceless.append(matrix)
    return traceless


def no_loop(matrices):
    """Select out matrices with loops between two vertices."""
    no_loop = []
    for matrix in matrices:
        test = True
        for i in range(len(matrix[0])):
            for j in range(i+1):
                if (matrix[i][j] != 0) and (matrix[j][i] != 0):
                    test = False
                    break
        if test:
            no_loop.append(matrix)
    return no_loop


def check_degree(matrices, three_N):
    """Check the degrees of the vertices
    (i.e. its effective one-, two- or three-body structure).
    """
    deg_ok = []
    for matrix in matrices:
        test = True
        for i in range(len(matrix[0])):
            degree = 0
            for j in range(len(matrix[0])):
                degree += matrix[i][j] + matrix[j][i]
            if (degree != 2) and (degree != 4):
                if (not three_N) or (degree != 6):
                    test = False
                    break
        if test:
            deg_ok.append(matrix)
    return deg_ok


def diagram_generation(n):
    """Generate the diagrams for the MBPT case."""
    seeds = seed(n)
    all = [[[0 if i != j else 1 for i in range(n)] for j in k] for k in seeds]
    traceless = no_trace(all)
    coeffs = [i for i in itertools.combinations_with_replacement(range(len(traceless)), 2)]
    double = []
    for coef in coeffs:
        matrix = copy.deepcopy(traceless[coef[0]])
        for i, line in enumerate(traceless[coef[1]]):
            for j, elem in enumerate(line):
                matrix[i][j] += elem
        double.append(matrix)
    doubleUniq = []
    for i in double:
        if i not in doubleUniq:
            doubleUniq.append(i)
    doubleUniq.sort(reverse=True)
    diagrams = []
    for el in doubleUniq:
        diagrams.append(np.array(el))
    return diagrams


def BMBPT_generation(p_order, three_N, norm):
    """Generate diagrams for BMBPT from bottom up."""
    # Begin by creating a zero oriented adjacency matric of good dimensions
    empty_mat = []
    for i in range(p_order):
        empty_mat.append([])
        for j in range(p_order):
            empty_mat[i].append(0)

    deg_max = 4
    if three_N:
        deg_max = 6

    temp_matrices = []
    temp_matrices.append(empty_mat)

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
            # Column is not iterated upon for the first vertex in operator diagrams
            if norm or (vertex != 0):
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
        deg_vertex_ok = []
        for matrix in temp_matrices:
            test = True
            degree = 0
            for i in range(0, p_order):
                degree += matrix[i][vertex] + matrix[vertex][i]
            if (degree != 2) and (degree != 4):
                if (not three_N) or (degree != 6):
                    test = False
            if test:
                deg_vertex_ok.append(matrix)
        matrices = copy.deepcopy(deg_vertex_ok)
        temp_matrices = copy.deepcopy(deg_vertex_ok)

    # Checks to exclude non-conform matrices
    good_degree = check_degree(matrices, three_N)
    mat_wo_loops = no_loop(good_degree)
    matricesUniq = []
    for i in mat_wo_loops:
        if i not in matricesUniq:
            matricesUniq.append(i)
    matricesUniq.sort(reverse=True)
    diagrams = []
    for el in matricesUniq:
        diagrams.append(np.array(el))
    return diagrams


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

### Graph part (computing, writing, drawing)
G = []
for diagram in diagrams:
    G.append(nx.from_numpy_matrix(diagram, create_using=nx.MultiDiGraph(), parallel_edges=True))
G1 = []
for diag in G:
    if (nx.number_weakly_connected_components(diag)) == 1:
        G1.append(diag)
G = G1
# Specific checks for loop diagrams and topologically identical diagrams in BMBPT
if theory == "BMBPT":
    G1 = []
    for diag in G:
        if nx.is_directed_acyclic_graph(diag):
            G1.append(diag)
    G = G1
    G1 = []
    for diag in G:
        test = True
        if not norm:
            # Account for different status of vertices in operator diagrams
            for node in diag:
                diag.node[node]['operator'] = False
            diag.node[0]['operator'] = True
            nm = nx.algorithms.isomorphism.categorical_node_match('operator', False)
        if G1 == []:
            G1.append(diag)
        else:
            for good_diag in G1:
                if norm:
                    if nx.is_isomorphic(diag, good_diag):
                        test = False
                        break
                else:
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

### Algebraic expressions:
### CAVEAT !!! This works only for MBPT


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


# Treatment of the algebraic expressions
### To be extended to BMBPT in the future
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
        # neq_lines=np.asarray(list(i for i in set(map(tuple, eq_lines)))).transpose()
        neq_lines = np.asarray(list(i for i in set(map(tuple, eq_lines))))
        n_sym = len(eq_lines)-len(neq_lines)
        #### CAVEAT !!! Valid only for *MBPT*
        nedges_eq.append(2**n_sym)
        # print "After neqlines"
        #### Loops

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
                numerator += "\Omega"
            # Attribute the good "type number" to each vertex
            numerator = numerator + "^{%i" % diag.out_degree(vertex) + "%i}_{" % diag.in_degree(vertex)
            # First add the qp states corresponding to propagators going out
            for prop in diag.out_edges_iter(vertex, keys=True):
                numerator = numerator + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state']
            # Add the qp states corresponding to propagators coming in
            for prop in diag.in_edges_iter(vertex, keys=True):
                numerator = numerator + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state']
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
                            denominator = denominator + " + E_{" + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
                    for prop in diag.out_edges_iter(subdiag, keys=True):
                        if subdiag.has_edge(prop[0], prop[1], prop[2]) is False:
                            denominator = denominator + " - E_{" + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
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
                for path in nx.all_simple_paths(diag_copy, source=vertex_i, target=vertex_j):
                    if len(path) > len(lgst_path):
                        lgst_path = path
                time_diag.add_path(lgst_path)
        TimeStructureIsTree = nx.is_arborescence(time_diag)
        # Determine the integral component in the Feynman expression
        integral = ""
        for vertex in range(1, norder):
            integral = integral + "\mathrm{d}\\tau_%i" % vertex
        if norder > 2:
            for vertex_i in range(1, norder):
                for vertex_j in range(1, norder):
                    if diag.has_edge(vertex_i, vertex_j):
                        integral = integral + "\\theta(\\tau_%i" % vertex_j + "-\\tau_%i) " % vertex_i
        for vertex in range(1, norder):
            integral = integral + "e^{-\\tau_%i (" % vertex
            for prop in diag.in_edges_iter(vertex, keys=True):
                integral = integral + " + E_{" + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
            for prop in diag.out_edges_iter(vertex, keys=True):
                integral = integral + " - E_{" + diag.edge[prop[0]][prop[1]][prop[2]]['qp_state'] + "}"
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
                    prop_multiplicity[diag.number_of_edges(vertex_i, vertex_j)-1] += 1
        for prop_id, multiplicity in enumerate(prop_multiplicity):
            if multiplicity == 1:
                sym_fact = sym_fact + "(%i!)" % (prop_id+1)
            elif multiplicity >= 2:
                sym_fact = sym_fact + "(%i!)" % (prop_id+1) + "^%i" % multiplicity
        if sym_fact != "":
            prefactor = "\\frac{" + prefactor + "}{" + sym_fact + "}\sum_{k_i}"
        else:
            prefactor = prefactor + "\sum_{k_i}"
        feynman_exp = prefactor + numerator + "\int_{0}^{\\tau}" + integral + "\n"
        if denominator != "":
            diag_exp = prefactor + "\\frac{ " + numerator + " }{ " + denominator + " }\n"
        else:
            diag_exp = prefactor + numerator + "\n"
        feynman_expressions.append(feynman_exp)
        diag_expressions.append(diag_exp)


def feynmf_generator(diag, theory, diag_name):
    """Generate the feynmanmp instructions corresponding to the diagram."""
    p_order = diag.number_of_nodes()
    diag_size = 20*p_order

    theories = ["MBPT", "BMBPT", "SCGF"]
    th_index = theories.index(theory)
    prop_types = ["half_prop", "prop_pm", "double_arrow"]
    prop = prop_types[th_index]

    file_name = diag_name + ".tex"

    feynmf_file = open(file_name, 'w')

    begin_file = "\parbox{%i" % diag_size + "pt}{\\begin{fmffile}{" + diag_name + "}\n\\begin{fmfgraph*}(%i" % diag_size + ",%i)\n" % diag_size
    end_file = "\end{fmfgraph*}\n\end{fmffile}}\n\n"

    # Draw the diagram only if there is one to be drawn...
    if p_order >= 2:
        feynmf_file.write(begin_file)

        # Set the position of the top and bottom vertices
        feynmf_file.write("\\fmftop{v%i}\\fmfbottom{v0}\n" % (p_order-1))
        if (theory == "BMBPT") and (not norm):
            feynmf_file.write("\\fmfv{d.shape=square,d.filled=full,d.size=3thick}{v0}\n")
        else:
            feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v0}\n")
        feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" % (p_order-1))
        # Set the position of intermediate vertices if needed
        if p_order > 2:
            feynmf_file.write("\\fmf{phantom}{v0,v1}\n")
            for vertex in range(1, p_order-2):
                feynmf_file.write("\\fmf{phantom}{v%i" % vertex + ",v%i}\n" % (vertex+1))
                feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" % vertex)
            feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" % (p_order-2))
            feynmf_file.write("\\fmf{phantom}{v%i," % (p_order-2) + "v%i}\n" % (p_order-1))
            feynmf_file.write("\\fmffreeze\n")

        # Recover the oriented adjacency matrix of the diagram
        oriented_adj_mat = []
        for i in range(0, p_order):
            oriented_adj_mat.append([])
            for j in range(0, p_order):
                oriented_adj_mat[i].append(0)
        for line in nx.generate_edgelist(diag, data=False):
            i = int(line[0])
            j = int(line[2])
            oriented_adj_mat[i][j] += 1

        # Loop over all elements of the matrix to draw associated propagators
        for i in range(0, p_order):
            for j in range(0, p_order):
                # For directly consecutive vertices
                if (abs(i-j) == 1) and (oriented_adj_mat[i][j] != 0):
                    if oriented_adj_mat[i][j] == 1:
                        if oriented_adj_mat[j][i] != 1:
                            feynmf_file.write("\\fmf{" + prop + "}{v%i," % i + "v%i}\n" % j)
                        else:
                            feynmf_file.write("\\fmf{" + prop + ",right=0.5}{v%i," % i + "v%i}\n" % j)
                    else:
                        feynmf_file.write("\\fmf{" + prop + ",right=0.5}{v%i," % i + "v%i}\n" % j)
                        feynmf_file.write("\\fmf{" + prop + ",left=0.5}{v%i," % i + "v%i}\n" % j)
                        if oriented_adj_mat[i][j] == 3:
                            feynmf_file.write("\\fmf{" + prop + "}{v%i," % i + "v%i}\n" % j)
                        elif oriented_adj_mat[i][j] >= 4:
                            feynmf_file.write("\\fmf{" + prop + ",right=0.75}{v%i," % i + "v%i}\n" % j)
                            feynmf_file.write("\\fmf{" + prop + ",left=0.75}{v%i," % i + "v%i}\n" % j)
                            if oriented_adj_mat[i][j] >= 5:
                                feynmf_file.write("\\fmf{" + prop + ",right=0.9}{v%i," % i + "v%i}\n" % j)
                                if oriented_adj_mat[i][j] == 6:
                                    feynmf_file.write("\\fmf{" + prop + ",left=0.9}{v%i," % i + "v%i}\n" % j)

                # For more distant vertices
                elif (i != j) and (oriented_adj_mat[i][j] != 0):
                    if (oriented_adj_mat[i][j] == 1) and (oriented_adj_mat[j][i] == 2):
                        feynmf_file.write("\\fmf{" + prop + ",right=0.75}{v%i," % i + "v%i}\n" % j)
                    else:
                        feynmf_file.write("\\fmf{" + prop + ",right=0.6}{v%i," % i + "v%i}\n" % j)
                        if oriented_adj_mat[i][j] != 1:
                            feynmf_file.write("\\fmf{" + prop + ",left=0.6}{v%i," % i + "v%i}\n" % j)
                            if oriented_adj_mat[i][j] != 2:
                                feynmf_file.write("\\fmf{" + prop + ",right=0.75}{v%i," % i + "v%i}\n" % j)
                                if oriented_adj_mat[i][j] != 3:
                                    feynmf_file.write("\\fmf{" + prop + ",left=0.75}{v%i," % i + "v%i}\n" % j)
                                    if oriented_adj_mat[i][j] != 4:
                                        feynmf_file.write("\\fmf{" + prop + ",right=0.9}{v%i," % i + "v%i}\n" % j)
        feynmf_file.write(end_file)
    else:
        print "Perturbative order too small"


# Writing a feynmp file for each graph
msg = 'Generate diagrams feymanmf instructions ?'
pdraw = raw_input("%s (y/N) " % msg).lower() == 'y'
if pdraw:
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    for i in range(0, numdiag):
        diag_name = 'diag_%i' % i
        feynmf_generator(G[i], theory, diag_name)
        shutil.move(diag_name + '.tex', directory + "/Diagrams/" + diag_name + '.tex')


msg = 'Include diagrams in tex ?'
pdiag = raw_input("%s (y/N) " % msg).lower() == 'y'
### Write everything down in a nice LaTeX file
header = "\documentclass[10pt,a4paper]{article}\n \usepackage[utf8]{inputenc}\n\usepackage{braket}\n\usepackage{graphicx}\n"
header = header + "\usepackage[english]{babel}\n\usepackage{amsmath}\n\usepackage{amsfonts}\n\usepackage{amssymb}\n"
if pdiag:
    header = header + "\usepackage[force]{feynmp-auto}\n"
land = False
if norder > 3:
    msg = 'Expressions may be long rotate pdf ?'
    land = raw_input("%s (y/N) " % msg).lower() == 'y'
if land:
    header = header + "\usepackage[landscape]{geometry}\n"

header = header + "\\title{Diagrams and algebraic expressions at order %i" % norder + " in " + theory + "}\n"
header = header + "\\author{RDL, JR, PA, MD, AT}\n"
latex_file = open(directory + '/result.tex', 'w')
latex_file.write(header)
begdoc = "\\begin{document}\n"
enddoc = "\\end{document}"
begeq = "\\begin{equation}\n"
endeq = "\\end{equation}\n"
latex_file.write(begdoc)
latex_file.write("\maketitle\n")
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

if not pdiag or not pdraw:
    for i_diag in range(0, numdiag):
        if theory == "MBPT":
            diag_exp = "\dfrac{1}{%i}" % nedges_eq[i_diag]+phases[i_diag]+"\sum{\dfrac{"+mat_els[i_diag]+"}{"+denoms[i_diag]+"}}\n"
        elif theory == "BMBPT":
            diag_exp = diag_expressions[i_diag]
            feynman_exp = feynman_expressions[i_diag]
            latex_file.write("Diagram %i:\n" % (i_diag+1))
            latex_file.write(begeq)
            latex_file.write(feynman_exp)
            latex_file.write(endeq)
        latex_file.write(begeq)
        latex_file.write(diag_exp)
        latex_file.write(endeq)
    latex_file.write(enddoc)
else:
    if theory == "BMBPT":
        latex_file.write("\section{Two-body diagrams}\subsection{Two-body energy canonical diagrams}\n")
    for i_diag in range(0, numdiag):
        if theory == "BMBPT":
            if (i_diag == nb_2_HF) and (not norm):
                latex_file.write("\subsection{Two-body canonical diagrams for a generic operator only}\n")
            elif i_diag == nb_2_HF + nb_2_EHF:
                latex_file.write("\subsection{Two-body non-canonical diagrams}\n")
            if three_N:
                if i_diag == nb_2:
                    latex_file.write("\section{Three-body diagrams}\n\subsection{Three-body energy canonical diagrams}\n")
                elif (i_diag == nb_2 + nb_3_HF) and (not norm):
                    latex_file.write("\subsection{Three-body canonical diagrams for a generic operator only}\n")
                elif i_diag == nb_2 + nb_3_HF + nb_3_EHF:
                    latex_file.write("\subsection{Three-body non-canonical diagrams}\n")
        latex_file.write("Diagram %i:\n" % (i_diag+1))
        if theory == "MBPT":
            diag_exp = "\dfrac{1}{%i}" % nedges_eq[i_diag] + phases[i_diag] + "\sum{\dfrac{" + mat_els[i_diag] + "}{" + denoms[i_diag] + "}}\n"
        elif theory == "BMBPT":
            diag_exp = diag_expressions[i_diag]
            feynman_exp = feynman_expressions[i_diag]
            latex_file.write(begeq)
            latex_file.write(feynman_exp)
            latex_file.write(endeq)
        latex_file.write(begeq)
        latex_file.write(diag_exp)
        latex_file.write(endeq)
        latex_file.write('\n\\begin{center}\n')
        diag_file = open(directory+"/Diagrams/diag_%i.tex" % i_diag)
        latex_file.write(diag_file.read())
        latex_file.write('\\end{center}\n\n')
    latex_file.write(enddoc)
latex_file.close()

msg = 'Compile pdf ?'
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
