#!/bin/bash /usr/bin/python
# -*- coding: utf-8 -*-

import copy
import os
import numpy as np
import random
from collections import deque
import matplotlib.pyplot as plt
import networkx as nx
import random
import sys
import fileinput
import itertools
#from joblib import Parallel, delayed
import multiprocessing
from datetime import datetime
import string
import shutil


print "#####################"
print "# Automatic Diagram #"
print "#     Generator     #"
print "#    RDL,JR,PA,MD   #"
print "#####################"

print "Parallel Mode"
num_cores = multiprocessing.cpu_count()
print "There is %i" % num_cores + " core(s) available"
norder = int(raw_input('Order of the diagrams ?\n'))

directory = 'Order-%i'% norder
if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(directory+"/Diagrams"):
    os.makedirs(directory+"/Diagrams")


#Generate all 1-magic square of dimension n
def seed(n):
    return [k for k in itertools.permutations(range(n),n)]

#Select matrices with full 0 diagonal
def g(matrices):
    traceless = []
    for matrix in matrices:
        test = True
        for i,n in enumerate(matrix):
            if n[i] == 1:
                test = False
                break
        if test:
            traceless.append(matrix)
    return traceless

#Select matrices with full 0 first line
def v0_specs(matrices):
    v0_specs_ok = []
    for matrix in matrices:
        test = True
        line = matrix[0]
        for i in range(len(line)):
            if line[i] != 0:
                test = False
                break
        if test:
            v0_specs_ok.append(matrix)
    return v0_specs_ok

def diagram_generation(n):
    seeds = seed(n)
    all = [[[0 if i != j else 1 for i in range(n)] for j in k] for k in seeds]
    traceless = g(all)
    coeffs = [i for i in itertools.combinations_with_replacement(range(len(traceless)),2)]
    double = []
    for coef in coeffs:
        matrix = copy.deepcopy(traceless[coef[0]])
        for i,line in enumerate(traceless[coef[1]]):
            for j,elem in enumerate(line):
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

print "Running"
start_time = datetime.now()
diagrams = diagram_generation(norder)
numdiag = len(diagrams)
print "Number of possible diagrams, ",numdiag

i = 0
with open(directory+"/Diagrams.list", "w") as f:
    for diagram in diagrams:
        f.write("Diagram n: %i" % i)
        np.savetxt(f,diagram)
        #diagram.tofile(f,"")
        f.write("\n")
        i += 1

### Graph part (computing, writing, drawing)
G=[]
for diagram in diagrams:
    G.append(nx.from_numpy_matrix(diagram,create_using=nx.MultiDiGraph(),parallel_edges=True))
G1=[]
for diag in G:
    #if (not nx.is_strongly_connected(diag)):
    if((nx.number_strongly_connected_components(diag)) == 1):
        G1.append(diag)
G=G1
numdiag = len(G)
print "Time ellapsed: ",datetime.now() - start_time
print "Number of connected diagrams, ",numdiag


### Algebraic expressions:
### CAVEAT !!! This works only for MBPT

def line_label_h(n):
    labels=list(string.ascii_lowercase)
    labels=labels[0:15]
    return labels[n]
def line_label_p(n):
    labels=list(string.ascii_lowercase)
    labels=labels[15:-1]
    return labels[n]

def mat_elements(irow):
    return


###

mat_els = []
denoms = []
phases = []
nedges_eq = []
for diag in G:
    type_edg =[]
    braket = ''
    #Beware of the sign convention !!!
    incidence = - nx.incidence_matrix(diag,oriented=True).todense()
    nrow = diag.number_of_nodes()
    ncol = diag.number_of_edges()
    n_holes= 0
    diffcols = set()
    for col in range(ncol):
        flat = list(incidence[:,col].A1)
        if(flat.index(1) < flat.index(-1)):
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
            if (incidence[row,col] == 1):
                if (type_edg[col] == 'h'):
                    bra = bra + line_label_h(col)
                else:
                    bra = bra + line_label_p(col)
            if (incidence[row,col] == -1):
                if (type_edg[col] == 'h'):
                    ket = ket + line_label_h(col)
                else:
                    ket = ket + line_label_p(col)
        ###################################
        braket = braket + '\\braket{'+bra+'|H|'+ket+'}'
    mat_els.append(braket)
    denom = ''
    for row in range(1,nrow):
        denom = denom + '('
        for col in range(ncol):
            val_test = incidence[0:row,col].sum()
            if (val_test == 1):
                if (type_edg[col] == 'h'):
                    denom=denom+' +E_'+line_label_h(col)
                else:
                    denom=denom+' +E_'+line_label_p(col)
            if (val_test == -1):
                if (type_edg[col] == 'h'):
                    denom=denom+'-E_'+line_label_h(col)
                else:
                    denom=denom+'-E_'+line_label_p(col)
        denom = denom+')'
    if ('( +' in denom):
        denom = denom.replace('( +','(')
    denom = denom.strip(' ')
    denoms.append(denom)
    phases.append('(-1)^{%i' % n_holes + '+l}')
    #print incidence
    eq_lines=np.array(incidence.transpose())
    #neq_lines=np.asarray(list(i for i in set(map(tuple,eq_lines)))).transpose()
    neq_lines=np.asarray(list(i for i in set(map(tuple,eq_lines))))
    n_sym = len(eq_lines)-len(neq_lines)
    #### CAVEAT !!! Valid only for *MBPT*
    nedges_eq.append(2**n_sym)
    #print "After neqlines"
    #### Loops


## Function generating the feynmanmf instructions
def feynmf_generator(diag,theory,diag_name):
    p_order = diag.number_of_nodes()
    diag_size = 20*p_order

    theories = ["MBPT","BMBPT","SCGF"]
    th_index = theories.index(theory)
    prop_types = ["half_prop","prop_pm","double_arrow"]
    prop = prop_types[th_index]

    file_name = diag_name + ".tex"

    feynmf_file = open(file_name,'w')

    begin_file = "\parbox{%i" %diag_size +"pt}{\\begin{fmffile}{" + diag_name + "}\n\\begin{fmfgraph*}(%i" %diag_size + ",%i)\n" %diag_size
    end_file = "\end{fmfgraph*}\n\end{fmffile}}\n\n"

    if p_order >= 2:
        feynmf_file.write(begin_file)
        feynmf_file.write("\\fmftop{v%i}\\fmfbottom{v0}\n" %(p_order-1))
        feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v0}\n")
        feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" %(p_order-1))

        if p_order > 2:
            feynmf_file.write("\\fmf{phantom}{v0,v1}\n")
            for vertex in range(1,p_order-2):
                feynmf_file.write("\\fmf{phantom}{v%i" %vertex + ",v%i}\n" %(vertex+1))
                feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" %vertex)

            feynmf_file.write("\\fmfv{d.shape=circle,d.filled=full,d.size=3thick}{v%i}\n" %(p_order-2))
            feynmf_file.write("\\fmf{phantom}{v%i,"%(p_order-2) + "v%i}\n" %(p_order-1))
            feynmf_file.write("\\fmffreeze\n")

        oriented_adj_mat = []
        for i in range(0,p_order):
            oriented_adj_mat.append([])
            for j in range(0,p_order):
                oriented_adj_mat[i].append(0)

        for line in nx.generate_edgelist(diag,data=False):
            i = int(line[0])
            j = int(line[2])
            oriented_adj_mat[i][j] += 1

        for i in range(0,p_order):
            for j in range(0,p_order):
                if (abs(i-j) == 1) and (oriented_adj_mat[i][j] != 0): ## Vertex consecutifs
                    if oriented_adj_mat[i][j] == 1:
                        if oriented_adj_mat[j][i] !=1:
                            feynmf_file.write("\\fmf{" + prop + "}{v%i," %j + "v%i}\n" %i)
                        else:
                            feynmf_file.write("\\fmf{" + prop + ",right=0.5}{v%i," %j + "v%i}\n" %i)
                    else:
                        feynmf_file.write("\\fmf{" + prop + ",right=0.5}{v%i," %j + "v%i}\n" %i)
                        feynmf_file.write("\\fmf{" + prop + ",left=0.5}{v%i," %j + "v%i}\n" %i)

                        if oriented_adj_mat[i][j] == 3:
                            feynmf_file.write("\\fmf{" + prop + "}{v%i," %j + "v%i}\n" %i)
                        elif oriented_adj_mat[i][j] == 4:
                            feynmf_file.write("\\fmf{" + prop + ",right=0.75}{v%i," %j + "v%i}\n" %i)
                            feynmf_file.write("\\fmf{" + prop + ",left=0.75}{v%i," %j + "v%i}\n" %i)

                elif (i != j) and (oriented_adj_mat[i][j] != 0): ## Vertex non consecutifs non diagonaux
                    if (oriented_adj_mat[i][j] == 1) and (oriented_adj_mat[j][i] == 2):
                        feynmf_file.write("\\fmf{" + prop + ",right=0.75}{v%i," %j + "v%i}\n" %i)
                    else:
                        feynmf_file.write("\\fmf{" + prop + ",right=0.6}{v%i," %j + "v%i}\n" %i)
                        if oriented_adj_mat[i][j] != 1:
                            feynmf_file.write("\\fmf{" + prop + ",left=0.6}{v%i," %j + "v%i}\n" %i)
                            if oriented_adj_mat[i][j] != 2:
                                feynmf_file.write("\\fmf{" + prop + ",right=0.75}{v%i," %j + "v%i}\n" %i)
        feynmf_file.write(end_file)
    else:
        print "Perturbative order too small"

## Writing a dot file for each graph
msg = 'Generate diagrams feymanmf instructions ?'
pdraw = raw_input("%s (y/N) " % msg).lower() == 'y'
if (pdraw):
    shutil.copy('feynmp.mp', directory + '/feynmp.mp')
    shutil.copy('feynmp.sty', directory + '/feynmp.sty')
    for i in range(0,numdiag):
        diag_name = 'diag_%i' %i
        feynmf_generator(G[i],"MBPT",diag_name)
        shutil.move(diag_name +'.tex', directory + "/Diagrams/" + diag_name + '.tex')


msg = 'Include diagrams in tex ?'
pdiag = raw_input("%s (y/N) " % msg).lower() == 'y'
### Latexisation
header = "\documentclass[10pt,a4paper]{article}\n \usepackage[utf8]{inputenc}\n\usepackage{braket}\n\usepackage{graphicx}\n"
header = header + "\usepackage[english]{babel}\n\usepackage{amsmath}\n\usepackage{amsfonts}\n\usepackage{amssymb}\n"
if pdiag:
    header = header + "\usepackage[force]{feynmp-auto}\n"
land = False
if (norder > 3):
    msg = 'Expressions may be long rotate pdf ?'
    land = raw_input("%s (y/N) " % msg).lower() == 'y'
if (land):
    header = header + "\usepackage[landscape]{geometry}\n"

header = header + "\\title{Diagrams and algebraic expressions at order %i" % norder +" in MBPT}\n"
latex_file = open(directory + '/result.tex','w')
latex_file.write(header)
begdoc ="\\begin{document}\n"
enddoc ="\\end{document}"
begeq = "\\begin{equation}\n"
endeq = "\\end{equation}\n"
latex_file.write(begdoc)
latex_file.write("\maketitle\n")
latex_file.write("\\graphicspath{{Diagrams/}}")
if (not pdiag or not pdraw):
    for i_diag in range(0,numdiag):
        diag_exp = "\dfrac{1}{%i}" % nedges_eq[i_diag]+phases[i_diag]+"\sum{\dfrac{"+mat_els[i_diag]+"}{"+denoms[i_diag]+"}}\n"
        latex_file.write(begeq)
        latex_file.write(diag_exp)
        latex_file.write(endeq)
    latex_file.write(enddoc)
else:
    for i_diag in range(0,numdiag):
        diag_exp = "\dfrac{1}{%i}" % nedges_eq[i_diag]+phases[i_diag]+"\sum{\dfrac{"+mat_els[i_diag]+"}{"+denoms[i_diag]+"}}\n"
        latex_file.write(begeq)
        latex_file.write(diag_exp)
        latex_file.write(endeq)
        latex_file.write('\n\\begin{center}\n')
        diag_file = open(directory+"/Diagrams/diag_%i.tex" %i_diag)
        latex_file.write(diag_file.read())
        latex_file.write('\\end{center}\n')
    latex_file.write(enddoc)
latex_file.close()

msg = 'Compile pdf ?'
pdfcompile = raw_input("%s (y/N) " % msg).lower() == 'y'
if (pdfcompile):
    os.chdir(directory)
    os.system("pdflatex result.tex")
    if pdiag:
        os.system("pdflatex result.tex")
        for i_diag in range(0,numdiag):
            os.unlink("diag_%i.1" %i_diag)
            os.unlink("diag_%i.mp" %i_diag)
            os.unlink("diag_%i.log" %i_diag)
    print "Result saved in "+directory +'/result.pdf'
