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
from networkx.drawing.nx_agraph import write_dot
import fileinput
import itertools
#from joblib import Parallel, delayed
import multiprocessing
from datetime import datetime
import string


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
sizegraph = norder*100
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

def line_label(n):
    labels=list(string.ascii_lowercase)
    return labels[n]

def mat_elements(irow):
    return 

    
###

mat_els = []
denoms = []
phases = []
nedges_eq = []
for diag in G:
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
        diffcols.add(repr(flat))

    for row in range(nrow):
        ket = ''
        bra = ''
        for col in range(ncol):
        ######### Mtrx Elements ###########
            if (incidence[row,col] == 1):
                bra = bra + line_label(col) 
            if (incidence[row,col] == -1):
                ket = ket + line_label(col)
        ###################################
        braket = braket + '\\braket{'+bra+'|H|'+ket+'}'
    mat_els.append(braket)
    denom = ''
    for row in range(1,nrow):
        denom = denom + '('
        for col in range(ncol):
            val_test = incidence[0:row,col].sum()
            if (val_test == 1):
                denom=denom+' +E_'+line_label(col)
            if (val_test == -1):
                denom=denom+'-E_'+line_label(col)
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


## Optimizing the position of each vertex for the set of diagrams
for i in range(0,numdiag):
    pos = []
    for vertex in range(0,norder):
        if (vertex == 0):
           pos.append("%i" % vertex +' [pos = "0,0",shape=circle,fixedsize=true,width=1] \n ')
        else:
           position = sizegraph*(vertex)/(norder-1)
           pos.append( "%i" % vertex +' [pos = "0,%i"' % position + ',shape=circle,fixedsize=true,width=1]\n ')

## Writing a dot file for each graph
msg = 'Generate and plot diagrams ?'
pdraw = raw_input("%s (y/N) " % msg).lower() == 'y'
if (pdraw):
    for i in range(0,numdiag):
        write_dot(G[i],directory+'/Diagrams/diag_%i.dot' % i)

## Function to replace a specific line of a textfile
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

## Plotting features
if (pdraw):
    for i in range(0,numdiag):
        line = "digraph  {\n"+'splines=true \n sep = 1\n overlap = false;ratio="fill";margin=0;\n'+' node[label =""]\n'.join(pos)
        replace_line(directory+'/Diagrams/diag_%i.dot' % i,0,line)

## Printing
    for i in range(0,numdiag):
        os.system("dot " +directory+"/Diagrams/diag_%i.dot" % i + " -Tpng -o"+directory+"/Diagrams/diag_%i.png" %i)
        ## "Pretty but misleading"
        #os.system("neato " +directory+"/Diagrams/diag_%i.dot" % i + " -n -Tpng -o"+directory+"/Diagrams/diag_%i.png" %i)



msg = 'Include diagrams in tex ?'
pdiag = raw_input("%s (y/N) " % msg).lower() == 'y'
### Latexisation
header = "\documentclass[10pt,a4paper]{article}\n \usepackage[utf8]{inputenc}\n\usepackage{braket}\n\usepackage{graphicx}\n"
header = header + "\usepackage[english]{babel}\n\usepackage{amsmath}\n\usepackage{amsfonts}\n\usepackage{amssymb}\n"
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
        latex_file.write('\\begin{center}\n')
        latex_file.write("\\includegraphics[scale=0.15]{diag_%i" %i_diag +".png}\n")
        latex_file.write('\\end{center}\n')
    latex_file.write(enddoc)
latex_file.close()

msg = 'Compile pdf ?'
pdfcompile = raw_input("%s (y/N) " % msg).lower() == 'y'
if (pdfcompile):
    os.chdir(directory) 
    os.system("pdflatex result.tex")
    print "Result saved in "+directory +'/result.pdf'
    


