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


print "#####################"
print "# Automatic Diagram #"
print "#     Generator     #"
print "#    RDL,JR,PA,MD   #"
print "#####################"

print "Parallel Mode"
num_cores = multiprocessing.cpu_count()
print "There is %i" % num_cores + " core(s) available"
norder = int(raw_input('Order of the diagrams ?\n'))

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
with open("Diagrams.list", "w") as f:
    for diagram in diagrams:
        f.write("Diagram n: %i" % i)
        #np.savetxt(f,diagram)
        diagram.tofile(f,"")
        f.write("\n")
        i += 1

### Graph part (computing, writing, drawing)        
G=[]
sizegraph = norder*100
for diagram in diagrams:
    G.append(nx.from_numpy_matrix(diagram,create_using=nx.MultiDiGraph()))
print diagrams[1]
G1=[]
for diag in G:
    #if (not nx.is_strongly_connected(diag)):
    if((nx.number_strongly_connected_components(diag)) == 1):
        G1.append(diag)
G=G1
numdiag = len(G)
print "Time ellapsed: ",datetime.now() - start_time
print "Number of connected diagrams, ",numdiag

## Optimizing the position of each vertex for the set of diagrams
for i in range(0,numdiag):
    pos = []
    for vertex in range(0,norder):
        if (vertex == 0):
           pos.append("%i" % vertex +' [pos = "0,0",shape=circle,fixedsize=true,width=1.2] \n ')
        else:
           position = sizegraph*(vertex)/(norder-1)
           pos.append( "%i" % vertex +' [pos = "0,%i"' % position + ',shape=circle,fixedsize=true,width=1]\n ')

## Writing a dot file for each graph
for i in range(0,numdiag):
    write_dot(G[i],'diag_%i.dot' % i)

## Function to replace a specific line of a textfile
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

## Plotting features
for i in range(0,numdiag):
    line = "digraph  {\n"+'splines=true \n sep = 1\n overlap = false;ratio="fill";margin=0;\n'+' node[label =""]\n'.join(pos)
    replace_line('diag_%i.dot' % i,0,line)

## Printing
for i in range(0,numdiag):
    os.system("neato diag_%i.dot" % i + " -n -Tpng -o diag_%i.png" %i)

