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


def partition(n):
    """ A function which computes the numbers of partitions and return them into an list """
    part_of = []
    part_of.append([()])
    part_of.append([(1,)])
    for num in range(2, n+1):
        ptitions = set()
        for i in range(num):
            for partition in part_of[i]:
                ptitions.add(tuple(sorted((num - i, ) + partition)))
        part_of.append(list(ptitions))
    return part_of[n]


nbody=int(raw_input('N-body interaction where n is ?\n'))
adjcoeff = (partition(nbody))
## Pretty printing but useless af
for k in range(0,len(adjcoeff)):
    if (len(adjcoeff[k]) == 1):
        adjcoeff[k] += (0,)
print ('There is %i' % len(adjcoeff) + ' partitions')
print adjcoeff
numdiag = len(adjcoeff)
print numdiag

def swap(xs, a, b):
    xs[a], xs[b] = xs[b], xs[a]


def derange(xs):
    for a in xrange(1, len(xs)):
        b = random.choice(xrange(0, a))
        swap(xs, a, b)
    return xs

norder = int(raw_input('Order of the diagrams ?\n'))
adjmtrx = []
## Initialisation
print 'numdiag'
print numdiag
for i in range (0,numdiag):
    adjmtrx.append(np.zeros((norder, norder)))


icount = 0
for element in adjcoeff:
    #element = list(element)
    element = deque(element)
    element.extendleft([0] *(norder - len(element)))
    for line in range(0,norder):
            adjmtrx[icount][line] = element 
    icount += 1

print adjmtrx
print 'Tri en cours'
for diagram in adjmtrx:
    while True:
        sumcol = []
        for col in range(0,norder):
            sumcol.append(diagram[:,col].sum())
        t1 = all(va == nbody for va in sumcol)
        t2 = (diagram.trace() == 0)
        t3 = (diagram[line].sum() == nbody)
        #print "Somme col",t1,"Trace nulle ?",t2,"Somme sur ligne ?",t3
        if (not t1 or not t2 or not t3):
            sumcol = []
            for line in range(0,norder):
                random.shuffle(diagram[line])
            continue
        else:
            print "Done !"
            break
print 'original'
print adjmtrx[0]
diaginit = copy.deepcopy(adjmtrx[0])
###Column permutations
while True:
    #diagpermut = copy.deepcopy(adjmtrx[1])
    diagpermut = copy.deepcopy(adjmtrx[0])
    print diagpermut
    sumcol2 =[]
    t1 = all(va == nbody for va in sumcol2)
    t2 = (diagpermut.trace() == 0)
    t3 = (diagpermut[line].sum() == nbody)
    eq1 = np.array_equal(diagpermut,diaginit)
    print eq1
    #if (not t1 or not t2 or not t3 or eq1):
    if (not eq1):
        print 'plop'
        print diagpermut
        for col in range(0,norder):
            random.shuffle(diagpermut[:,col])
        print 'apres'
        print diagpermut
        continue
    else:
        break


### Graph part (computing, writing, drawing)        
G=[]
for diagram in adjmtrx:
    G.append(nx.from_numpy_matrix(diagram,create_using=nx.MultiDiGraph()))

for i in range(0,numdiag):
    write_dot(G[i],'diag_%i.dot' % i)

for i in range(0,numdiag):
    for line in fileinput.FileInput('diag_%i.dot' % i,inplace=1):
        if "digraph" in line:
            line=line.replace(line,line+'splines=true \n sep = 1\n overlap = false;ratio="fill";margin=0;\n' + 
            ' node[label =""]\n 0 [pos = "0,0",shape=circle,fixedsize=true,width=1.2] \n 1 [pos = "0,150",shape=circle,fixedsize=true,width=1]\n 2 [pos = "0,300",shape=circle,fixedsize=true,width=1.2] \n')
        print line,

for i in range(0,numdiag):
    os.system("neato diag_%i.dot" % i + " -n -Tpng -o diag_%i.png" %i)


