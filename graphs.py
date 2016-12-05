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
numpart = len(adjcoeff)

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
print 'numpart'
print numpart


for partition in adjcoeff:
    adjmtrx.append((partition*norder+tuple([0]*(norder**2-len(partition*norder)))))


icount = 0
#for element in adjcoeff:
#    #element = list(element)
#    element = deque(element)
#    element.extendleft([0] *(norder - len(element)))
#    print element
#


class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

def perm_unique(elements):
    eset=set(elements)
    listunique = [unique_element(i,elements.count(i)) for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

def perm_unique_helper(listunique,result_list,d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  perm_unique_helper(listunique,result_list,d-1):
                    yield g
                i.occurrences+=1


#poss_lines = [[] for _ in range(norder)]
#poss_lines.remove(0)





    
#    for line in range(0,norder):
#            adjmtrx[icount][line] = element 
#    icount += 1

#print adjmtrx
print "Size,",len(adjmtrx[0])
print adjmtrx[0]
diagrams = []
print 'Tri en cours'
for diagram in adjmtrx:
    permuts = perm_unique(diagram)
    for permutation in permuts:
        matpermut = np.array(permutation).reshape(norder,norder)
        sumcol = []
        sumrow = []
        for col in range(0,norder):
            sumcol.append(matpermut[:,col].sum())
            sumrow.append(matpermut[col,:].sum())
        t1 = all(va == nbody for va in sumcol)
        t2 = (matpermut.trace() == 0)
        t3 = all(va == nbody for va in sumrow)
        if (t1 and t2 and t3):
            print "Done !"
            print matpermut
            diagrams.append(matpermut)
print 'original'
print adjmtrx[0]
print "Number of diagrams, ",len(diagrams)
diaginit = copy.deepcopy(adjmtrx[0])
###Column permutations
#while True:
#    #diagpermut = copy.deepcopy(adjmtrx[1])
#    diagpermut = copy.deepcopy(adjmtrx[0])
#    print diagpermut
#    sumcol2 =[]
#    t1 = all(va == nbody for va in sumcol2)
#    t2 = (diagpermut.trace() == 0)
#    t3 = (diagpermut[line].sum() == nbody)
#    eq1 = np.array_equal(diagpermut,diaginit)
#    print eq1
#    #if (not t1 or not t2 or not t3 or eq1):
#    if (not eq1):
#        print 'plop'
#        print diagpermut
#        for col in range(0,norder):
#            random.shuffle(diagpermut[:,col])
#        print 'apres'
#        print diagpermut
#        continue
#    else:
#        break
#

### Graph part (computing, writing, drawing)        
G=[]
for diagram in diagrams:
    G.append(nx.from_numpy_matrix(diagram,create_using=nx.MultiDiGraph()))

for i in range(0,numpart):
    write_dot(G[i],'diag_%i.dot' % i)

for i in range(0,numpart):
    for line in fileinput.FileInput('diag_%i.dot' % i,inplace=1):
        if "digraph" in line:
            line=line.replace(line,line+'splines=true \n sep = 1\n overlap = false;ratio="fill";margin=0;\n' + 
            ' node[label =""]\n 0 [pos = "0,0",shape=circle,fixedsize=true,width=1.2] \n 1 [pos = "0,150",shape=circle,fixedsize=true,width=1]\n 2 [pos = "0,300",shape=circle,fixedsize=true,width=1.2] \n')
        print line,

for i in range(0,numpart):
    os.system("neato diag_%i.dot" % i + " -n -Tpng -o diag_%i.png" %i)


