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
import itertools
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


print adjcoeff,type(adjcoeff[1])

for item in itertools.combinations_with_replacement(adjcoeff, r=norder):
    flattupl = sum(item,())
    adjmtrx.append(flattupl+tuple([0]*(norder**2-len(flattupl))))

print adjmtrx

#for partition in adjcoeff:
#    adjmtrx.append((partition*norder+tuple([0]*(norder**2-len(partition*norder)))))


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
numdiag = len(diagrams)
print "Number of diagrams, ",numdiag

### Graph part (computing, writing, drawing)        
G=[]
sizegraph = norder*100 
for diagram in diagrams:
    G.append(nx.from_numpy_matrix(diagram,create_using=nx.MultiDiGraph()))


for i in range(0,numdiag):
    pos = []
    for vertex in range(0,norder):
        if (vertex == 0):
           pos.append("%i" % vertex +' [pos = "0,0",shape=circle,fixedsize=true,width=1.2] \n ')
        else:
           position = sizegraph*(vertex)/(norder-1)
           pos.append( "%i" % vertex +' [pos = "0,%i"' % position + ',shape=circle,fixedsize=true,width=1]\n ')
for i in range(0,numdiag):
    write_dot(G[i],'diag_%i.dot' % i)

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

for i in range(0,numdiag):
    
    line = "digraph  {\n"+'splines=true \n sep = 1\n overlap = false;ratio="fill";margin=0;\n'+' node[label =""]\n'.join(pos)
    replace_line('diag_%i.dot' % i,0,line)
#    for line in fileinput.FileInput('diag_%i.dot' % i,inplace=1):
#        p += 1
#        if "digraph" in line:
#            line=line.replace(line,line+'splines=true \n sep = 1\n overlap = false;ratio="fill";margin=0;\n' + 
#            ' node[label =""]\n'.join(pos)) 
#            print line


for i in range(0,numdiag):
    #os.system("neato diag_%i.dot" % i + " -n -Tpng -o diag_%i.png" %i)
    os.system("neato diag_%i.dot" % i + " -n -Tpng -o diag_%i.png" %i)


