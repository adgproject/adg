### Goal ###
The aim of this project is to provide a tool generating and calculating diagrams
for a given formalism. This formalism reduces to a set of rules which constrains
the generated adjacency matrices' structure.

### Status ###
As for now, the code is capable of handling two differents formalisms, i.e. MBPT
and BMBPT. On the long run, its extension to SCGF might be considered.
    - For MBPT, the code generates all canonical (i.e. HF) diagrams at any given
    order.
    - For BMBPT, the code generates all diagrams for norm or operator kernels at
    any order. Its is capable of handling kernekl operator Feynman expressions
    up to a sign and a symmetry factor, as well as the Goldstone expressions for
    some specific diagrams, up to the same restrictions. Extension of the
    expression generation ot norm diagrams will be considered in a near future.

### Dependencies ###
In order to run the code, you will need a Python install > 2.7.1
  - Python libraries:
  	* networkx (for graph treatment)
	* numpy

### Use ###
(To be modified for faster runs)
To run the program and generate the diagrams
	python adg.py

(Last README update: April 28th 2017)
