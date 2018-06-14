Generate diagrams with ADG
==========================

Run ADG
--------

To run the program and generate BMBPT diagrams at order 4 for example, use

.. code:: bash

  adg -o 4 -t BMBPT -d -c

where the ```-o``` flag is for the order, ```-t``` for the type of theory,
```-d``` indicates you want the diagrams to be drawn and ```-c``` that you want
ADG to compile the LaTeX output.

You can alternatively run the program in interactive mode by typing

.. code:: bash

  adg -i

Finally, to obtain more information on all the available flags, use

.. code:: bash

  adg -h


CLI options
-----------

Generic options:
*****************

-o, --order         order of the diagrams [1-9]
-t, --theory        theory of interest: MBPT or BMBPT
-i, --interactive   execute ADG in interactive mode

BMBPT options:
**************

-can, --canonical      consider only canonical diagrams
-3N, --with_three_body  use two and three-body forces for BMBPT diagrams
-dt, --draw_tsds        draw Time-Structure Diagrams

MBPT option:
************

-cd, --cd_output  produce output for C. Drischler's framework

Run management options:
***********************

-d, --draw_diags  draw the diagrams using FeynMF
-c, --compile     compile the LaTeX output file with PDFLaTeX


Output files
------------

The output of the program is stored in a folder named after the theory, and a
subfolder named after the order, i.e. for our previous example, results would
be stored under ``/BMBPT/Order-4``. In the case of BMBPT, suffixes are added
for diagrams computed using three-body forces.

The main output file of the program, called ``result.tex``, is a LaTeX file
containing the expressions of the diagrams along other basic infos on their
structure, and, if flag ``-d`` has been used, drawing instructions. The file
is automatically compiled and produces a PDF file ``result.pdf`` when using the
``-c`` file.

A list of the adjacency matrices associated with the diagrams is printed
separately in the ``adj_matrices.list`` file to allow for an easy use with
another many-body diagrams code.

In the case of a MBPT calculations, it is possible to produce output
specifically tailored for C. Drischler's automated calculations framework by
using the ``-cd`` flag. The associated output files use ``CD_`` as a prefix.
