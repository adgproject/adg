Install ADG on your computer
============================

Install
--------
To install ADG, download the source files and run

.. code:: bash

  pip install <project_folder>

or alternatively

.. code:: bash

  python setup.py install


If you want to install ADG in ``develop`` mode, then run

.. code:: bash

  pip install -e <project_folder>


Dependencies
------------
In order to run the code, you will need a Python2 install >= 2.7.1 and the
following Python libraries:

  - networkx >= 2.0
  - numpy
  - scipy
  - future

If you want ADG to compile the LaTeX output file, you will need a Latex install
with the PDFLaTeX compiler and the feynmp and feynmp-auto packages installed,
which are standard packages in most recent distributions.
