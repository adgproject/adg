Install ADG on your computer
============================

Install
--------
To install ADG starting, download the source files and run

.. code:: bash

  pip install <project_folder>


If you want to install ADG in ``develop`` mode, then run

.. code:: bash

  pip install -e <project_folder>


Dependencies
------------
In order to run the code, you will need a Python install >= 2.7.1

  - Python libraries:

  	* networkx >= 2.0
    
    * numpy

If you want ADG to compile the LaTeX output file, you will need a Latex install
with the PDFLaTeX compiler and the feynmp and feynmp-auto packages installed,
which are standard packages in most recent distributions.
