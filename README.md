# ADG - Automatic Diagram Generator

[![Documentation Status](https://readthedocs.org/projects/adg/badge/?version=master)](https://adg.readthedocs.io/en/master/?badge=master)
[![Build Status](https://travis-ci.com/adgproject/adg.svg?branch=master)](https://travis-ci.com/adgproject/adg)
[![PyPI version](https://img.shields.io/pypi/v/adg.svg)](https://pypi.org/project/adg/)
[![Python version](https://img.shields.io/pypi/pyversions/adg)](https://pypi.org/project/adg/)

## Description
ADG is a tool generating diagrams and producing their expressions for given
many-body formalisms. Diagrammatic rules from the formalism are combined with
graph theory objects to produce diagrams and expressions in a fast, simple and
error-safe way.

The only input consists in the theory and order of interest, and the N-body
character of the operators of interest. The main output is a LaTeX file
containing the diagrams, their associated expressions and additional
informations that can be compiled by ADG if needed. Other computer-readable
files may be produced as well.

## Status
As for now, the code is capable of handling four different formalisms, i.e.
Many-Body Perturbation Theory (MBPT), Bogoliubov Many-Body Perturbation
Theory (BMBPT), Projected Bogoliubov Many-Body Perturbation Theory (PBMBPT)
and Bogoliubov In-Medium Similarity Renormalization Group (BIMSRG).
  - For MBPT, the code generates all Hartree-Fock energy diagrams at any given
    order along with their expression and additional information
    (conjugate diagram, excitation level...).
  - For (P)BMBPT, the code generates all diagrams for a generic observable
    commuting with the Hamiltonian, along with their time-dependent and
    time-integrated expressions.
  - For BIMSRG, the code generates all diagrams and expressions at any given
    truncation order for the two operators as well as the commutator itself.
    The traditional BIMSRG(n) truncation order corresponds to truncating both
    operators as well as the commutator at the same rank.

## Future developments
ADG is currently being extended to diagrams and expressions generation for
Gorkov Self-Consistent Green's Functions (GSCGF).

## Install
The easiest way to install the latest stable version of ADG is to use
```
pip install adg
```
Updating after the release of a new version can be done via
```
pip install --upgrade adg
```

To install ADG after downloading the source files, run
```
pip install <project_folder>
```
or alternatively
```
python setup.py install
```
If you want to install ADG in ```develop``` mode, then run
```
pip install -e <project_folder>
```

## Dependencies
In order to run the code, you will need a Python install >= 2.7.1
  - Python libraries:
  	* networkx >= 2.0
    * numpy
    * scipy
    * future
    * more-itertools

If you want ADG to compile the LaTeX output file, you will need a Latex install
with the PDFLaTeX compiler and the feynmp and feynmp-auto packages installed,
which are standard packages in most recent distributions.


## Use
To run the program and generate BMBPT diagrams at order 4 for example, use
```
adg -o 4 -t BMBPT -d -c
```
where the ```-o``` flag is for the order, ```-t``` for the type of theory,
```-d``` indicates you want the diagrams to be drawn and ```-c``` that you want
ADG to compile the LaTeX output.

You can alternatively run the program in interactive mode by typing
```
adg -i
```

Finally, to obtain more information on all the available flags, use
```
adg -h
```

## Documentation

An extensive on-line documentation is available at https://adg.rtfd.io/.
Alternatively, the documentation can be generated from entering the ```doc```
folder and using
```
make html
```

## Package structure

The main folder contains the README file of the package and the setup.py file
used for installing the package.

### adg folder

This folder contains the various Python files with the functions used by ADG.

### doc folder

This folder contains the documentation of the package in reStructuredText,
as well as the Makefile and the conf.py files used by Sphinx for an automatic
generation of the documentation. Part of the documentation is stored in the adg
subfolder.

The compiled documentation is stored in the build subfolder, in a folder named
after the type of documentation. Especially, the manpages of ADG distributed
with the package are stored in the subfolder /build/man.

### examples folder

This folder contains very simple scripts that can be used to launch test
calculations automatically. The precompiled output corresponding to the tests is
available in the sample_output subfolder, then organized in the same way an
actual calculation output would be stored.

## Citing
If you use ADG in your research work, we kindly ask you to cite the following
papers:
  - P. Arthuis, T. Duguet, A. Tichai, R.-D. Lasseri and J.-P. Ebran,
    [Comput. Phys. Commun. **240**, 202-227](https://doi.org/10.1016/j.cpc.2018.11.023) (2019).
  - P. Arthuis, A. Tichai, J. Ripoche and T. Duguet,
    [Comput. Phys. Commun. **261**, 107677](https://doi.org/10.1016/j.cpc.2020.107677) (2021).
  - A. Tichai, P. Arthuis, H. Hergert and T. Duguet,
    [arXiv:2102.10889](https://arxiv.org/abs/2102.10889) (2021).

## License
ADG is licensed under GNU General Public License version 3 (see LICENSE.txt).
```
Copyright (C) 2018-2021 ADG Dev Team
Pierre Arthuis - TU Darmstadt & ExtreMe Matter Institute EMMI, GSI, Darmstadt (previously University of Surrey & Irfu, CEA, Université Paris-Saclay & CEA, DAM, DIF)
Thomas Duguet - Irfu, CEA, UPSaclay & KU Leuven, IKS
Jean-Paul Ebran - CEA, DAM, DIF
Heiko Hergert - NSCL/FRIB Laboratory & Department of Physics and Astronomy, Michigan State University
Raphaël-David Lasseri - ESNT, Irfu, CEA, UPSaclay (previously IPN, CNRS/IN2P3, UPSud, UPSaclay)
Julien Ripoche - CEA, DAM, DIF
Alexander Tichai - MPI fuer Kernphysik, Heidelberg & TU Darmstadt & ExtreMe Matter Institute EMMI, GSI, Darmstadt (previsously ESNT, Irfu, CEA, Université Paris-Saclay)
```
