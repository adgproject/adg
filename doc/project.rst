The ADG Project
===============

Description
-----------
ADG is a tool generating diagrams and producing their expressions for given
many-body formalisms. Diagrammatic rules from the formalism are combined with
graph theory objects to produce diagrams and expressions in a fast, simple and
error-safe way.

The only input consists in the theory and order of interest, and the N-body
character of the operators of interest. The main output is a LaTeX file
containing the diagrams, their associated expressions and additional
informations that can be compiled by ADG id needed. Other computer-readable
files may be produced as well.


Status
------
As for now, the code is capable of handling three different formalisms, i.e.
Many-Body Perturbation Theory (MBPT), Bogoliubov Many-Body Perturbation
Theory (BMBPT). and Projected Bogoliubov Many-Body Perturbation Theory (PBMBPT).
  - For MBPT, the code generates all Hartree-Fock energy diagrams at any given
    order along with their expression and additional information
    (conjugate diagram, excitation level...).
  - For (P)BMBPT, the code generates all diagrams for a generic observable
    commuting with the Hamiltonian, along with their time-dependent and
    time-integrated expressions.

Future developments
-------------------
ADG is currently being extended to diagrams and expressions generation for
Gorkov Self-Consistent Green's Functions (GSCGF).
