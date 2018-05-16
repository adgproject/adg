The ADG Project
===============

Description
-----------
ADG is a tool generating diagrams and producing their expressions for given
many-body formalisms. Diagrammatic rules from the formalism are combined with
graph theory objects to produce diagrams and expressions in a fast, simple and
error-safe way.

Status
------
As for now, the code is capable of handling two differents formalisms, i.e.
Many-Body Perturbation Theory (MBPT) and Bogoliubov Many-Body Perturbation
Theory (BMBPT).
  - For MBPT, the code generates all canonical (i.e. HF) diagrams at any given
    order along with their expression and additional information
    (conjugate diagram, excitation level...).
  - For BMBPT, the code generates all diagrams for norm or operator kernels,
    and the time-dependent and time-integrated expressions for all operator
    diagrams.

Future developments
-------------------
Extensions under discussions are norm diagrams for BMBPT as well as diagrams
and expressions generation for Gorkov Self-Consistent Green's Functions (GSCGF).
