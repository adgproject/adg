# ADG - Automatic Diagram Generator

## Description
ADG is a tool generating diagrams and producing their expressions for given
many-body formalisms. Diagrammatic rules from the formalism are combined with
graph theory objects to produce diagrams and expressions in a fast, simple and
error-safe way.

## Status
As for now, the code is capable of handling two differents formalisms, i.e.
Many-Body Perturbation Theory (MBPT) and Bogoliubov Many-Body Perturbation
Theory (BMBPT).
  - For MBPT, the code generates all canonical (i.e. HF) diagrams at any given
    order along with their expression and additional information
    (conjugate diagram, excitation level...).
  - For BMBPT, the code generates all diagrams for norm or operator kernels,
    and the time-dependent and time-integrated expressions for all operator
    diagrams.

## Future developments
Extensions under discussions are norm diagrams for BMBPT as well as diagrams
and expressions generation for Gorkov Self-Consistent Green's Functions (GSCGF).

## Install
To install ADG starting, download the source files and run
```
pip install <project_folder>
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

## Citing
If you use ADG in your research work, we kindly ask you to cite the following
paper: []

## License
ADG is licensed under
