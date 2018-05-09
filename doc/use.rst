Generate diagrams with ADG
==========================

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
