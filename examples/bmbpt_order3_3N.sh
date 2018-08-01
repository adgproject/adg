#! /usr/bin/env bash

# Run the ADG program for BMBPT at 3rd order
# with standard flags using 3-body observable and 3-body forces

adg -t BMBPT -o 3 -3NF -nobs 3 -d -dt -c
