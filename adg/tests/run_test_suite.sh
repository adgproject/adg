#!/bin/bash

python -m pytest --doctest-modules -v ../*.py ./*.py --cov-report=html --cov=adg -vv
