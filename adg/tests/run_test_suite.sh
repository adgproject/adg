#!/bin/bash

py.test --doctest-modules -v ../*.py ./*.py --cov-report=html --cov=adg -vv
