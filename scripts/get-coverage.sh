#!/bin/bash

micromamba activate offqcdata-test

python get-low-coverage-torsions.py > logs/get-low-coverage-torsions.log 2>&1
python get-coupled-torsions.py > logs/get-coupled-torsions.log 2>&1
