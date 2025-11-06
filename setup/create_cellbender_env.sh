#!/bin/bash

conda create -n cellbender-test python=3.9
source activate cellbender-test
conda install -c anaconda pytables
pip install torch
pip install -e /mnt/z/Caelen/CellBender