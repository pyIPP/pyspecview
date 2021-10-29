#!/bin/bash -f
 
source /etc/profile.d/modules.sh

module purge 
module load anaconda/3 aug_sfutils/0.4.8

export LD_LIBRARY_PATH=/afs/ipp/home/t/todstrci/SuiteSparse/lib

export PYTOMO=/afs/ipp/home/g/git/python/tomo_dev
export PYSPECVIEW=/afs/ipp/home/g/git/python/specview_dev

$PYSPECVIEW/specview.py $@
