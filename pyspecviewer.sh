#!/bin/bash
module purge  


module load omfit 
  
echo ' If you find any bugs, please sent email to odstrcilt@fusion.gat.com'

export PYTHONPATH=/fusion/projects/codes/pyspecview/pyfftw:$PYTHONPATH

export PATH=/fusion/projects/codes/pyspecview/:${PATH}

python  /fusion/projects/codes/pyspecview/pyspecview.py  $@


 
