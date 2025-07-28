 
 

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fusion/projects/codes/pyspecview/FFTW3/lib/
##export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fusion/projects/codes/pyspecview/FFTW3/include/
#export C_INCLUDE_PATH=$C_INCLUDE_PATH:/fusion/projects/codes/pyspecview/FFTW3/lib/
#export C_INCLUDE_PATH=$C_INCLUDE_PATH:/fusion/projects/codes/pyspecview/FFTW3/include/
#export LIBRARY_PATH=$LIBRARY_PATH:/fusion/projects/codes/pyspecview/FFTW3/lib/
#export LIBRARY_PATH=$LIBRARY_PATH:/fusion/projects/codes/pyspecview/FFTW3/include/

 
#export PYTHON=/fusion/usc/opt/python/2.7.11/bin/python2.7 
#export PYTHONPATH="${PYTHONPATH}:/fusion/usc/lib"
#export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/fusion/usc/opt/python/2.7.11/lib"
export PATH="${PATH}:/fusion/projects/codes/pyspecview/"

module purge  
#module load mdsplus 
#module load fftw
module load omfit/unstable


echo ' If you find any bugs, please sent email to odstrcilt@fusion.gat.com'

python  /fusion/projects/codes/pyspecview/pyspecview.py  $@
