#!/bin/bash -f

echo 12/11/2019 was done an upgrade to Python3.7, Qt5 and Matplolib 3. Despite a detailed testing there can be still remaining bugs. Let me know if you find any.
echo Tomas Odstrcil
echo 
echo In case of any problem, please contact: tomas.odstrcil@gmail.com or git@ipp.mpg.de for local issues

source /etc/profile.d/modules.sh

module purge 
module load mkl intel ffmpeg anaconda/3/2019.03

export LD_LIBRARY_PATH=${MKL_HOME}/lib/intel64_lin:/afs/ipp-garching.mpg.de/home/t/todstrci/SuiteSparse/lib/
export C_INCLUDE_PATH=${LD_LIBRARY_PATH}:/afs/ipp-garching.mpg.de/home/t/todstrci/SuiteSparse/lib/

#load mencoder, used only in nogui.py

export PATH=${PATH}:/afs/@cell/common/soft/visualization/mencoder/svn-2012-11-15/amd64_sles11/bin

python /afs/ipp/home/g/git/python/specviewer/pyspecview.py  $argv
