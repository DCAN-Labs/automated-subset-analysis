#! /usr/bin/env bash

/home/exacloud/lustre1/fnl_lab/code/external/GUIs/MATLAB/R2018a/bin/mcc -v -m -R -singleCompThread -R -nodisplay -o MultiShadedBars MultiShadedBars.m -a  /home/users/conan/asa-sbatch-2020-01-02/src/shadedbars.m -a /home/users/conan/asa-sbatch-2020-01-02/src/loadcorr.m

#add MCR_CACHE_ROOT to all run scripts for Exahead processing
sed -i '/exe_dir=`dirname "$0"`/a if [ ! -d $TMPDIR/$USER ]; then\n    mkdir $TMPDIR/$USER\nfi\nexport MCR_CACHE_ROOT=$TMPDIR/$USER' run_MultiShadedBars.sh
