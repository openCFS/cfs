#!/bin/bash -l
#
# Note that the PBS statements are no comments but options for qsub 
# (equivalent to command line options)
# 
# allocate 1 node with 4 cores for given time hours
#PBS -l nodes=1:ppn=4,walltime=06:00:00
#
# send mail for begin (b), end (e), abort (a), e.g. -m abe
#PBS -m e
#PBS -M fabian.wein@fau.de
#
# job name wold be -N my_nme
# stdout and stderr files would be -o job33.out -e job33.err
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME
# if there is no cd ... following generate_qsub_script() from cfs_utils.py will add it
# the last line will get the cfs call by generate_qsub_script
