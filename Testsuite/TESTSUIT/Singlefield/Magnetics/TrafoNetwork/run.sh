#!/bin/sh
#cfs_master=/home/kroppert/Devel/CFS_BIN/hbfemGit_opt/bin/cfs
#hbcfs=/home/kroppert/Devel/CFS_BIN/hbfemGit_opt/bin/cfs
#gmsh=/share/programs/bin/gmsh


# mesh
if [ 1 -eq 0 ]; then
  $gmsh -3 SimpleTrafo.geo
fi

# harmonic run
if [ 1 -eq 0 ]; then
  $cfs_master -p TrafoNetwork-harmonic.xml harmonic
fi

# run the cases with low, medium and high f
NHARMONICS=5
#MAT='analytic-bilinear'
MAT='analytic-s-3'

trans=1 # switch for transient
mh=1 # switch for multiharmonic

# low frequency
ICOIL=3.0
FCOIL=1.0

CASE=$MAT'_f'$FCOIL'-i'$ICOIL

if [ $trans -eq 1 ]; then
  JOB='trans_'$CASE
  sed 's|FCOIL|'$FCOIL'|g' TrafoNetwork-transient.xml | sed 's|ICOIL|'$ICOIL'|g' | sed 's|MAT|'$MAT'|g' > $JOB'.xml'
  $cfs_master $JOB
fi
if [ $mh -eq 1 ]; then
  JOB='mh'$NHARMONICS'_'$CASE
  sed 's|FCOIL|'$FCOIL'|g' TrafoNetwork-multiharmonic.xml | sed 's|ICOIL|'$ICOIL'|g' | sed 's|NHARMONICS|'$NHARMONICS'|g' | sed 's|MAT|'$MAT'|g' > $JOB'.xml'
  $hbcfs $JOB
fi


# medium frequency
ICOIL=3.3
FCOIL=2.5

CASE=$MAT'_f'$FCOIL'-i'$ICOIL

if [ $trans -eq 1 ]; then
  JOB='trans_'$CASE
  sed 's|FCOIL|'$FCOIL'|g' TrafoNetwork-transient.xml | sed 's|ICOIL|'$ICOIL'|g' | sed 's|MAT|'$MAT'|g'  > $JOB'.xml'
  $cfs_master $JOB
fi
if [ $mh -eq 1 ]; then
  JOB='mh'$NHARMONICS'_'$CASE
  sed 's|FCOIL|'$FCOIL'|g' TrafoNetwork-multiharmonic.xml | sed 's|ICOIL|'$ICOIL'|g' | sed 's|NHARMONICS|'$NHARMONICS'|g' | sed 's|MAT|'$MAT'|g' > $JOB'.xml'
  $hbcfs $JOB
fi


# high frequency
ICOIL=7.0
FCOIL=10.0

CASE=$MAT'_f'$FCOIL'-i'$ICOIL

if [ $trans -eq 1 ]; then
  JOB='trans_'$CASE
  sed 's|FCOIL|'$FCOIL'|g' TrafoNetwork-transient.xml | sed 's|ICOIL|'$ICOIL'|g'  | sed 's|MAT|'$MAT'|g' > $JOB'.xml'
  $cfs_master $JOB
fi
if [ $mh -eq 1 ]; then
  JOB='mh'$NHARMONICS'_'$CASE
  sed 's|FCOIL|'$FCOIL'|g' TrafoNetwork-multiharmonic.xml | sed 's|ICOIL|'$ICOIL'|g' | sed 's|NHARMONICS|'$NHARMONICS'|g' | sed 's|MAT|'$MAT'|g' > $JOB'.xml'
  $hbcfs $JOB
fi
