#!/bin/bash

SOURCES_DIR=$1
DEST_DIR="$2/pygmentized"
ATTACH_DIR="$2/attachments"
STYLE=tango
OPTIONS="-f latex -O linenos=True,linenostep=5,style=$STYLE"
PYG=pygmentize

mkdir -p $DEST_DIR

$PYG -f latex -O full,style=$STYLE $SOURCES_DIR/hdf5_file_format.xml > $DEST_DIR/dummy.tex
# Grep the definitions for the fancy verbatim environments from dummy.tex (for older pygmentize)
fgrep newcommand $DEST_DIR/dummy.tex > $DEST_DIR/pygmentize.tex
# Extract the definitions (for newer pygmentize)
cat $DEST_DIR/dummy.tex | sed -e '/./{H;$!d;}' -e 'x;/def/!d;' >> $DEST_DIR/pygmentize.tex
#echo $SOURCES_DIR
#echo $DEST_DIR

function Pygmentize {
  INPUT=$1
  OUTPUT=$2
 
  if [ ! -d $DEST_DIR ]; then
      mkdir -p $DEST_DIR
  fi
  if [ ! -d $ATTACH_DIR ]; then
      mkdir -p $ATTACH_DIR
  fi

  $PYG $OPTIONS $SOURCES_DIR/$INPUT > $DEST_DIR/dummy.tex
  cat $DEST_DIR/dummy.tex | sed 's/commandchars/frame=lines,framesep=2mm,fontsize=\\\small,commandchars/' > $DEST_DIR/$OUTPUT
  
  ATTACHFILE=$(echo $OUTPUT | sed 's/\.tex/\.txt/')
  cp $SOURCES_DIR/$INPUT "$ATTACH_DIR/$ATTACHFILE"
}

Pygmentize hdf5_file_format.xml hdf5_file_format.tex
Pygmentize crosspoint_test.c gmsh_cp_test.tex
Pygmentize photo.c gmsh_mesh_from_pic.tex
Pygmentize gid_batch.tcl gid_batch.tex
Pygmentize gid_pw_server.sh gid_pw_server.tex
Pygmentize pv_analytical_sol.py pv_analytical_sol.tex
Pygmentize pv_scipy_hankel1.py pv_scipy_hankel1.tex
Pygmentize pv_matlab_besselh.m pv_matlab_besselh.tex
Pygmentize icem_to_ansys.in.txt icem_to_ansys.in.tex
Pygmentize Pfeifenfussmodell_modificated_swept.in.txt Pfeifenfussmodell_modificated_swept.in.tex
Pygmentize mkmesh.txt mkmesh.tex
Pygmentize mkmesh.mlib.txt mkmesh.mlib.tex
Pygmentize ansys_surf_vol_mesh.txt ansys_surf_vol_mesh.tex
Pygmentize cplreader_meshgen.xml cplreader_meshgen.tex
Pygmentize qs_obtain_sources.sh qs_obtain_sources.tex
Pygmentize git_checkout_cfs.sh git_checkout_cfs.tex
Pygmentize qs_configure_build.sh qs_configure_build.tex
Pygmentize pv_tria_local_global_mapping.py pv_tria_local_global_mapping.tex
Pygmentize pv_xi_eta.py pv_xi_eta.tex
