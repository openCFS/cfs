#!/bin/sh

SOURCES_DIR=$1
DEST_DIR=$2
STYLE=tango
OPTIONS="-f latex -O linenos=True,linenostep=5,style=$STYLE"
PYG=pygmentize

mkdir -p $DEST_DIR

$PYG -f latex -O full,style=$STYLE $SOURCES_DIR/cfs_interpolation.xml > $DEST_DIR/dummy.tex
# Grep the definitions for the fancy verbatim environments from dummy.tex (for older pygmentize)
fgrep newcommand $DEST_DIR/dummy.tex > $DEST_DIR/pygmentize.tex
# Extract the definitions (for newer pygmentize)
cat $DEST_DIR/dummy.tex | sed -e '/./{H;$!d;}' -e 'x;/def/!d;' >> $DEST_DIR/pygmentize.tex
#echo $SOURCES_DIR
#echo $DEST_DIR

function Pygmentize {
  INPUT=$1
  OUTPUT=$2
 
  $PYG $OPTIONS $SOURCES_DIR/$INPUT > $DEST_DIR/dummy.tex
  cat $DEST_DIR/dummy.tex | sed 's/commandchars/frame=lines,framesep=2mm,fontsize=\\\small,commandchars/' > $DEST_DIR/$OUTPUT
}

Pygmentize cfs_interpolation.xml cfs_interpolation.tex
Pygmentize gengrids_spherical_shells.xml gengrids_spherical_shells.tex
Pygmentize gengrids_rectlin.xml gengrids_rectlin.tex
Pygmentize gengrids_refelems.xml gengrids_refelems.tex
Pygmentize gengrids_refelems.xml gengrids_refelems.tex
Pygmentize gengrids_fields.xml gengrids_fields.tex
Pygmentize check_cons_interp.m check_cons_interp.tex
Pygmentize sumup.m sumup.tex