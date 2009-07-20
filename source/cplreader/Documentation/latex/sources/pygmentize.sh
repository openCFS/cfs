#!/bin/sh

SOURCES_DIR=$1
DEST_DIR=$2
STYLE=borland
OPTIONS="-f latex -O linenos=True,linenostep=5,style=$STYLE"
PYG=pygmentize

mkdir -p $DEST_DIR

$PYG -f latex -O full,style=$STYLE $SOURCES_DIR/cfs_interpolation.xml > $DEST_DIR/dummy.tex
fgrep newcommand $DEST_DIR/dummy.tex > $DEST_DIR/pygmentize.tex
#echo $SOURCES_DIR
#echo $DEST_DIR

$PYG $OPTIONS $SOURCES_DIR/cfs_interpolation.xml > $DEST_DIR/cfs_interpolation.tex
$PYG $OPTIONS $SOURCES_DIR/gengrids_spherical_shells.xml > $DEST_DIR/gengrids_spherical_shells.tex
$PYG $OPTIONS $SOURCES_DIR/gengrids_rectlin.xml > $DEST_DIR/gengrids_rectlin.tex
$PYG $OPTIONS $SOURCES_DIR/gengrids_refelems.xml > $DEST_DIR/gengrids_refelems.tex
$PYG $OPTIONS $SOURCES_DIR/gengrids_refelems.xml > $DEST_DIR/gengrids_refelems.tex
$PYG $OPTIONS $SOURCES_DIR/gengrids_fields.xml > $DEST_DIR/gengrids_fields.tex
$PYG $OPTIONS $SOURCES_DIR/check_cons_interp.m > $DEST_DIR/check_cons_interp.tex
$PYG $OPTIONS $SOURCES_DIR/sumup.m > $DEST_DIR/sumup.tex