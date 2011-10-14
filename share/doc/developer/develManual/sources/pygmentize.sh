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
Pygmentize qs_obtain_sources.sh qs_obtain_sources.tex
Pygmentize git_checkout_cfs.sh git_checkout_cfs.tex
Pygmentize qs_configure_build.sh qs_configure_build.tex
Pygmentize kreuzinger.sh kreuzinger.tex
Pygmentize platform_defaults_mac.cmake platform_defaults_mac.tex
Pygmentize macports_installed_pckgs.txt mac_ports_installed_pckgs.tex
Pygmentize nacs_win.txt nacs_win.tex
Pygmentize nacs_win_build_howto.txt nacs_win_build_howto.tex
Pygmentize mpcci_reactivate_lse24_lic.txt mpcci_reactivate_lse24_lic.tex
Pygmentize org.eclipse.ui.editors.prefs org.eclipse.ui.editors.prefs.tex
Pygmentize org.eclipse.core.resources.prefs org.eclipse.core.resources.prefs.tex
Pygmentize emacs_settings.sh emacs_settings.tex
Pygmentize qs_bootstrap.sh qs_bootstrap.tex
