#!/bin/sh

# Script to generate a graph showing the internal dependencies of all
# targets of CFS++.
# This script uses the built-in graphviz support of CMake. Modify the
# settings in CMakeGraphVizOptions.cmake to exclude targets from the
# graph. For further information see the documentation therein.

usage() {
    echo "./gendepgraph.sh [options]"
    echo "Generate a graph showing the internal dependencies of all targets."
    echo "Possible options:"
    echo "  -f EXT [--format=EXT]  : Format of graph file (e.g. fig|gif|imap|pdf|png|ps|svg). See http://www.graphviz.org/content/output-formats for more extensions. [default: pdf]"
    echo "  -l LEN [--length=LEN]  : Length for unflattening the graph. The minimum length of leaf edges is staggered between 1 and LEN. [default: 1]"
    echo "  -c [--colored]     : Edges of graph are colored."
    echo ""
}

# Initialize
FORMAT=pdf
UNFLAT_LENGTH=1
COLORED=0

# Command line parsing
# Note that we use `"$@"' to let each command-line parameter expand to a
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
OPTIND=1
TEMP=`getopt -o f:l:c --long format:,length:,colored \
     -n 'gendepgraph.sh' -- "$@"`

if [ $? != 0 ] ; then usage ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Parse options
while true ; do
  case "$1" in
    -h|--help|?) usage ; shift ;;
    -f|--format) FORMAT="$2" ; shift 2 ;;
    -l|--length) UNFLAT_LENGTH="$2" ; shift 2 ;;
    -c|--colored) COLORED=1 ; shift ;;
    --) shift ; break ;;
    *) echo "Internal error!"; exit 1 ;;
  esac
done

unset CDPATH
CURRENTDIR="$(pwd)"

# Get directory where this script is stored
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Copy settings file to release folder
if [ ! -d "$SCRIPTDIR/../../release" ]; then
  mkdir "$SCRIPTDIR/../../release"
fi
cp "$SCRIPTDIR/CMakeGraphVizOptions.cmake" "$SCRIPTDIR/../../release"
cd "$SCRIPTDIR/../../"

# Get name of project
PROJECTNAME="${PWD##*/}"
FILENAME="$PROJECTNAME""_dependencies"

# Call CMAKE with graphviz option to generate dot files
cd "$SCRIPTDIR/../../release"
cmake .. --graphviz="$SCRIPTDIR/../../dependencygraph/$FILENAME"
cd "$SCRIPTDIR/../../dependencygraph"

# Colorize nodes and edges using HSV
if [ "$COLORED" = "1" ]; then
  gawk -i inplace '
  BEGIN {
    hue = 0;
    sat = 1;
    value = 1;
    satdrain = 0;
    colormap[0] = 0;
    node_regex = "\"node([0-9]*)\" \\[ label=(.*)\\];";
    edge_regex = "\"node([0-9]*)\" -> \"node([0-9]*)\"";
  }
  {
    if ( match( $0, node_regex, currentnodeid ) ) {

      # Colorize current node
      gsub(node_regex, substr(currentnodeid[0],1,length(currentnodeid[0])-2)" color=\""hue/360" "sat" "value"\" ]");
      
      # Set colormap
      colormap[currentnodeid[1],0] = hue;
      colormap[currentnodeid[1],1] = sat;
      colormap[currentnodeid[1],2] = value;
      
      hue += 6.1803398875; # golden ratio assures we never get a hue value modulo 360 twice!
      if (hue >= 360) {
        sat -= .2;
        value -= .2;
        hue -= 360;
      };
      
      # Do not use too bright or too dark colors
      if (sat < .15 || value < .3) {
        satdrain += .2;
        sat = 1 - satdrain;
        value = 1;
      };
      
      if (sat < .15)
        print "Warning: Not enough colors available." > "/dev/stderr";
    }
    
    # Colorize current edge in the same color as the start node
    match($0, edge_regex, currentedge);
    gsub(edge_regex, "& [ color=\""colormap[currentedge[1],0]/360" "colormap[currentedge[1],1]" "colormap[currentedge[1],2]"\" ]");
    print;
  }
  ' "$FILENAME"
fi

# Unflatten graph and write to file
unflatten -f -l$UNFLAT_LENGTH "$FILENAME" | dot -T$FORMAT -x -Gratio=0.70710678118 -o"$FILENAME".$FORMAT

cd "$CURRENTDIR"



## To color a node use
#sed -i 's/label="NODENAME"/label="NODENAME" color="COLOR"/' "$1"

## Node colors can flow through the graph by
#dot "$FILENAME" | gvcolor | unflatten -f -l$UNFLAT_LENGTH | dot -T$FORMAT -x -Gratio=.5 -o"$FILENAME".$FORMAT

## Scale the pdf to a poster (with overlap for printing on A4 sheets)
# Scale to A3
#pdfposter -O 2 -s 0.258 -x -1 "$FILENAME".$FORMAT "$FILENAME"_poster.$FORMAT
# Scale to A2
#pdfposter -O 2 -s 0.365 -x -1 "$FILENAME".$FORMAT "$FILENAME"_poster.$FORMAT
