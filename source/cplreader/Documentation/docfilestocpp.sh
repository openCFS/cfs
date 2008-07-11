#!/bin/sh -x

# Generate ASCII hexdumps of PNG and HTML files for creation of C++ source code.
OUTFILE="$1"

echo "%%%%%%%%%%%% OUTFILE: $OUTFILE" > $HOME/out.txt

rm -f $OUTFILE

echo "#include \"GetDocumentation.hh\"" >> $OUTFILE
echo "" >> $OUTFILE
echo "namespace CoupledField {" >> $OUTFILE
echo "" >> $OUTFILE


echo "  void GetDocumentation(std::map<std::string, std::vector<unsigned char> >& doc) {" >> $OUTFILE
echo "" >> $OUTFILE

DOCFILES=`ls *.png`
DOCFILES="$DOCFILES `ls *.html`"
COUNTER=1

for DOCFILE in $DOCFILES; do
    echo "    static unsigned char doc_file$COUNTER[] =" >> $OUTFILE
    echo "    {"  >> $OUTFILE
    hexdump -C $DOCFILE | awk '{fields=NF-2; if(fields > 16) fields=16; printf("      "); for (i=0;i<fields;i++) { printf("0x%s", $(i+2)); if(fields != 16 && i==fields-1) printf("\n    };"); else printf(", "); } print ""; }' >> $OUTFILE
    
    echo "    std::copy( doc_file$COUNTER, doc_file$COUNTER + sizeof(doc_file$COUNTER), std::back_inserter(doc[\"$DOCFILE\"]) );" >> $OUTFILE
    echo "" >> $OUTFILE

    COUNTER=$((COUNTER + 1))
done

echo "  }" >> $OUTFILE

echo "} // end of namespace" >> $OUTFILE
