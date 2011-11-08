#!/bin/sh

RSTFILE="$1"
TMPFILE=$(mktemp)
echo "$RSTFILE $TMPFILE"

dd bs=4c skip=0 count=1 obs=4c if="$RSTFILE" of="$TMPFILE" 2> /dev/null
FILE_NUMBER=$(hexdump -C "$TMPFILE" | head -1 | sed -e 's/00000000  //' -e 's/ //g' | cut -d'|' -f1)
echo $FILE_NUMBER

dd bs=4c skip=105 count=1 obs=4c if="$RSTFILE" of="$TMPFILE" 2> /dev/null
ID=$(hexdump -C "$TMPFILE" | head -1 | sed -e 's/00000000  //' -e 's/ //g' | cut -d'|' -f1)
echo $ID

HexIntToDec() {
  HEXINT=$1
  DUMMY=$(echo $HEXINT | tr '[:lower:]' '[:upper:]' | sed 's/0*//')
  echo "ibase=16; $DUMMY" | bc
}

HexIntToDec $ID

ANSYS_VERSION=$(dd bs=4c skip=11 count=1 obs=4c if="$RSTFILE" 2> /dev/null)
echo "ANSYS Version: $ANSYS_VERSION"

dd bs=4c skip=12 count=1 obs=4c if="$RSTFILE" of="$TMPFILE" 2> /dev/null
REL_DATE_HEX=$(hexdump -C "$TMPFILE" | head -1 | sed -e 's/00000000  //' -e 's/ //g' | cut -d'|' -f1)
echo -n "ANSYS Release Date: "; HexIntToDec $REL_DATE_HEX

ANSYS_MACHINE_ID=$(dd bs=4c skip=13 count=3 obs=4c if="$RSTFILE" 2> /dev/null)
echo "ANSYS Machine ID: $ANSYS_MACHINE_ID"

dd bs=4c skip=5 count=1 obs=4c if="$RSTFILE" of="$TMPFILE" 2> /dev/null
FILE_DATE_HEX=$(hexdump -C "$TMPFILE" | head -1 | sed -e 's/00000000  //' -e 's/ //g' | cut -d'|' -f1)
echo -n "File Date: "; HexIntToDec $FILE_DATE_HEX

dd bs=4c skip=4 count=1 obs=4c if="$RSTFILE" of="$TMPFILE" 2> /dev/null
FILE_TIME_HEX=$(hexdump -C "$TMPFILE" | head -1 | sed -e 's/00000000  //' -e 's/ //g' | cut -d'|' -f1)
echo -n "File Time: "; HexIntToDec $FILE_TIME_HEX

#dd bs=4c skip=6 count=1 obs=4c if="$RSTFILE" of="$TMPFILE" 2> /dev/null
#UNITS_HEX=$(hexdump -C "$TMPFILE" | head -1 | sed -e 's/00000000  //' -e 's/ //g' | cut -d'|' -f1)
#echo -n "Units: "; HexIntToDec $UNITS_HEX

ANSYS_JOB_NAME=$(dd bs=4c skip=16 count=2 obs=4c if="$RSTFILE" 2> /dev/null)
echo "ANSYS Jobname: $ANSYS_JOB_NAME"

ANSYS_PRODUCT_NAME=$(dd bs=4c skip=18 count=2 obs=4c if="$RSTFILE" 2> /dev/null)
echo "ANSYS Product Name: $ANSYS_PRODUCT_NAME"

ANSYS_SPECIAL_VERSION_LABEL=$(dd bs=4c skip=20 count=1 obs=4c if="$RSTFILE" 2> /dev/null)
echo "ANSYS Special Version Label: $ANSYS_SPECIAL_VERSION_LABEL"

TITLE=$(dd bs=4c skip=41 count=20 obs=4c if="$RSTFILE" 2> /dev/null)
echo "Title: $TITLE"

SUBTITLE=$(dd bs=4c skip=62 count=20 obs=4c if="$RSTFILE" 2> /dev/null)
echo "First Subtitle: $SUBTITLE"

rm -f "$TMPFILE"
