#!/bin/sh

freecores=1

ncores=$(nproc)
ncores=$((ncores-freecores))
ncores=6

if [ "$ncores" -lt "0" ]; then
  echo "Too much free cores requested. Please lower number of cores to be left free."
  exit
fi

filename=${1##*/} # Split string at / and take last part
filesize=$(stat --printf="%s" $1)

blocksize=$(($filesize/$ncores))
blocksize=$(($blocksize>0?$blocksize:1))

logfile=$filename.log
rm $logfile
rm catalogues/detailed_stats_$filename

parallel --will-cite --pipepart -a $1 --block $blocksize --joblog $logfile --jobs $ncores --resume-failed --header '.*\n' "cat > {#}; ./callMatlab.sh {#}"

cat catalogues/detailed_stats_[123456789]* >> catalogues/detailed_stats_$filename

rm [123456789]*
rm catalogues/detailed_stats_[123456789]*

#matlab -nodesktop -nodisplay -nosplash -r "try;writeHeader('catalogues/detailed_stats_$filename');catch ME;disp(ME.message);end;exit;">$filename.out 2>$filename.err

