#!/bin/sh

rm my_log

ncores=$(nproc)

filename=${1##*/}
filesize=$(stat --printf="%s" $1)

blocksize=$(($filesize/$ncores))
blocksize=$(($blocksize>0?$blocksize:1))

parallel --will-cite --pipepart -a $1 --block $blocksize --joblog my_log --jobs 200% --resume-failed --header '.*\n' "cat > {#}; ./callMatlab.sh {#}"

cat catalogues/detailed_stats_[123456789]* >> catalogues/detailed_stats_$filename

rm [123456789]*
rm catalogues/detailed_stats_[123456789]*

matlab -nodesktop -nodisplay -nosplash -r "sortCatalogue('catalogues/detailed_stats_$filename','$1');exit;">$1.out 2>$1.err

