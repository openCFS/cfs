#!/bin/sh

# computes the coding score 
# excepts a ctest build output as input, i.e. from `ctest -D ExperimentalBuild --track Gitlab -V > log.out`

log=$1

echo "end of $1"
echo "-----------------------------------------"
tail -n 5 $1
echo "-----------------------------------------"

errors=$( grep '[0-9]\+ \(or more \)\?Compiler errors' $1 | sed 's|[ ]\+\([0-9]\+\).\+|\1|' )
warnings=$( grep '[0-9]\+ \(or more \)\?Compiler warnings' $1 | sed 's|[ ]\+\([0-9]\+\).\+|\1|' )
insertions=$( (git diff --stat HEAD~1 HEAD | tail -n 1 | grep insertion || echo 'a, 0 insertions') | sed 's|.\+, \([0-9]\+\) insertion.\+|\1|')
deletions=$( (git diff --stat HEAD~1 HEAD | tail -n 1 | grep deletion || echo 'a, 0 deletions') | sed 's|.\+, \([0-9]\+\) deletion.\+|\1|')

echo "$errors errors, $warnings warnings, $insertions insertions, $deletions deletions"
# compute and print the score:
# 0 for more errors and warnings than new lines
# (1-($errors+$warnings)/($insertions+$deletions))*100 otherwise
awk "BEGIN {printf \"coding_score=%.2f\n\", ( ($errors+$warnings)>=($insertions+$deletions)) ? 0.0 : (1-($errors+$warnings)/($insertions+$deletions))*100 }"

exit 0
