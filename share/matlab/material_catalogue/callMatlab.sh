#!/bin/sh
# This is a helper for the parallel generation of a material catalogue
# (c.f. generateCatalogueParallel.sh).

echo $1

matlab -nodesktop -nodisplay -nosplash -singleCompThread -r "callCatalogueGeneration('$1')">$1.out 2>$1.err
