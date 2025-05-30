# 3D Acoustic Channel harmonic simulation using tree stepping (restart aborted simulation)

This example tests the restart functionality in combination with tree stepping for harmonic simulations. The restart is performed by using 

`cfs -r -p input.xml jobname`. 

This way, CFS looks for an result file in the 'results_hdf5' folder and continues the simulation based on the 'input.xml'. Start and End Frequencies have to be the same as in the original simulation. NumFreqs can be increased based on the 'maxTreeLevel' $n$ of the original simulation (default is 10) up to 

$`2^n+1 \overset{n=10}{=} 1025`$

with the default parameter.