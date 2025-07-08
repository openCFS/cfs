# Test case for CFSDat time derivative

## WARNING: time derivative is identical with analytical solution but the result is shifted one time step!
Due to the time shift, the test case is not checked in yet.

## Testcase:
- analytical field is generated with *cfs* in *srcData*
- time derivative is computed with *cfsdat*
- data for comparison is extracted with *paraview* state file *eval.pvsm* and stored in *dpdt0.0.csv*
- comparison is done with *eval.py*
