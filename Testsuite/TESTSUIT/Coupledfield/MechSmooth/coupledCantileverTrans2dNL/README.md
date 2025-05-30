# nonlinear mech-smooth coupling

## Description

Test to check subiterations from a nonlinear mechanics problem inside the iterations stemming from an iteratively coupled problem.

This test is similar to the testcase cantileverTrans2dNL for mechanics where a cantilever is excited with a force. Compared to the original test we excite with a force density since currently only nodal forces are implemented and the mesh generation with Trelis only gives us surfaces or lines within CFS. Nevertheless, the input values are adjusted in order to obtain similar results compared to the original test. The deformation itself is large, hence, a nonlinear solution strategy is applied. The only difference is that in this test an air domain is iterativel coupled where the smoothPDE is solved. The setup consists of a simple forward coupling (mechanics->smooth) and should test the functionality of the sub-iterations from the nonlinear problem inside the iterations stemming from the coupling.

## File structure

coupledCantileverTrans2dNL.* ... files used for the computation of the testcase.
Comparison.pvsm ... Paraview state file for the comparison of the reference solutions of the original problem and the coupled one.



