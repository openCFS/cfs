# Testing the implemented noPenetration BC based on the Lagrange multiplier

## Description

With this test the implementation of the noPenetration BC based on the scalar Lagrange multiplier for the LinFlow PDE gets tested. This test consists of a channel which is excited with a velocity. The velocity normal to the side boundaries is zero. 

## File structure

PlaneWaveNoPenetrationScalarLM.* ... files used for the computation/evaluation of the testcase. The channel is rotated and the scalarVelocityConstraint is used as BC. A Paraview state file is added to compare the reference with the results of the rotated channel computation.
not_rotated ... here all files for the computation of a reference solution can be found. The channel is not rotated and the velocity normal to the boundary is zero (noPenetration). 

