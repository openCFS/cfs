# mech-smooth-linFlow coupling

## Description

Test to check the full (iterative) LinFlow-Smooth-Mech coupling.

This test consists of a metal plate connected with a volume filled with air. The plate is excited with a pressure BC.
First the mechanicPDE is calculated using the pressure BC. Here, Dirichlet BCs are used to restrict movement only to the x-axis. The quantity 'mechDisplacement' is used to couple the mechanic domain to the smoothPDE at the interface. Afterwards, the fluidMechLinPDE is calculated and coupled via the 'smoothVelocity'. As a last step the fluidMechLinPDE is recoupled to the mechanicPDE. This is done by applying a 'fluidMechNormalSurfaceStress' as traction in the mechanicPDE which counteracts the pressure with which the plate is excited. In order to achieve a convergent solution, we define an iterative coupling with corresponding stopping criteria within the couplingList.

## File structure

CompressedChamber2D.* ... files used for the computation of the testcase. A matlab file (eval_hist.m) as well as a state file for paraview (CompressedChamber2D.pvsm) are also included in the folder.



