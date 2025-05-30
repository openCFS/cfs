# mech-smooth-linFlow-electrostatic coupling

## Description

Test to check the full (iterative) LinFlow-Smooth-Mech-Electrostatic coupling.

This test consists of a metal plate connected with a volume filled with air. The plate is excited with a pressure BC as well as an electric force that is changing due to the movement of the domain.
First, the electrostaticPDE is used to calculate the additional force acting on the plate. Then, the mechanicPDE is calculated using the pressure BC as well as the electrostatic force acting on the plate. Here, Dirichlet BCs are used to restrict movement only to the x-axis. The quantity 'mechDisplacement' is used to couple the mechanic domain to the smoothPDE at the interface. Afterwards, the fluidMechLinPDE is calculated and coupled via the 'smoothVelocity'. As a last step the fluidMechLinPDE is recoupled to the mechanicPDE. This is done by applying a 'fluidMechNormalSurfaceStress' as traction in the mechanicPDE which counteracts the forces with which the plate is excited. In order to achieve a convergent solution, we define an iterative coupling with corresponding stopping criteria within the couplingList.

## File structure

ElecForceCompressedChamber2D.* ... files used for the computation of the testcase. A matlab file (eval_hist.m) for the evaluation of the equlibrium also included in the folder.



