# mech-smooth-linFlow coupling

## Description

Test to check the Mech-Smooth coupling via contact.

This test consists of a metal plate connected with a volume filled with air. The plate is excited with a pressure BC.
First the mechanicPDE is calculated using the pressure BC. Here, Dirichlet BCs are used to restrict movement only to the x-axis. The quantity 'mechDisplacement' is used to couple the mechanic domain to the smoothPDE at the interface. For the smoothPDE a contact law is described via the mathParser. This displacement dependent contact force scales linearly with the displacement and starts to act when the gap between the two compressed surfaces is below a certain threshold. Since the system is undamped, the plate simply bounces back. In order to achieve a convergent solution, we define an iterative coupling with corresponding stopping criteria within the couplingList.

## File structure

CompressedChamberContact2D.* ... files used for the computation of the testcase.



