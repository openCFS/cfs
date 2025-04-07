# 2D transient analysis of a channel with updated geometry

## Description

This test consists of a channel with three sections. The middle section is predeformed via the smoothPDE in a static step. Afterwards, the static deformation is used to compute the transient solution of the LinFlow-PDE in the deformed geometry. Due the static precomputation the initial velocity caused by the smoothPDE is zero. Hence, this testcase is used to check the geometry update of the integratos used in the LinFlow-PDE but also the correct initialization of the GLM-scheme based on an initial state.

## File structure

ChannelUpdateGeom.* ... files used for the computation of the testcase. Here the geometry gets updated by the smoothPDE where a sine shaped dent is prescribed with a certain dent height.
Predeformed_reference ... here all files for the computation of a reference solution an a predeformed mesh can be found (the .jou is almost the same, but here we specify the dent height for the meshing process).

## Test definition

The velocity computed by the smoothPDE based on the static initial state should be zero but differs slightly from zero due to machine precision. Hence, the relL2diff comparison might lead to large errors. Therefore, only the absL2diff test is activated to avoid the test to more or less randomly fail.
