# Testing the implemented MaxwellSlipBC based on the vectorial Lagrange multiplier

## Description

With this test the implementation of the MaxwellSlipBC based on the vectorial Lagrange multiplier for the LinFlow PDE gets tested. This test consists of a channel which is excited with a velocity. After a short inflow region, where a simple noPenetration BC is used in order to let the flow develop, a MaxwellSlipBC is applied on the right side. Due to the increased viscosity as well as the adapted scaling factor, the scaled strain is high enough to cause some slip. We use a finer mesh (reference folder) where we numerically check if the tangential slip velocity is calculated correctly. In this case the tangential velocity (uy) has to be equal to -(duy/dx+dux/dy)*C1*meanFreePath. 
In three sequence steps the formulation of the MaxwellSlipBC differs ('hollowIncompressibleStress', 'incompressibleStress' and 'compressibleStress').

## File structure

MaxwellSlipBC2D.* ... files used for the computation/evaluation of the testcase. A Paraview state file is added to compare the h5ref-file with the finely meshed reference solution. Another Paraview state file file compares the three sequence Steps.
reference ... here all files for the computation of a reference solution can be found. The channel is finely meshed so that the element result (normal derivative) used for the comparison/evaluation is close to the actual value at the boundary.

