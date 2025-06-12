# Testing the implemented MaxwellSlipBC based on the vectorial Lagrange multiplier for moving domains

##  Description

With this test the implementation of the MaxwellSlipBC based on the vectorial Lagrange multiplier for the LinFlow PDE on moving domains gets tested. This test consists of a channel which is excited with two timely seperated velocity pulses. In between those pulses the domain changes leading to different stress gradients at the boundaries and therefore different slip velocities.
The channel itself is modeled as one half of a symmetric channel with a hole in a small intermediate chamber. The walls of this chamber move in positive x-direction after the first pulse, leading to a smaller channel for the fluid. Hence, the normal stress on the boundary where the MaxwellSlipBC lives is increased, which leads to a higher slip velocity
Due to the increased viscosity as well as the adapted scaling factor, the scaled strain is high enough to cause some slip for both cases. 

## File structure

MaxwellSlipBC2D.* ... files used for the computation/evaluation of the testcase. A Paraview state file is added to compare solution (slip velocity) at to time steps (t = 125e-6 s: v_slip = 0.22 m/s; t = 900e-6 s: v_slip = 0.64 m/s)
For this test we have to use the .cdb file for the mesh input because the .h5ref file causes problems.


