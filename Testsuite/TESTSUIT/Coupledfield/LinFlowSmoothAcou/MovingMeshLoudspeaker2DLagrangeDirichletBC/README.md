## Analysis of a channel with updated geometry

###  Description

This test consists of a channel with three sections. A LinFlow section is coupled with an acoustic one which adjoins the PML (perfectly matched layer). The geometry is updated via smoothPDE. The velocity and gridvelocity are assigned to the first part of the Linflow section as boundary condition. This results in an oscillation and wave in the channel. Therefore, the geometry of the testcase stands in analogy to loudspeakers. 
A characteristic of this testcase is that the maximum amplitude of the oscillation is greater than one grid box. 

This test is similar to the original one (MovingMeshLoudspeaker2D) - the only difference is, that we use Lagrange multiplier based BCs (noPenetration BCs on all horizontal surfaces (including moving ones) as well as a noSlip BC for the excitation). With this testcase we can test the BCs on moving domains.



