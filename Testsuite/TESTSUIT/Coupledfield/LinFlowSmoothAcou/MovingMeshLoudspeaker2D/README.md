## Analysis of a channel with updated geometry

###  Description

This test consists of a channel with three sections. A LinFlow section is coupled with an acoustic one which adjoins the PML (perfectly matched layer). The geometry is updated via smoothPDE. The velocity and gridvelocity are assigned to the first part of the Linflow section as boundary condition. This results in an oscillation and wave in the channel. Therefore, the geometry of the testcase stands in analogy to loudspeakers. 
A characteristic of this testcase is that the maximum amplitude of the oscillation is greater than one grid box. 



