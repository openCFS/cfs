Finite Element Space  ([back to main page](/source/CFS_Library_Documentation.md))
===

Capabilities of FESpace
==

* Manage equation numbers for the associated unknown
* Hold and manage associated reference elements
* Know how to incorporate boundary conditions

![](/share/doc/developer/pages/pics/FeSpaceStruct.png)


The PDE class is resonsible 

* BasePDE: Abstract base class 

![](/share/doc/developer/pages/pics/PDEsHierarchy.png)

>

Each PDE class contains at least the following methods

* CreateFeSpaces: FeSpace for the unknown physical quantity, see, e.g., [AcousticPDE-Source](/source/PDE/AcousticPDE.cc#L89)

>

Example link to [AcousticPDE](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1AcousticPDE.html); for the source see [AcousticPDE-Source](/source/PDE/AcousticPDE.cc#L55)
