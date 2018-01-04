Finite Element Space and FE  ([back to main page](/source/CFS_Library_Documentation.md))
===

Capabilities of FESpace
-----------------------

* Manage equation numbers for the associated unknown
* Hold and manage associated reference elements
* Know how to incorporate boundary conditions

![](/share/doc/developer/pages/pics/FeSpaceStruct.png)

>

Relation to geometric element
-----------------------------

* CFS++ distiguishes between geometric elements (usually linear and quadratic) as given by the preprocessor and computational element (the basis functions for approximating the unknown provided by the FESpace)
* Possibility to calculate with arbitrary order on any given mesh
* During element matrix computation both „types“ of elements are used

![](/share/doc/developer/pages/pics/GeometricVersCompElement.png)

>

[See doygen documentation](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1FeSpace.html);

FeFunction 
----------

This class represents a physical quantity used in the PDE

* For time dependent problems, also time derivatives are represented as **FeFunctions**
* Each (unknown) !FeFunction is associated to a unique !FeSpace,  which represents e.g. its equation numbers
* This class has methods to give the solution on a specific element or on the complete domain
* A **FeFunction** can also be passed to an integrator (e.g. for non-linear PDEs)

![](/share/doc/developer/pages/pics/FeFunctionExample1.png)


[See doygen documentation](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1FeFunction.html);
