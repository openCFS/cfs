BiLinear Forms ([back to main page](/source/CFS_Library_Documentation.md))
==============

A bilinear form, also called integrator, are the terms on the left hand side of a weak formulation; examples are


![](/share/doc/developer/pages/pics/BiLinearForms.png)


Link to [BiLinearForm-doxygen](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1BiLinearForm.html); 

>

The computation of an element matrix according to a bilinear form can be formulated generically by introducing **operators** and **coefficient functions**


![](/share/doc/developer/pages/pics/BBoperator11.png)


Thereby, a bilinear form is defined by plugging operators and coefficient functions to the generic integrators. Operators are classes, which return suiteable matrices according to the shape functions (FE basis functions) defined by the space.


