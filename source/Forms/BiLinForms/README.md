BiLinear Forms ([back to main page](/source/CFS_Library_Documentation.md))
==============

A bilinear form, also called integrator, are the terms on the left hand side of a weak formulation; examples are
* Source [BBInt](/source/Forms/BiLinForms/BBInt.hh) and [doygen documentation](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1BBInt.html)
* Source [BDBInt](/source/BiLinForms/BDBInt.hh) and [doxygen documentation](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1BDBInt.html)
* Source [ABInt](/source/BiLinForms/ABInt.hh) and [doxygen documentation](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1ADBInt.html)

<!--\[
\textrm{BBInt:}\ \int\limits_\Omega \, \alpha \nabla \varphi \cdot \nabla u\, \textrm{d}\Omega
\]

\[
\textrm{BDBInt:}\ \int\limits_\Omega \, \left( B \vec u'\right)^t\, D\, \left( B \vec u\right)^t\, \textrm{d}\Omega
\]

\[
\textrm{ABInt:}\ \int\limits_\Omega \, \alpha\,  \vec \psi\, \cdot \nabla u\,  \textrm{d}\Omega
\]
-->

![](/share/doc/developer/pages/pics/BiLinearForms.png)

Link to [BiLinearForm-doxygen](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1BiLinearForm.html);


The computation of an element matrix according to a bilinear form can be formulated generically by introducing **operators** and **coefficient functions**

<!--\[
\int\limits_\Omega \, \underbrace{\alpha}_{\textrm{CoefFunction}}\ \  \underbrace{\nabla \varphi}_{\textrm{BOperator}}\ \cdot \ \underbrace{\nabla u}_{\textrm{BOperator}}\ \ \textrm{d}\Omega
\]
-->

![](/share/doc/developer/pages/pics/BBoperator11.png)


Thereby, a bilinear form is defined by plugging operators (e.g., a [GradientOperator](source/Forms/Operators/GradientOperator.hh) and coefficient functions (e.g., [constnat material parameter](/source/Domain/CoefFunctions/CoefFunctionConst.hh) to the generic integrators. Operators are classes, which return suitable matrices according to the shape functions (FE basis functions) defined by the space.
