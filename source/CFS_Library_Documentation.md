Many of the classes in CFS++ are grouped into modules. These modules define the building blocks of any finite element program: 

* Computational Domain (CoefFunction, ElemMapping, Mesh, Results, ...) 
* [PDEs](/source/PDE/README.md) (StdPDE, SinglePDE, AcousticPDE, ElecPDE, MechPDE, ...)
* CoupledPDEs (DirectCoupledPDEs, IterCoupledPDEs, ...)
* Driver (MultiSequence, Static, Transient, Harmonic, EigenFrequency, ...)
* FeBasis (FeSpace for H1, Hcurl, Hdiv, L2; FeFunction, ...)
* DataInOut (SimInput, SimOutput, SimState, ParamHandling, ...)
* Forms (BiLinForms, LinForms, Operators, ...)
* Materials (AcousticMaterial, MechaniMaterial, MagneticMaterial, ...)
* Algebraic Solver OLAS (AlgSys, Solver, PreCond, Graph, ...)
* Optimization (Condition, Ersazmaterial, LevelSet, ...)
* Utils (Approxdata, AutoDiff, ...)

>

The directory structure of the source looks as follows (e.g., **eclips project**)

![](/share/doc/developer/pages/pics/Src_Directory.png)