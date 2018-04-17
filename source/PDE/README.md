PDE    ([back to main page](/source/CFS_Library_Documentation.md))
===

The PDE class is responsible to connect numerics with physics. Please not that for each PDE we generate an own class being derived from **SinglePDE** and it is stored in an own .hh and .cc file. Thereby, we have the following hierarchy:

* BasePDE: Abstract base class
* StdPDE: Base class for all single-field and direct-coupled field problems
* SinglePDE: Base class for all kinds of single field problems (PDEs)

The hierarchy for the PDE for electrostatics is as follows:
```mermaid
graph LR
    A[ElecPDE] --> B[SinglePDE]
    B --> C[StdPDE]
    C --> D[BasePDE]
```

<!--![](/share/doc/developer/pages/pics/PDEsHierarchy.png)-->

Each PDE class contains at least the following methods

* CreateFeSpaces: FeSpace for the unknown physical quantity, see, e.g., [AcousticPDE](/source/PDE/AcousticPDE.cc#L89)
* DefineIntegrators: Bilinear forms (integrators),  see, e.g., [ElecPDE](/source/PDE/ElecPDE.cc#L142)
* DefineNcIntegrators: defines the handling of non-conforming interfaces, see, e.g., [HeatPDE](/source/PDE/HeatPDE.cc#L655)
* DefineSurfaceIntegrators: defines surface bilinear forms (e.g., absorbing boundary conditions), see, e.g., [AcousticPDE](/source/PDE/AcousticPDE.cc#L994)
* DefineRhsLoadIntegrators: right hand side loads (e.g., pressure load, surface traction), see, e.g., [MechPDE](/source/PDE/MechPDE.cc#L1263)
* InitTimeStepping: Time-Stepping scheme, see, e.g., [MechPDE](/source/PDE/MechPDE.cc#L2131)
* DefineSolveStep: defines the solve step for the analysis, see, e.g., [HeatPDE](/source/PDE/HeatPDE.cc#L920)
* DefinePrimaryResults: Available primary results (physical quantity), see, e.g., [MechPDE-Source](/source/PDE/MechPDE.cc#L2141)
* DefinePostProcResults: define post-processing results, see, e.g., [MechPDE](/source/PDE/MechPDE.cc#L2183)

>

Example link to [AcousticPDE-doxygen](https://cfs.mdmt.tuwien.ac.at/docu/doxygen/html/classCoupledField_1_1AcousticPDE.html); for the source see [AcousticPDE-Source](/source/PDE/AcousticPDE.cc#L55)
