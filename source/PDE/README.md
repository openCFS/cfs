PDE
===

The PDE class is resonsible to connect numerics with physics. It defines
* CreateFeSpaces: FeSpace for the unknown physical quantity
* DefineIntegrators: Bilinear forms (integrators)
* DefineNcIntegrators: defines the handling of non-conforming interfaces
* DefineSurfaceIntegrators: defines surface bilinear forms (e.g., absorbing boundary conditions)                               
* DefineRhsLoadIntegrators: right hand side loads (e.g., pressure load, surface traction)
* InitTimeStepping: Time-Stepping scheme
* DefineSolveStep: defines the solve step for the analysis
* DefinePrimaryResults: Available primary results (physical quantity, e.g., mechanical displacement)
* DefinePostProcResults: define post-processing results (e.g., mechanical stress)


Please not that for each PDE we generate an own class being derived from {{{SinglePDE}}} and it has to be stored in an own .hh and .cc file. Thereby, we have the following hierarchy:
* BasePDE: Abstract base class 
* StdPDE: Base class for all single-field and direct-coupled field problems
* SinglePDE: Base class for all kinds of single field problems (PDEs)
