Tour through CFS++  ([back to main page](/source/CFS_Library_Documentation.md))
===

We perform a simple static mechanical computation and explain, which objects are created at which instant of time during the computation! For this purpose, we model a plate in 3D with one computational domain named *compDomain* and a surface *surfLoad*, where we apply a pressure load. The plate is fully supported by fixing all degrees of freedom at the boundary to zero.

```
<?xml version="1.0"?>
<cfsSimulation xmlns="http://www.cfs++.org" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.cfs++.org http://cfs-doc.mdmt.tuwien.ac.at/xml/CFS-Simulation/CFS.xsd">

  <fileFormats>
    <output>
      <hdf5/>
    </output>
    <materialData file="mat.xml" format="xml"/>
  </fileFormats>

  <domain geometryType="3d">
    <regionList>
      <region name="compDomain"/>
    </regionList>
    <surfRegionList>
      <surfRegion name="surfLoad"/>
    </surfRegionList>
  </domain>

  <sequenceStep>
    <analysis>
      <static/>
    </analysis>

    <pdeList>
      <mechanic subType="3d">
        <regionList>
          <region name="compDomain"/>
        </regionList>

        <bcsAndLoads>
          <fix name="fix-all">
            <comp dof="x"/>
            <comp dof="y"/>
            <comp dof="z"/>
          </fix>
          <pressure name="SurfLoad" value="1"/>
        </bcsAndLoads>
        <storeResults>
          <nodeResult type="mechDisplacement">
            <allRegions/>
          </nodeResult>
      </storeResults>
      </mechanic>
    </pdeList>
  </sequenceStep>

</cfsSimulation>
```
In the [main program](/source/main/CFS.cc) the following lines are executed:
```
int main(int argc, const char **argv)
{
  CFS cfs(argc, argv);
  int ret = cfs.Run();
  return ret;
}
```
So, we see that the object **CFS** is initiated and then the method **Run()** of this object is called, which performs the following main steps:
```
ReadXMLFile();
domain = new Domain( gridInputs, resultHandler, materialHandler, simState, paramNode_, infoNode );
domain->CreateGrid();
SolveProblem();
```
This means that the simulation xml-file is parsed and then the object **domain** is initiated with the following arguments:
* names of mesh filesile
* pointer to [result handler](/source/DataInOut/ResultHandler.hh)
* pointer to [material handler](/source/DataInOut/ParamHandling/MaterialHandler.hh)
* pointer to [simulatioon state handler](/source/DataInOut/SimState.hh)
* pointer to [simulation parameter node (simulation xml-file)](/source/DataInOut/ParamHandling/ParamNode.hh)
* pointer to [material parameter node (material xml-file)](/source/DataInOut/ParamHandling/ParamNode.hh)

In a next step, the method **CreateGrid()** of class **Domain** is called so that the data from the grid input files are read (we can handle more than one grid input file). Furthermore, an object of class [GridCFS](/source/Domain/Mesh/GridCFS/GridCFS.hh) is initiated. Finally, the method **SolveProblem()** of the class **CFS** is called, which performs the following calls
```
domain->PostInit();
domain->SolveProblem();
```
The first command checks which type of analysis (Static, Transient, Harmonic, EigenFrequency) is performed, or if we do a optimization problem. In addition, we check for a *multisequence* analysis (e.g., first sequence is a static one and second is a transient one, etc.). The classes behind are the [driver-classes](/source/Driver/BaseDriver.hh), and in our case, a [static driver](/source/Driver/StaticDriver.hh) is used. The last main step in *domain->PostInit()* is to call *driver->Init(restart)*, which calls **InitializePDEs()** a method being implemented in [SingleDriver](/source/Driver/SingleDriver.cc) and performing the following steps:
```
domain_->CreatePDEs( 1, info_->GetParent() );
ptPDE_ = domain_->GetBasePDE();
domain_->InitPDEs( 1 );
```
So, we see that **domain** is a quite powerful class, which first initiates all single PDEs being defined in the simulation xml-file (in our case it is just *mechanic*) and then triggers the PDEs to perform their tasks. The driver (in our case *StaticDriver*) also gets  a pointer to the PDE, so that it can directly trigger actions to be done in the PDE (see later).

The second call in the method **SolveProblem()** of the class **CFS** (see above) is *domain->SolveProblem()*, which executes in the class **Domain** the call *driver->SolveProblem()*. In our case, the method *Solveproblem()* of class **StaticDriver** is called, which performs the following three main steps implemented in the class **StdSolveStep** ([see StdSolveStep.cc](/source/Driver/SolveSteps/StdSolveStep.cc))
```
ptPDE_->GetSolveStep()->PreStepStatic();
ptPDE_->GetSolveStep()->SolveStepStatic();
ptPDE_->GetSolveStep()->PostStepStatic();
```
The fist command just triggers the method *algsys_->InitRHS()*, in our case just setting the vector of right hand side to zero. The second command (we have a linear problem, so *StepStaticLin()* is called) triggers all steps so that the stiffness matrix (which is in our case equal to the system, matrix) and right hand side is computed. Furthermore, the algebraic system is solved and the solution stored. Here are the main steps being executed ([see StdSolveStep.cc](/source/Driver/SolveSteps/StdSolveStep.cc))
```
assemble_->AssembleMatrices();
assemble_->AssembleLinRHS();
PDE_.SetRhsValues();
PDE_.SetBCs();
algsys_->ConstructEffectiveMatrix( NO_FCT_ID,matrix_factor_[NO_FCT_ID] );
algsys_->BuildInDirichlet();
algsys_->Solve();
algsys_->GetSolutionVal(solVec_);
```
In our case, *PostStepStatic()* is doing noting and we are back in *StaticDriver::SolveProblem()*. Here, the method *StoreResults* of class **StaticDriver** takes care that all results are written to the output files by executing
```
ptPDE_->WriteResultsInFile(stepNum, step_val);
ptPDE_->WriteGeneralPDEdefines();
handler_->FinishStep();
simState_->WriteStep(stepNum, step_val );
```

Finally, we show by the help of a sequence diagram
```mermaid
sequenceDiagram
User -> A: DoWork
activate A
A -> B: << createRequest >>
activate B
B -> C: DoWork
activate C
C --> B: WorkDone
destroy C
B --> A: RequestCreated
deactivate B
A -> User: Done
deactivate A
```

```mermaid
    participant CFS
    participant Domain
    participant GridCFS
    CFS-->>Domain: 1:Create
    CFS->>Domain: 2:call CreateGrid
    Domain-->>GridCFS: 3:Create
```
