#include <fstream>
#include <iostream>
#include <string>

#include "harmonicAdapt.hh"

#include <PDE/basepde.hh>
#include "basedriver.hh"


namespace CoupledField
{

  HarmonicAdaptSpaceDriver :: HarmonicAdaptSpaceDriver(Domain * adomain)
    : HarmonicDriver(adomain)
  {
#ifdef TRACE
    (*trace) << "entering HarmonicAdaptSpaceDriver::HarmonicDriver" << std::endl;
#endif
  
    flags->adaptSpace_ = true;

  }

  HarmonicAdaptSpaceDriver :: ~HarmonicAdaptSpaceDriver()
  {
#ifdef TRACE
    (*trace) << "entering HarmonicAdaptSpaceDriver::~HarmonicAdaptSpaceDriver" << std::endl;
#endif

  }

  void HarmonicAdaptSpaceDriver :: SolveProblem()
  {
#ifdef TRACE
    (*trace) << "entering HarmonicAdaptSpaceDriver::SolveProblem" << std::endl;
#endif

    Integer      level=0;
    Integer      pdenumber  = 0;                      // we solve only 1 equation
    bool      printMeshes;                         // output intermediate meshes or not
    Integer      numrepeat = 0;                       // counter of adaptive steps
    Integer      maxnumrepeat = 0;                    // maximum number of loops for rifenement
    BasePDE*     ptPDE;                               // pointer to PDE
  
    if (!conf->ifget("maxnumrepeat",maxnumrepeat,"SpaceAdaptivity"))
      maxnumrepeat=10;

    printMeshes = printMeshesOrNot();

    if (printMeshes)
      {
        ptMeshes_->Init(ptdomain_->GetGrid());
        PrintSeqMeshes();
      }
      
    ptPDE = ptdomain_->GetPDE(pdenumber);
    ptPDE->SolveStepHarmonic(level);   
  
    while(ptPDE->TestError(level) && numrepeat<maxnumrepeat)
      {
        ptPDE->RefineMesh();
      
        level++;                                       // new level of grid
      
        if (printMeshes)
          {
            ptMeshes_->Init(ptdomain_->GetGrid());
            ptPDE->WriteErrorInfo(ptMeshes_);         // info for old mesh
            PrintSeqMeshes();                         // print new mesh
          }

        ptdomain_->Update(level);

        ptPDE->SolveStepHarmonic(level);

        numrepeat++;
      }

    ptPDE->PostProcess(level);
    ptdomain_->PrintGrid(level);
    ptPDE->WriteResultsInFile();

  }

} // end of namespace
