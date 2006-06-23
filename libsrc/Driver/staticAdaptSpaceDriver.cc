#include <fstream>
#include <iostream>
#include <string>

#include "staticAdaptSpaceDriver.hh"

#include "DataInOut/GMV/outGMV.hh"
#include <PDE/basepde.hh>
#include "basedriver.hh"


namespace CoupledField
{

  StaticAdaptSpaceDriver :: StaticAdaptSpaceDriver(Domain * adomain)
    : StaticDriver(adomain)
  {
    ENTER_FCN( "StaticAdaptSpaceDriver::StaticDriver" );
  
    flags->adaptSpace_ = true;

  }

  StaticAdaptSpaceDriver :: ~StaticAdaptSpaceDriver()
  {
    ENTER_FCN( "StaticAdaptSpaceDriver::~StaticAdaptSpaceDriver" );

  }

  void StaticAdaptSpaceDriver :: SolveProblem()
  {
    ENTER_FCN( "StaticAdaptSpaceDriver::SolveProblem" );

    Integer level=0;
    Integer pdenumber  = 0;
    Integer numrepeat = 0;

    Integer maxnumrepeat;
    if (!conf->ifget("maxnumrepeat",maxnumrepeat,"SpaceAdaptivity"))
      maxnumrepeat=10;

    // print all stages of mesh 
    bool meshesInfo=printMeshesOrNot();
  
    if (meshesInfo)
      { // print info about mesh in mesh-info-file
        ptMeshes_->Init(ptdomain_->GetGrid());
        PrintSeqMeshes();
      }
   
    if (ptdomain_->GetNumPDE() <= 1) 
      {
        BasePDE * ptPDE =  ptdomain_->GetPDE(pdenumber);
        ptPDE->SolveStepStatic(level);   

        while(ptPDE->TestError(level) && numrepeat<maxnumrepeat)
          {
            ptPDE->RefineMesh();

            level++; // new level of grid

            if (meshesInfo)
              {
                ptMeshes_->Init(ptdomain_->GetGrid());
                ptPDE->WriteErrorInfo(ptMeshes_); // info for old mesh
                PrintSeqMeshes(); // print new mesh
              }

            ptdomain_->Update(level);

            ptPDE->SolveStepStatic(level);

            numrepeat++;
          }

        ptPDE->PostProcess(level);
        ptdomain_->PrintGrid(level);
        ptPDE->WriteResultsInFile();
      }
    else
      {
        Error("Adaptivity for Coupled Field Problem is not implemented",__FILE__,__LINE__);
     
        //     ptdomain_->GetCoupledPDE()->SolveStepStatic(level);
        //  ptdomain_->PrintGrid(level);
        // ptdomain_->GetCoupledPDE()->WriteResultsInFile();
      }
 
  }

}
