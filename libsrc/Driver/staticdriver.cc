#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"

#include "DataInOut/GMV/outGMV.hh"

namespace CoupledField
{

StaticDriver :: StaticDriver(Domain * adomain)
:BaseDriver(adomain)
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::StaticDriver" << std::endl;
#endif

}

StaticDriver :: ~StaticDriver()
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::~StaticDriver" << std::endl;
#endif

}

void StaticDriver :: SetupMatricesPDE(const Integer pdenumber, const Integer type)
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SetUpMatricesPDE" << std::endl;
#endif

  ptdomain_->GetPDE(pdenumber)->SetupMatrices(type);
}

void StaticDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SolveProblem" << std::endl;
#endif
  Integer level=0;
  Integer pdenumber  = 0;

  ptdomain_->GetPDE(pdenumber)->SolveStepStatic(level);

  ptdomain_->PrintGrid(level);
  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
}

void StaticDriver :: SolveProblemAdaptSpace()
{
#ifdef TRACE
  (*trace) << "entering SolveProblemAdapt::SolveProblemAdaptSpace " << std::endl;
#endif

  Integer level = 0;
  Integer pdenumber = 0;
  Integer numrepeat=0;

// create files for printing seq. of meshes  
  Boolean PrintMeshes=TRUE;
  // conf->is_there("sequence_of_meshes");
  ptMeshes_=new WriteResultsGMV("meshes");
  ptMeshes_->Init(ptdomain_->GetGrid());
  if (PrintMeshes) PrintSeqMeshes();
  //

  BasePDE * ptPDE;
  ptPDE=ptdomain_->GetPDE(pdenumber);

  Integer maxnumrepeat;
  conf->get("maxnumrepeat",maxnumrepeat,"SpaceAdaptivity");

  ptPDE->SetMatrixFactors();
  ptPDE->SolveStepStatic(level);

  if (InfoPrint)
    (*infofile) << " ---------- step 0 ----------------- " << std::endl;

  while (ptPDE->TestError() && numrepeat!=maxnumrepeat ) {
    ptPDE->RefineMesh();
        
    if (PrintMeshes) 
      {
	ptPDE->PrintMeshesInfo(ptMeshes_);
	PrintSeqMeshes();
      }

    ptdomain_->Update(level);

     ptPDE->SetMatrixFactors(); 
     ptPDE->SolveStepStatic(level);

    numrepeat++;

    if (InfoPrint)
      (*infofile) << " ---------- step " << numrepeat << " ----------------- " << std::endl;
  }
    
  ptdomain_->PrintGrid(level);
  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

  if (PrintMeshes) { 
	ptPDE->PrintMeshesInfo(ptMeshes_);
	delete ptMeshes_;
  }
}
} // end of namespace
