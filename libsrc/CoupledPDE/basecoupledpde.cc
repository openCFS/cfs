#include <fstream>
#include <iostream>
#include <string>

#include "basecoupledpde.hh"
#include <CoupledPDE/pdecoupling.hh>
#include <DataInOut/conffile.hh>
 
namespace CoupledField
{

BaseCoupledPDE::BaseCoupledPDE(std::vector<BasePDE*> & PDEs,
			       std::vector<PDECoupling*> & Couplings,
			       Grid *aptgrid, 
			       BCs *aptBCs, 
			       FileType *aInFile, 
			       WriteResults * aOutFile)
{
#ifdef TRACE
  (*trace) << "entering BaseCoupledPDE::BaseCoupledPDE" << std::endl;
#endif

  InFile_     = aInFile;
  OutFile_    = aOutFile;
  ptgrid_     = aptgrid;
  ptBCs_      = aptBCs;
  PDEs_       = PDEs;
  Couplings_  = Couplings;

  actlevel_ = 0;
  NumPDEs_ = PDEs.size();
  coupledpdename_ = "coupling";

  // get analysis type
  std::string analysis;
  conf->get("analysis", analysis);


  if (analysis=="static") 
    analysistype_ = STATIC;
  else if (analysis=="transient")
    analysistype_ = TRANSIENT;
  else
    Error("Analysis Type not supported",__FILE__,__LINE__);
}

BaseCoupledPDE::~BaseCoupledPDE()
{
#ifdef TRACE
  (*trace) << " entering BaseCoupledPDE::~BaseCoupledPDE() " << std::endl;
#endif

}

} // end of namespace
