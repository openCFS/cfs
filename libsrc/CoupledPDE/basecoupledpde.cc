#include <fstream>
#include <iostream>
#include <string>

#include "basecoupledpde.hh"
#include <CoupledPDE/pdecoupling.hh>
#include <DataInOut/ParamHandling/ConfFile.hh>

 
namespace CoupledField
{

  BaseCoupledPDE::BaseCoupledPDE(StdVector<BasePDE*> & PDEs,
				 StdVector<PDECoupling*> & Couplings,
				 Grid *aptgrid, 
				 BCs *aptBCs, 
				 FileType *aInFile, 
				 WriteResults * aOutFile)
{
  ENTER_FCN( "BaseCoupledPDE::BaseCoupledPDE" );

  InFile_     = aInFile;
  OutFile_    = aOutFile;
  ptgrid_     = aptgrid;
  ptBCs_      = aptBCs;
  PDEs_       = PDEs;
  Couplings_  = Couplings;

  actlevel_ = 0;
  NumPDEs_ = PDEs.GetSize();

  
  // get analysis type
  std::string analysis;
#ifndef XMLPARAMS
  conf->get("analysis", analysis);
#else
  params->Get( "type", analysis, "analysis" );
#endif

  if (analysis=="static") 
    analysistype_ = STATIC;
  else if (analysis=="transient")
    analysistype_ = TRANSIENT;
  else
    Error("Analysis Type not supported",__FILE__,__LINE__);

  coupledpdename_ = "CoupledPDE: ";

  for (Integer actPDE=0; actPDE < PDEs.GetSize()-1; actPDE++)
    coupledpdename_ += PDEs[actPDE] -> GetName() + "+";

  coupledpdename_ += PDEs[PDEs.GetSize()-1] -> GetName();
}

BaseCoupledPDE::~BaseCoupledPDE()
{
  ENTER_FCN( "BaseCoupledPDE::~BaseCoupledPDE()" );

}




  // time stepping params, provided by "TransientDriver"
  void BaseCoupledPDE::SetTimeSteppingParams(Integer numstep, Double firstdt, Integer isavebegin, 
					     Integer isaveend, Integer isaveincr)
  {
    ENTER_FCN( "BaseCoupledPDE::SetTimeSteppingParams" );

    numstep_    =  numstep; 
    firstdt_    =  firstdt;    
    isavebegin_ =  isavebegin; 
    isaveend_   =  isaveend;   
    isaveincr_  =  isaveincr;  
  }




// perform on every pde a pre step (before solving transient step
void BaseCoupledPDE::PreStepTrans(const Integer kstep, const Double asteptime,
				  const Integer level, const Boolean reset) 
{
  ENTER_FCN( "BaseCoupledPDE::PreStepTrans" );

    for (Integer i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->PreStepTrans(kstep, asteptime, level, reset);
};




// perform on every pde a post step (after solving transient step
void BaseCoupledPDE::PostStepTrans(const Integer kstep, const Double asteptime, 
				   const Integer level) 
{
  ENTER_FCN( "BaseCoupledPDE::PostStepTrans" );

    for (Integer i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->PostStepTrans(kstep,asteptime,level);
}



void BaseCoupledPDE::InitStepTransCoupled(Double stepTime) 
{
  ENTER_FCN( "BaseCoupledPDE::InitStepTransCoupled" );

    for (Integer i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->InitStepTransCoupled(stepTime);
}



void BaseCoupledPDE::InitTimeStepping(const Double dt)
{
  ENTER_FCN( "BaseCoupledPDE::InitTimeStepping" );

    for (Integer i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->InitTimeStepping(dt);
};


void BaseCoupledPDE::WriteGeneralPDEdefines()
{
  ENTER_FCN( "BaseCoupledPDE::WriteGeneralPDEdefines" );

    for (Integer i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->WriteGeneralPDEdefines();
}


// ======================================================
// POSTPROC SECTION
// ======================================================

// Do Postprocessing as descriped in conf file
void BaseCoupledPDE::PostProcess(const Integer level) 
{
  ENTER_FCN( "BaseCoupledPDE::PostProcess" );

    for (Integer i=0; i<PDEs_.GetSize(); i++)
      PDEs_[i]->PostProcess(level);
}




} // end of namespace
