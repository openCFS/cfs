#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include "blocknodeEQN.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "smoothPDE.hh" 

namespace CoupledField
{

SmoothPDE::SmoothPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
  ENTER_FCN( "SmoothPDE::SmoothPDE" );
  
  pdename_ = "smooth";
  pdematerialclass_ = "piezo";
  firstTurn_ = TRUE;

  // No time step algorithm for this PDE
  isAlwaysStatic_ = TRUE;

  
  
  
#ifndef XMLPARAMS
    conf->getstr("subtype", subType_, pdename_ );

    if (subType_ == "3d")
      {
	dofspernode_ = 3;
	Info->PrintF("", "=== 3D PROBLEM\n");
      }
    else if (subType_ == "axi")
      {
	isaxi_ = TRUE;
	dofspernode_ = 2;
	Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
      }
    else if (subType_ == "planeStrain" )
      {
	dofspernode_ = 2;
	Info->PrintF("", "=== PLAIN STRAIN PROBLEM\n");
      }
    else
      {
	std::string errmsg = "Subtype " + subType_ + " is not defined for";
	errmsg += " PDEs of type " + pdename_ + '\n';
	Info->Error( errmsg, __FILE__, __LINE__ );
      }
#else

    // Get problem geometry and PDE subtype
    params->Get( "subType", subType_, pdename_ );
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode_ = 3;
          }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 2;
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      dofspernode_ = 2;
    }
    else
      {
	std::string errmsg = "Subtype '" + subType_;
	errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
	errmsg += "geometry '" + probGeo + "'\n";
	Info->Error( errmsg, __FILE__, __LINE__ );
      }
#endif

    // =====================================================================
    // set solution information
    // =====================================================================
    
  solTypes_ = SMOOTH_DISPLACEMENT;
  solDofs_ = dofspernode_;
  
  method_ = "mechanic";
#ifndef XMLPARAMS
  conf->ifget("method", method_, pdename_ );
#else

#endif

  //is a nonlinear PDE, since in each iteration, we have to setup the matrices new!
  nonLin_ = TRUE;
}


void SmoothPDE::DefineIntegrators(const Integer level)
{
  ENTER_FCN( "SmoothPDE::DefineIntegerators" );

  Boolean nonLin = FALSE;

  // Weigthening factors for smoothening 
  factor_.Resize(numElems_);
   for (Integer i=0; i<factor_.GetSize(); i++)
     factor_[i] = 1.0;

  for (int actSD = 0; actSD < subdoms_.GetSize(); actSD++)
    {
      // ==============  add stiffness ===========================================

      MaterialData actSDMat(materialData_[actSD]);

      // ==============  add "standard" stiffness ===============================
      BaseForm * bilinearStiff;
#ifndef XMLPARAMS
      if (subType_ == "plainStrain")
#else
      if (subType_ == "planeStrain")
#endif
	bilinearStiff = new smoothPlainStrainInt(actSDMat);
      else if (subType_ == "3d")
	bilinearStiff = new smooth3DInt(actSDMat);
      else if (subType_ == "axi")
	bilinearStiff = new SmoothAxiInt(actSDMat);
      else 
	Error("Unknown subtype in smooth PDE! ",__FILE__,__LINE__);

      assemble_->AddIntegrator(bilinearStiff, subdoms_[actSD], STIFFNESS, nonLin);
    }
}


void SmoothPDE::InitCoupling(PDECoupling * coupling)
{
  ENTER_FCN( "SmoothPDE::Initcoupling" );

  pdeIsCoupled_ = TRUE;
  ptCoupling_   = coupling; 

  // input couplings
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
    {
      // check for input mechanic displacement
      if (ptCoupling_->GetInputQuantity(i) == MECH_DISPLACEMENT)
	{
	  numDirichletBCs_ += (dofspernode_ * ptCoupling_->GetInputNumNodes(i));
	}

    }

  // output couplings
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      // check for output displacement
      if (ptCoupling_->GetOutputQuantity(i) == SMOOTH_DISPLACEMENT)
	{
	ptCoupling_->CreateCouplingVector(i,isComplex_); 
	}
    }

  // now overwrite number of Dirichlet BCs due to coupling 
  assemble_->SetNumDirichlet(numDirichletBCs_);
}

void SmoothPDE::PreStepStatic(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset)
{
  ENTER_FCN( "SmoothPDE::PreStepStatic" );

  algsys_->InitRHS();
  algsys_->InitSol();
  assemble_->InitMatrices();

  assemble_->SetReassemble();

}


void SmoothPDE::StepStaticNonLin(const Integer kstep, const Double aTime,
				const Integer level, const Boolean reset)
{
  ENTER_FCN( "SmoothPDE::StepStaticNonLin" );

  Integer job = 1;
  Double * ptsol;

  assemble_->AssembleMatrices(level);
  assemble_->AssembleSrcRHS(level);
  
  SetBCs(level,0);

#ifdef USE_OLAS
  algsys_->BuildInDirichlet();
  algsys_->SetupPrecond(job);
#else
  algsys_->CalcPrecond(job);
#endif

  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

   // save solution
  sol_->CopyFromAlgSysDataPointer(ptsol);
  //  sol_->SetAlgSysDataPointer(ptsol);
}


void SmoothPDE:: PostStepStatic(const Integer kstep, const Double asteptime,
				const Integer level)
{
  ENTER_FCN( "SmoothPDE::PostStepStatic" );
  
  iterCoupledCounter_++;
}




void SmoothPDE::CalcOutputCoupling()
{
  ENTER_FCN( "SmoothPDE::CalcOutputCoupling" );

  SolutionType quantity;
  StdVector<Integer> * couplingnodes;
  CFSVector * values;

  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);

      switch(ptCoupling_->GetOutputType(i))
	{

	case NODE:
	  
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  ptCoupling_->GetOutputValues(i, values);

	  if (quantity == SMOOTH_DISPLACEMENT)
	    {
	      sol_->NodeSolutionToCoupling(*values,*couplingnodes);
	    }
	  
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}

    }
}

void SmoothPDE::WriteResultsInFile(Integer stepOffset,
				   Double timeOffset)
{
  ENTER_FCN( "SmoothPDE::WriteResultsInFile" );
  
  NodeStoreSol<Double> * solTransient;
  
  Integer actStep = laststepcalc_ + stepOffset;
  Double actTime = lasttimecalc_ + timeOffset;

  if (analysistype_ == STATIC ||
      analysistype_ == HARMONIC) {
    solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
  }
  else
    Error("SmoothPDE: Only static and transient results can be written out",
	  __FILE__, __LINE__);
}



Integer SmoothPDE::GetNrBCDof(const std::string & dofStartString)
{
  ENTER_FCN( "SmoothPDE::GetNrBCDof" );

  Integer nrActDof;
  
  if (dofStartString == "ux")
    nrActDof = 1;
  else 
    if (dofStartString == "uy")
      nrActDof = 2;
    else
      if (dofspernode_ == 3)
	if (dofStartString == "uz")
	  nrActDof = 3;
	else
	  {
	    Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
	    std::cerr << dofStartString << std::endl;
	  }
      else
	{
	  std::cerr << dofStartString << std::endl;
	  Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
	}
  
  return nrActDof;
}


Integer SmoothPDE::GetBCDof(const std::string dofString)
{
  ENTER_FCN( "SmoothPDE::GetBCDof" );

  if (dofString == "ux")
    return 1;
  else
    if (dofString == "uy")
      return 2;
    else
      if (dofString == "uz")
	return 3;
      else
	Error("The direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
}
  

Boolean SmoothPDE::HasOutput(SolutionType output)
{
  ENTER_FCN( "SmoothPDE::HasOutput" );
  
  if (output == SMOOTH_DISPLACEMENT)
    return TRUE;

  return FALSE;
}

void SmoothPDE::ReadStoreResults()
{
  ENTER_FCN( "SmoothPDE::ReadStoreRestuls");
}


} // end namespace CoupledField
