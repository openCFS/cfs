#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
#include "newmark.hh"
#include "newmarkdamp.hh"
#include <DataInOut/WriteInfo.hh>


namespace CoupledField
{

AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
			 WriteResults *aptOut)
:BasePDE(aptgrid,aptbcs,aptFileType,aptOut,aptTimeFunc)
{
  ENTER_FCN( "AcousticPDE::AcousticPDE" );

  dofspernode_=1;

  pdename_    ="acoustic";
  pdematerialclass_ = "fluid";

  isaxi_ = FALSE;
  std::string subtype;
  conf->ifget("subtype",subtype,pdename_);
  if (subtype == "axi")
    isaxi_ = TRUE;

  laststepcalc_=0;

  conf->getsubdompde(subdoms_,pdename_);

  AssignPDENodeNumbers(mesh2PDENode_, pde2MeshNode_, subdoms_);  
  AssignPDEElemNumbers(mesh2PDEElem_, pde2MeshElem_, subdoms_);
  numPDENodes_ = pde2MeshNode_.size();
  numElems_ = pde2MeshElem_.size();

  size_ = numPDENodes_;

  // Initalize solution class
  sol_->SetNumSolutions(1);
  sol_->SetSolutionType(ACOU_POTENTIAL);
  sol_->SetNumNodes(numPDENodes_);
  sol_->SetNumDofs(dofspernode_);
  sol_->Init(0.0);

  
  
  with_fracdamping_=FALSE;
  std::string dampstr;
  conf->ifget("damping",dampstr,pdename_);

  if (dampstr == "fractional")
    {
      Info->PrintF(pdename_,"\n Attenuation according to power law, number of memory is %g\n\n", frac_memory_);
       with_fracdamping_ = TRUE;
       conf->get("frac_memory",frac_memory_,pdename_);
       damping_type_ = FRACTIONAL;
    }

  with_absBCs_=FALSE;
  conf->ifgetliststr("bnd_for_absBCs",bnd_absBCs_,pdename_); 
  if (bnd_absBCs_.size() > 0)
    {
      with_absBCs_ = TRUE;
      if ( damping_type_ == NONE) damping_type_ = ABCDamp;
    }

  ReadBCs(pdename_);


  if (analysistype_==HARMONIC)
    {
      conf->get("frequency", freq_, pdename_);
      solIm_.SetNumSolutions(1);
      solIm_.SetSolutionType(ACOU_POTENTIAL);
      solIm_.SetNumNodes(numPDENodes_);
      solIm_.SetNumDofs(dofspernode_);
      solIm_.Init(0);

    }

  // set analysis parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, bnd_absBCs_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&mesh2PDENode_);

#ifdef USE_OLAS
  assemble_->SetMatrixEntryType(DOUBLE);
  assemble_->SetMatrixStorageType(SPARSE_NONSYM);
#else
  assemble_->SetMatrixType(RSCALAR);
#endif

  assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
  assemble_->SetPtrBCs(ptBCs_);
  assemble_->SetPtr2Sol(sol_);
  assemble_->SetPtr2TimeFnc(ptTimeFunc_);

  if (with_absBCs_ || with_fracdamping_)
    assemble_->NeedDampingMatrix();

  ReadMaterialData();
   
  DefineIntegrators(actlevel_);

  ReadSavings();
  if (savederiv1_)
    {
      sol_der1Array_.SetNumSolutions(1);
      sol_der1Array_.SetNumNodes(numPDENodes_);
      sol_der1Array_.SetSolutionType(ACOU_POTENTIAL_DERIV1);
      sol_der1Array_.SetNumDofs(1);
      sol_der1Array_.Init(0);
    }

  if (savederiv2_)
    {
      sol_der2Array_.SetNumSolutions(1);
      sol_der2Array_.SetNumNodes(numPDENodes_);
      sol_der2Array_.SetSolutionType(ACOU_POTENTIAL_DERIV2);
      sol_der2Array_.SetNumDofs(1);
      sol_der2Array_.Init(0);
    }
  
}


void AcousticPDE::DefineIntegrators(const Integer level)
{
  ENTER_FCN( "AcousticPDE::DefineIntegerators" );

  Boolean nonLin = FALSE;

  for (Integer actSD = 0; actSD < subdoms_.size(); actSD++)
    {
      Double density = materialData_[actSD].GetDensity();
      Double compressibility = materialData_[actSD].GetCompressibility();

      //stiffness integrator
      BaseForm * bilinearStiff = new LaplaceInt(density,isaxi_);	  
      assemble_->AddIntegrator(bilinearStiff, subdoms_[actSD], STIFFNESS, nonLin);

      //mass integrator
      Double coeffmass  = density*density/compressibility;
      BaseForm * bilinear_mass  = new MassInt(coeffmass, dofspernode_, isaxi_);
      assemble_->AddIntegrator(bilinear_mass, subdoms_[actSD], MASS, nonLin);
    }

  //surface-integration
  // BEGIN DAMPING MATRIX PART: Absorbing boundaries
  if (with_absBCs_ && analysistype_!=HARMONIC) 
    {  
      for (Integer actSD = 0; actSD < bnd_absBCs_.size(); actSD++)
	{
	  //currently hard-coded!!
	  Double density = materialData_[0].GetDensity();
	  Double compressibility = materialData_[0].GetCompressibility();
	  Double coeffdamp = density/((sqrt(compressibility/density)));
	  
	  BaseForm * bilinear_damp = new MassInt(coeffdamp,dofspernode_, isaxi_);
	  assemble_->AddSurfIntegrator(bilinear_damp,  bnd_absBCs_[actSD], DAMPING, nonLin);
	}
    }
  }


// ======================================================
// SOLVING SECTION
// ======================================================

void AcousticPDE :: InitTimeStepping(const Double dt)
{
  ENTER_FCN( "AcousticPDE::InitTimeStepping" );

  if (with_fracdamping_)
    {
      //currently the parameter y is taken from the first subdomain
      //=> currently just one subdomain makes sense
      //      Double y = materialData_[0].GetDampingBeta();
      //      TS_alg_ = new NewmarkDamp(pdename_, algsys_, dofspernode_, numPDENodes_, damping_type_,
      //			      frac_memory_,y);
    }
  else
    TS_alg_ = new Newmark(pdename_, algsys_, dofspernode_, numPDENodes_, damping_type_);

  TS_alg_->Init(matrix_factor_, dt);

}




// ======================================================
// COUPLING SECTION
// ======================================================


void AcousticPDE::InitCoupling(PDECoupling * Coupling)
{
  ENTER_FCN( "AcousticPDE::InitCoupling" );
  
  pdeIsCoupled_ = TRUE;
  ptCoupling_   = Coupling;
  
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      if (ptCoupling_->GetOutputQuantity(i) == "acousticforce")
	{
	  // Intialize the memory of the coupling values
	  ptCoupling_->CreateStoreSol(i,ACOU_FORCE,isComplex_);
	}
    }

  iterCoupledCounter_ = 0;
}




void AcousticPDE::CalcOutputCoupling()
{
  ENTER_FCN( "AcousticPDE::CalcOutputCoupling" );

  Integer dof;
  std::string quantity;
  std::vector<Elem*> * couplingElems = NULL;
  std::vector<Elem*> * interfaceVolElems = NULL;
  std::vector<Integer> * couplingNodes = NULL;
  std::vector<MaterialData*> * couplingMaterials = NULL;
  BaseStoreSol * values = NULL;
  
  TRY_CAST

  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);
      ptCoupling_->GetOutputValues(i, values);
      PTRCAST(values,StoreSol<Double>,temp);

      switch(ptCoupling_->GetOutputType(i))
	{
	case NODE:	  
	  if (quantity == "acousticforce")
	    {
	      ptCoupling_->GetOutputElements(i, couplingElems);
	      ptCoupling_->GetOutputNodes(i, couplingNodes);
	      ptCoupling_->GetOwnMaterials(i, couplingMaterials);
	      ptCoupling_->GetInputNeighbourElems(i, interfaceVolElems);
	      dof = ptCoupling_->GetOutputDof(i);
	    
	      CalcMechCouplingRHS(couplingElems, *couplingNodes, couplingMaterials, *temp, dof, interfaceVolElems);

	    }	  
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
    }

  CATCH_CAST
}

void AcousticPDE::CalcMechCouplingRHS(std::vector<Elem*> * couplingElems, 
				      std::vector<Integer> & couplingNodes,
				      std::vector<MaterialData*> * couplingMaterials,
				      StoreSol<Double>& elemCouplingSols,
				      Integer couplingdof,
				      std::vector<Elem*> * interfaceVolElems)
{
  ENTER_FCN( "AcousticPDE::CalcMechCouplingRHS" );

  Double density=0;
   
  elemCouplingSols.Init(0.0);

  for (Integer actElem=0; actElem<couplingElems->size(); actElem++)
    {
      Elem * actCoupleElem = (*couplingElems)[actElem];
      
      BaseFE * ptElem          = actCoupleElem->ptElem;
      Vector<Integer> connecth = actCoupleElem->connect;
      
      Matrix<Double> ptCoord; 
      GetElemCoords(connecth, ptCoord, actlevel_);

      Boolean found = FALSE;
      
      // get correct density belonging to the neighbouring element of the interface
      for (Integer actSD = 0; actSD < subdoms_.size(); actSD++)
	{  
	  if ((*interfaceVolElems)[actElem]->namesd ==  subdoms_[actSD])
	    {
	      density = materialData_[actSD].GetDensity();
	      found = TRUE;	      
	    }
	}
      
      if (found ==FALSE) 
	{
	  mycout << "Could not found correct density to compute acoustic pressure forces!!!!!"
		 << myendl;
	  mycout << "Take density of acoustic subdomain 1" << myendl;
	  density = materialData_[0].GetDensity();
	}
  
          
      BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
      Matrix<Double> elemmat;
      bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
      delete bilinear_mass;	  


      Vector<Integer> connect_PDE;
      Mesh2PDENode(connect_PDE, connecth, mesh2PDENode_);
      
      Vector<Double> sol;
      GetDerivSolVecOfElement(sol, connect_PDE);
      
      Vector<Double> forceOnElem = elemmat * sol;
      // force has to be added on RHS with negative sign
      forceOnElem *= -1;


      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      Vector<Double> n;
      CalcLineNormalVec(n, *actCoupleElem, *(*interfaceVolElems)[actElem]); // points outward own domain

      
      for (Integer actNode=0; actNode<ptCoord.GetSizeRow(); actNode++)
	{
	  Integer nodePos = 0;      
	  
	  while(connecth[actNode] != couplingNodes[nodePos] && nodePos < couplingNodes.size()) 
	    nodePos++;
	  
	  for (Integer actDof=0; actDof < couplingdof ; actDof++)  
	    elemCouplingSols(nodePos,actDof) += forceOnElem[actNode] * n[actDof];
	}      
    }
}




Boolean AcousticPDE::HasOutput(std::string output)
{
  ENTER_FCN( "AcousticPDE::HasOutput" );

  if (output == "acousticforce")
    return TRUE;

  return FALSE;
}



// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void AcousticPDE::WriteResultsInFile()
{
  ENTER_FCN( "AcousticPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
  int commrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  if (!commrank) 
  	{
#endif
  StoreSol<Double> sol_mesh, solder1_mesh, solder2_mesh, solIm_mesh;
 
  if (analysistype_==HARMONIC)
    {
      sol_->TransformNodeSolution(sol_mesh,pde2MeshNode_,ptgrid_,actlevel_);
      sol_->TransformNodeSolution(solIm_mesh,pde2MeshNode_,ptgrid_,actlevel_);      
      outFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw realpart,");
      outFile_->WriteNodeSolution(solIm_mesh,laststepcalc_,lasttimecalc_,"fluid potential, cw imagpart, ");
    }
  else
    {  

      if (savesol_)
	  sol_->TransformNodeSolution(sol_mesh,pde2MeshNode_,ptgrid_,actlevel_);

      if (savederiv1_) 
	{
	  sol_der1Array_.SetSolVector(ACOU_VELOCITY,getS1());
	  sol_der1Array_.TransformNodeSolution(solder1_mesh,pde2MeshNode_,ptgrid_,actlevel_);
	}

      if (savederiv2_)
	{
	  sol_der2Array_.SetSolVector(ACOU_VELOCITY,getS2());
	  sol_der2Array_.TransformNodeSolution(solder2_mesh,pde2MeshNode_,ptgrid_,actlevel_);
	}
      
      if (outFile_->IsGMV())
	{
	  if (savesol_)
	    outFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"vp");
	  if (savederiv1_)
	    outFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"vp_1der");
	  if (savederiv2_)
	    outFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"vp_2der");
	}
      else
	{
	  if (savesol_)
	    outFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential");
	  

	  if (savederiv1_)
	    outFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv.");
	  if (savederiv2_)
	    outFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv.");
	}
    }
    #ifdef PARALLEL
    }//!commrank
    #endif
}




}


