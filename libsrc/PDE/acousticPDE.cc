#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

namespace CoupledField
{

AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
			 WriteResults *aptOut)
:BasePDE(aptgrid,aptbcs,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::AcousticPDE " << std::endl;
#endif

  dofspernode_=1;

  pdename_    ="acoustic";
  pdematerialclass_ = "fluid";

  laststepcalc_=0;

  conf->getsubdompde(subdoms_,pdename_);

  AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);  
  NumPDENodes_ = PDE2MeshNode_.size();

  size_ = NumPDENodes_;

  sol_.reshape(dofspernode_, NumPDENodes_);
  sol_.init();

  sol_der1_.reshape(dofspernode_, NumPDENodes_);
  sol_der1_.init();
  
  sol_der1_old_.reshape(dofspernode_, NumPDENodes_);
  sol_der1_old_.init();
  
  sol_der1_.reshape(dofspernode_, NumPDENodes_);
  sol_der1_.init();
  
  sol_der2_.reshape(dofspernode_, NumPDENodes_);
  sol_der2_.init();

  sol_old_.reshape(dofspernode_, NumPDENodes_);
  sol_old_.init();
  
  sol_der2_old_.reshape(dofspernode_, NumPDENodes_);
  sol_der2_old_.init();
   
  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  conf->ifget("alpha_NM",alpha_,pdename_); 
  conf->ifget("beta_NM",beta_,pdename_); 
  conf->ifget("gamma_NM",gamma_,pdename_);

  with_absBCs_=FALSE;
  std::string absBCs="no";
  conf->ifget("absorbingBCs",absBCs,pdename_);
  if (absBCs == "yes")
    {
      with_absBCs_ = TRUE;
      //Error("Currently no absorbing BCs supported",__FILE__,__LINE__);
    }

  SetMatrixFactors();
  ReadBCs(pdename_);

}


void AcousticPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_ = RSCALAR;
  GraphType_  = NODEGRAPH; 

  SystemMatrix_     = TRUE;
  StiffnessMatrix_  = TRUE;
  MassMatrix_       = TRUE;

  if (with_absBCs_)
    DampingMatrix_  = TRUE;
}


void AcousticPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;       // factor for stiffness matrix

  matrix_factor_[1] = 0.0; 
  if (with_absBCs_)
    matrix_factor_[1] = 1.0*a1_; // factor for damping matrix

  matrix_factor_[2] = 0.0;       // factor for convection matrix
  matrix_factor_[3] = 1.0*a0_;   // factor for mass matrix
}


void AcousticPDE::ComputeRHS(const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acoustic3dPDE::ComputeRHS" << std::endl;
#endif

  Integer n;
  Integer node;

  Vector<Double> coeffMass, coeffDamp;
  Vector<Double> elemvec;
  Array<Double> temp;
  std::list<Integer> nodes_hd;
  Integer i;

  coeffMass = (sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_);
   
  algsys_->UpdateRHS(MASS,coeffMass.get());

  // damping matrix part
  if (with_absBCs_) 
    {
      coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;   
 
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());

      coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;  
       
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());
   }
   
}

void AcousticPDE::SolveStepTrans(const Integer kstep, const Double asteptime, 
				   const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=FALSE;

  if (laststepcalc_==kstep && kstep!=0) Recalc=TRUE;
  else laststepcalc_= kstep;

  Double * ptsol;

  Integer update,job;

  if (kstep==0)
    {
      update = 0;
      job = 1;
      SetupMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      algsys_->InitRHS();
      algsys_->InitMatrix(SYSTEM);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      ComputeRHS(lasttimecalc_);
    }
  else
    {
      update = 1;
      job    = 3;
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
    };

  SetBCs(level,update,lasttimecalc_);
  algsys_->CalcPrecond(job);
  algsys_->Solve();
  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k = 0;
  
  for (Integer i=0; i<NumPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptsol[k++];
  
  if (InfoPrint)
    (*infofile) << "maxnode:" <<  ptgrid_->GetMaxnumnodes(level) << std::endl;

  // calculation of derivatives of solution
  CalcDerSol(); 
}

void AcousticPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::WriteResultsInFile" << std::endl;
#endif

  Array<Double> sol_mesh, solder1_mesh, solder2_mesh;
  Array<Double> solArray_, sol_der1Array_, sol_der2Array_;
  
  solArray_ = sol_;
  sol_der1Array_ = sol_der1_;
  sol_der2Array_ = sol_der2_;

  
  TransformNodeSolution(sol_mesh,sol_,PDE2MeshNode_);
  TransformNodeSolution(solder1_mesh,sol_der1_,PDE2MeshNode_);
  TransformNodeSolution(solder2_mesh,sol_der2_,PDE2MeshNode_);

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"vp");
//       OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"vp_1der");
//       OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteNodeSolution(sol_mesh,laststepcalc_,lasttimecalc_,"fluid potential");
//       OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
//       OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }

}

void AcousticPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::CalcParameters" << std::endl;
#endif

  a2_=1.0/(beta_*dt);
  a0_=a2_*(1/dt);
  a1_=gamma_*a2_;
  a3_=0.5/(beta_)-1.0;
  a4_=gamma_/beta_-1.0;
  a5_=0.5*dt*(a4_-1.0);
  a6_=dt*(1-gamma_);
  a7_=gamma_*dt;
}

void AcousticPDE :: CalcDerSol()
{
#ifdef TRACE
  (*trace) << " entering  AcousticPDE :: CalcDerSol() " << std::endl;
#endif

   sol_der2_=(sol_ - sol_old_)*a0_ - (sol_der1_old_)*a2_ - sol_der2_old_*a3_;
   sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;
}

void AcousticPDE::SaveSolAsPrevStep()
{
#ifdef TRACE
  (*trace) << " entering  AcousticPDE::SaveSolAsPrevStep() " << std::endl;
#endif

  sol_old_=sol_;
  sol_der1_old_=sol_der1_;
  sol_der2_old_=sol_der2_;  
}


void AcousticPDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
    (*trace) << "entering AcousticPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;
  Matrix<Double> ptCoord;


  BaseFE * ptEl;
  Double coeffstiff, coeffmass;
  Vector<Integer> connecth, connect_PDE;
  std::vector<Elem*> elemssd;

  Integer i, j;
 
  for (i=0; i<subdoms_.size(); i++)
    {
      //compute material coefficient
      const Double density = materialData_[i].GetDensity();
      const Double compressibility = materialData_[i].GetCompressibility();

      coeffmass  = density*density/compressibility;
      coeffstiff = density;


      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl = elemssd[j]->ptElem;
    
	  BaseForm * bilinear_mass  = new MassInt(ptEl, coeffmass);
	  BaseForm * bilinear_stiff = new LaplaceInt(ptEl, coeffstiff);

	  connecth=elemssd[j]->connect;
	  GetElemCoords(connecth, ptCoord, level); 

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif     
	  
	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), STIFFNESS);

	  // mass part
	  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), MASS);

	  delete bilinear_stiff;
	  delete bilinear_mass;

	}
    }

  // BEGIN DAMPING MATRIX PART
  if (with_absBCs_) {
    std::vector<Elem*>  DomainBnd;
    Double coeffdamp;
    std::cout<<"Aqui en damping"<<std::endl;
 
    conf->getliststr("bnd_for_absBCs",bnd_absBCs_,"acoustic");
    DomainBnd=ptBCs_->getFacesBC(bnd_absBCs_[0],level);

    for (j=0; j< DomainBnd.size(); j++)
      {
	const Double density = materialData_[0].GetDensity();
	const Double compressibility = materialData_[0].GetCompressibility();
	coeffdamp = density/((sqrt(compressibility/density)));
      
	ptEl=DomainBnd[j]->ptElem;
	// Here MassInt is used to calculate the damping matrix from the surface elements
	BaseForm * linear_damp = new MassInt(ptEl,coeffdamp);

	Integer ii;
	Integer elsize=(DomainBnd[j]->connect).size();
	connecth.Resize(elsize);
	for (ii=0; ii<elsize; ii++)
	  connecth[ii]=(DomainBnd[j]->connect)[ii];

	GetElemCoords(connecth, ptCoord, level);
	// CHANGE connecth
	Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	linear_damp->CalcElementMatrix(ptCoord, elemmat);

#ifdef DEBUG
	(*debug) << "Connection array  " << std::endl;
	(*debug)  << connecth  << std::endl;

	(*debug) << "Dampingmatrix, ElementNumber  " <<   j << std::endl;
	(*debug) << elemmat << std::endl;
#endif


	algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), DAMPING);

	delete linear_damp;

      }
  }
  // END DAMPING MATRIX PART



#ifdef TRACE
    (*trace) << "Leaving AcousticPDE::SetupMatrices" << std::endl;
#endif
  }


}


