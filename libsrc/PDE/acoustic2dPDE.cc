#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acoustic2dPDE.hh" 
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

namespace CoupledField
{

Acoustic2dPDE::Acoustic2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, 
			     TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::Acoustic2dPDE " << std::endl;
#endif

  dofspernode_=1;
  pdename_    ="Acoustic2d";
 
  ptgrid_=aptgrid; 

  laststepcalc_=0;
  size_=ptgrid_->GetMaxnumnodes(0);

  sol_.Resize(size_);
  sol_.Init(0);

  sol_der1_.Resize(size_);
  sol_der1_.Init(0);

  sol_der1_old_.Resize(size_);
  sol_der1_old_.Init(0);

  sol_der1_.Resize(size_);
  sol_der1_.Init(0);

  sol_der2_.Resize(size_);
  sol_der2_.Init(0);

  sol_old_.Resize(size_);
  sol_old_.Init(0);

  sol_der2_old_.Resize(size_);
  sol_der2_old_.Init(0);

  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  conf->ifget("alpha_NM",alpha_,pdename_); 
  conf->ifget("beta_NM",beta_,pdename_); 
  conf->ifget("gamma_NM",gamma_,pdename_);

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

  with_absBCs_=FALSE;
  std::string absBCs="no";
  conf->ifget("absorbingBCs",absBCs,pdename_);
  if (absBCs == "yes")
    with_absBCs_ = TRUE;

  preComputeRHS();
}

void Acoustic2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;       // factor for stiffness matrix

  matrix_factor_[1] = 0.0; 
  if (with_absBCs_)
    matrix_factor_[1] = 1.0*a1_; // factor for damping matrix

  matrix_factor_[2] = 0.0;       // factor for convection matrix
  matrix_factor_[3] = 1.0*a0_;   // factor for mass matrix
}

void Acoustic2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, 
				    Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = RSCALAR; 
  graphtype  = NODEGRAPH; 

  matrixsystype[0] = SYSTEM;      // memory for the system matrix
  matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix

  matrixsystype[2] = 0;
  if (with_absBCs_)
    matrixsystype[2] = DAMPING;   // memory for the damping matrix

  matrixsystype[3] = 0;           // memory for the convection matrix
  matrixsystype[4] = MASS;        // memory for the mass matrix


  numdofpernode  = dofspernode_;
  numdirichlets  = 1;
  numconstraints = 0;
}

void Acoustic2dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;

  // This is a smaller matrix since it is just for linear 1D elements.
  Matrix<Double> elemmatbnd;
  Point<2> * ptCoord; 

  BaseElem * ptEl;

  Vector<Double> coeffm, coeffst, coeffdamp;
  CalcCoeff(coeffm, coeffst, coeffdamp);

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j< elemssd.size(); j++)
	{
	  ptEl=elemssd[j]->ptElem;

	  BaseForm<2> * bilinear_stiff = new LaplaceInt<2>(ptEl,dofspernode_);
	  BaseForm<2> * bilinear_mass  = new MassInt<2>(ptEl,dofspernode_);

	  Integer ii;
	  Integer elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];
  
	  ptCoord=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level); 

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;

	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif     

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, 
				      STIFFNESS);

	  // mass part
	  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffm[i];

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;

	  (*debug) << elemmat << std::endl;
#endif
      
	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_,
				      MASS);
  
	  delete bilinear_stiff;
	  delete bilinear_mass;
	  delete [] ptCoord;
	}
    }

  // BEGIN DAMPING MATRIX PART
  if (with_absBCs_) 
    {
      std::vector<Elem*>  DomainBnd;
      conf->getliststr("bnd_for_absBCs",bnd_absBCs_,pdename_);
      DomainBnd=ptBCs->getEdgesBC(bnd_absBCs_[0],level);

      for (j=0; j< DomainBnd.size(); j++)
	{
	  ptEl=DomainBnd[j]->ptElem;

	  BaseForm<2> * linear_damp = new DampInt<2>(ptEl,1);

	  Integer ii;
	  Integer elsize=(DomainBnd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(DomainBnd[j]->connect)[ii];

	  ptCoord=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  linear_damp->CalcElemMatrix(ptCoord, elemmatbnd);
	  elemmatbnd*=coeffdamp[0];

#ifdef DEBUG
	  (*debug) << "Connection array  " << std::endl;
	  (*debug)  << connecth  << std::endl;

	  (*debug) << "Dampingmatrix, ElementNumber  " <<   j << std::endl;
	  (*debug) << elemmatbnd << std::endl;
#endif

	  ptalgsys_->PutElemMatAlgSys(elemmatbnd.getinarray(), connecth.get(), connecth.size(), 
				      as_sysid_, as_sysid_, DAMPING);
	  delete linear_damp;
	  delete [] ptCoord;
	}
    }
  // END DAMPING MATRIX PART

#ifdef TRACE
  (*trace) << "Leaving Acoustic2dPDE::SetupMatrices" << std::endl;
#endif
}

void Acoustic2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val,val_tfunc;
  Double time=atime;

  val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,level);

  // set dirichlet boundary conditions
  Integer i,j=0;
  std::list<Integer> nodes;
  Integer sizebc=bcs_hd_.size();
  val = 0.0;

  for (i=0; i< bcs_hd_.size(); i++)
    {
#ifdef DEBUG
      (*debug) << " Dirichlet BC" << std::endl;
#endif
      nodes=ptBCs->GetNodesLevel(bcs_hd_[i]);
   
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  if (update==1)
	    {
#ifdef DEBUG
	      (*debug) << " node: " << node << " val: " << val << std::endl;
#endif

	      ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, SYSTEM);
	    }
	  else
	    { 
#ifdef DEBUG
	      (*debug) << " node: " << node << " val: " << val << std::endl;
#endif

	      ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, SYSTEM);
	    }
	}  
    }

  //inhomogeneous boundary conditions
  for (i=0; i<bcs_id_.size(); i++)
    {
      
#ifdef DEBUG
      (*debug) << " Inhomog. Dirichlet BC" << std::endl;
#endif

      nodes=ptBCs->GetNodesLevel(bcs_id_[i], level);

      val=val_tfunc*val_id_[i];

      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
          if (update==1)
            {

#ifdef DEBUG
	      (*debug) << " node: " << node << " val: " << val << std::endl;
#endif

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, SYSTEM);
            }
          else
            {

#ifdef DEBUG
	      (*debug) << " node: " << node <<  " val: " << val << std::endl;
#endif
              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, SYSTEM);
            }
	}
    }

#ifdef TRACE
  (*trace) << "leaving Acoustic2dPDE::SetBCs" << std::endl;
#endif
}

void Acoustic2dPDE::preComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::preComputeRHS" << std::endl;
#endif

  SetRHSFnc = FALSE;

  std::string isthererhs="no";
  conf->ifget("load_force",isthererhs,pdename_);

  if (isthererhs=="yes" ) {
    conf->getliststr("rhs_surfaces",rhs_surfaces_,pdename_);
    std::string rhs;
    conf->get("rhs_fnc",rhs,pdename_);
    ptRHSFnc_=FncReader(rhs);
    conf->get("rhs_arg_fnc",arg_rhs_,pdename_);
    SetRHSFnc=TRUE;
    directionFnc_.resize(2);
    conf->get("direction_rhs_x",directionFnc_[0],pdename_);
    conf->get("direction_rhs_y",directionFnc_[1],pdename_);
    }    
 
}

void Acoustic2dPDE::ComputeRHS(const Double atime, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::ComputeRHS" << std::endl;
#endif
  Integer n;
  Integer node;
  Integer matrix_id;
  Vector<Double> coeffMass;
  Vector<Double> coeffDamp;
  std::list<Integer> nodes_hd;
  Integer i;

  coeffMass = sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,MASS,coeffMass.get());

  // damping matrix part

  if (with_absBCs_) 
    {

      coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;
      ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,DAMPING,coeffDamp.get());

      coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;  
      ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,DAMPING,coeffDamp.get());
    }


  Integer level=0;          ////
  if (SetRHSFnc) {         // for standing waves

    Double density, compress;    // get density & compressity
    Integer i,matnum;
    for (i=0; i<subdoms_.size(); i++)
      {
	conf->get(subdoms_[i],matnum,"list_subdomains");
	MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");
      }

    Integer j,node1,node2;
    Double valfnc,integrShFnc,val,multiplier;
    Point<2> * ptCoordNodes;
    std::vector<Double> normal;
    normal.resize(2);
    Vector<Integer> connecth;
    std::vector<Elem*> edgesBC;  // vector of 1D-elements from mesh-file

    edgesBC=ptBCs->getEdgesBC(rhs_surfaces_[0],level);
    valfnc = ptRHSFnc_(atime*arg_rhs_*2*3.1416); // value of fnc at the timestep

    for (j=0; j< edgesBC.size(); j++) { // loop over surface elements
      ptCoordNodes=new Point<2>[2];

      connecth=edgesBC[j]->connect;
      ptgrid_->GetCoordNodesElem(connecth,ptCoordNodes,level);  
      calcNormal2Line(normal,ptCoordNodes[0],ptCoordNodes[1]); 

      multiplier=ScalarMult(normal,directionFnc_); // n * V /direction fnc/
      integrShFnc = edgesBC[j]->ptElem->getIntVal(ptCoordNodes);
      val=integrShFnc*multiplier*valfnc*density; 

      Double * valArr=new Double[2];
      valArr[0]=val;
      valArr[1]=val;
      // 2 - number of nodes in element
      ptalgsys_->AddElementRHS(valArr,connecth.get(),2,as_sysid_);

      delete [] valArr;
      delete [] ptCoordNodes;
    }    
  } // end of part: standing waves      

}

void Acoustic2dPDE::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, 
				   const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::SolveStepTrans" << std::endl;
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
      SetupMatrices(level,ptBCs);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
      ptalgsys_->ResetRHS(as_sysid_);
      ComputeRHS(lasttimecalc_,ptBCs);
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      ptalgsys_->ResetRHS(as_sysid_);
      ptalgsys_->ResetMatrix(0,0,1);
      ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
      ComputeRHS(lasttimecalc_,ptBCs);
    }
  else
    {
      update = 1;
      job    = 3;
      ptalgsys_->ResetRHS(as_sysid_);
      ComputeRHS(lasttimecalc_,ptBCs);
    };

  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Integer i;
  for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    sol_[i]=ptsol[i];

  if (InfoPrint)
    (*infofile) << "maxnode:" <<  ptgrid_->GetMaxnumnodes(level) << std::endl;

  // calculation of derivatives of solution
  CalcDerSol(); 
}

void Acoustic2dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::WriteResultsInFile" << std::endl;
#endif

  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"vp");
      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"vp_1der");
      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"vp_2der");
    }
  else
    {
      OutFile_->WriteSolution(sol_,laststepcalc_,lasttimecalc_,"fluid potential");
      OutFile_->WriteSolution(sol_der1_,laststepcalc_,lasttimecalc_,"fluid potential, 1st deriv., ");
      OutFile_->WriteSolution(sol_der2_,laststepcalc_,lasttimecalc_,"fluid potential, 2nd deriv., ");
    }

}

void Acoustic2dPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Acoustic2dPDE::CalcParameters" << std::endl;
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

void Acoustic2dPDE :: CalcDerSol()
{
#ifdef TRACE
  (*trace) << " entering  Acoustic2dPDE :: CalcDerSol() " << std::endl;
#endif

  sol_der2_=(sol_ - sol_old_)*a0_ - (sol_der1_old_)*a2_ - sol_der2_old_*a3_;
  sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;
}

void Acoustic2dPDE::SaveSolAsPrevStep()
{
#ifdef TRACE
  (*trace) << " entering  Acoustic2dPDE::SaveSolAsPrevStep() " << std::endl;
#endif

  sol_old_=sol_;
  sol_der1_old_=sol_der1_;
  sol_der2_old_=sol_der2_;  
}


void Acoustic2dPDE ::CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff, Vector<Double> & coeffdamp)
{
  if (!MatFile_) Error("You didn't specialize material file. Check your config-file.");

  coeffmass.Resize(subdoms_.size());
  coeffstiff.Resize(subdoms_.size());
  if (with_absBCs_) coeffdamp.Resize(subdoms_.size());

  Integer i,matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      // read density and compress with material number matnum
      Double density, compress;
      MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");

      coeffmass[i]  = density*density/compress;
      coeffstiff[i] = density;
      if (with_absBCs_)
      coeffdamp[i]  = density/((sqrt(compress/density)));
    }
}

}


