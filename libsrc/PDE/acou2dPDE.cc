#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acou2dPDE.hh" 
#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

namespace CoupledField
{

Acou2dPDE::Acou2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::Acou2dPDE " << std::endl;
#endif

  dofspernode_=1;
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

  s_oldold_.Resize(size_);
  s_oldold_.Init(0);

  conf->get("alpha_NM",alpha_,"Acoustic"); 
  conf->get("beta_NM",beta_,"Acoustic"); 
  conf->get("gamma_NM",gamma_,"Acoustic");

  conf->getsubdompde(subdoms_,"Acoustic");
  ReadBCs("Acoustic");

  without_absBCs_=conf->get_option("without_absorbingBCs");

  WriteErrorMap_=TRUE;
  WriteMarkedElements_=TRUE;
  GridIsRefined_=FALSE;
  preComputeRHS();
}



void Acou2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;
  if (without_absBCs_)
    matrix_factor_[1] = 0.0;
  else
    matrix_factor_[1] = 1.0*a1_; 
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 1.0*a0_;
}

void Acou2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::SpecifyMatrices" << std::endl;
#endif

  matrixtype = 1;     // NOCLASS = 0
                      // RSPARSE = 1
                      // CSPARSE = 2
                      // RBLOCK  = 3
                      // CBLOCK  = 4
                      // RFULL   = 5
                      // CFULL   = 6
                      // MIXED   = 7

  /* matrixsystype: NOTYPE     = 0
                    SYSTEM     = 1
                    STIFFNESS  = 2
                    DAMPING    = 3
                    CONVECTION = 4
                    MASS       = 5
  */

  matrixsystype[0] = 1;   // memory for the system matrix
  matrixsystype[1] = 2;   // memory for the stiffness matrix
  if (without_absBCs_)
  matrixsystype[2] = 0;   // memory for the damping matrix
  else
     matrixsystype[2] = 3;
  matrixsystype[3] = 0;   // memory for the convection matrix
  matrixsystype[4] = 5;   // memory for the mass matrix

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3
                 // VOLUME = 4

  numdofpernode  = 1;
  numdirichlets = 1;
  numconstraints = 0;
}

void Acou2dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::SetupMatrices" << std::endl;
#endif

  Matrix<Double> elemmat;

  // This is a smaller matrix since it is just for linear 1D elements.
  Matrix<Double> elemmatbnd;
  Point<2> * ptCoord; 

  BaseElem * ptEl;

  Integer matrix_stiff=2;
  Integer matrix_mass=5;
  Integer matrix_damp=3;

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

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_stiff);

	  // mass part
	  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

	  elemmat*=coeffm[i];

#ifdef DEBUG
	  (*debug) << "Massmatrix, ElementNumber  " << i << std::endl;

	  (*debug) << elemmat << std::endl;
#endif
      
	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_,matrix_mass);

  
	  delete bilinear_stiff;
	  delete bilinear_mass;
	  delete [] ptCoord;
	}
    }

  // BEGIN DAMPING MATRIX PART
  if (!without_absBCs_) {
  std::vector<Elem*>  DomainBnd;
  conf->getliststr("bnd_for_absBCs",bnd_absBCs_,"Acoustic");
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

      ptalgsys_->PutElemMatAlgSys(elemmatbnd.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_damp);
      delete linear_damp;
      delete [] ptCoord;
    }
  }
      // END DAMPING MATRIX PART

#ifdef TRACE
  (*trace) << "Leaving Acou2dPDE::SetupMatrices" << std::endl;
#endif
}

void Acou2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val,val_tfunc;
  Double time=atime;

  //system matrix: id = 1
  Integer matrix_id = 1;

  val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,level);

  // set dirichlet boundary conditions
  Integer i,j=0;
  std::list<Integer> nodes;
  Integer sizebc=bcs_hd_.size();

  val = 0;
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

	      ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
	    }
	  else
	    { 
#ifdef DEBUG
      (*debug) << " node: " << node << " val: " << val << std::endl;
#endif

	      ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
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

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
            }
          else
            {

#ifdef DEBUG
	      (*debug) << " node: " << node <<  " val: " << val << std::endl;
#endif
              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, matrix_id);
            }
	}
    }

#ifdef TRACE
  (*trace) << "leaving Acou2dPDE::SetBCs" << std::endl;
#endif
}

void Acou2dPDE::preComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::preComputeRHS" << std::endl;
#endif

  SetRHSFnc = FALSE;

  //  if (conf->is_there("rhs_surfaces")) {
  std::string isthererhs;
  conf->get("load_force",isthererhs,"Acoustic");

  if (isthererhs=="yes" ) {
    conf->getliststr("rhs_surfaces",rhs_surfaces_,"Acoustic");
    std::string rhs;
    conf->get("rhs_fnc",rhs,"Acoustic");
    ptRHSFnc_=FncReader(rhs);
    conf->get("rhs_arg_fnc",arg_rhs_,"Acoustic");
    SetRHSFnc=TRUE;
    directionFnc_.resize(2);
    conf->get("direction_rhs_x",directionFnc_[0],"Acoustic");
    conf->get("direction_rhs_y",directionFnc_[1],"Acoustic");
    }    
 
}

void Acou2dPDE::ComputeRHS(const Double atime, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::ComputeRHS" << std::endl;
#endif
  Integer n;
  Integer node;
  Integer matrix_id;
  Vector<Double> coeffMass;
  Vector<Double> coeffDamp;
  std::list<Integer> nodes_hd;
  Integer i;

  matrix_id = 5;
  coeffMass = sol_old_*a0_+sol_der1_old_*a2_+sol_der2_old_*a3_;

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffMass.get());

  // damping matrix part
  if (!without_absBCs_) {

  matrix_id = 3; // Setting id for Damping matrix	
  coeffDamp = -sol_der1_old_-sol_der2_old_*a6_;
  
  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffDamp.get());

  matrix_id = 3; // Setting id for Damping matrix	
  coeffDamp = sol_old_*a1_+sol_der1_old_*a2_*a7_+sol_der2_old_*a7_*a3_;
  
  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id,coeffDamp.get());
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

void Acou2dPDE::ComputeRHS4RecoverySol()
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::ComputeRHS4RecoverySol" << std::endl;
#endif
  Integer level=0;
  BaseElem * ptEl;
  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Point<2> * ptCoord;

  Integer isd,j;
  for (isd=0; isd<subdoms_.size(); isd++) // loop over all subdomains
    {
      // get vector of Elem of subdomain with color: subdoms_[isd]
      ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);

      for (j=0; j< elemssd.size(); j++) // loop over elements of subdomain
	{
	  ptEl=elemssd[j]->ptElem;

	  BaseForm<2> * rhs=new LinearForm<2>(ptEl);
	  Integer ii;

	  Integer elsize=(elemssd[j]->connect).size();
	  connecth.Resize(elsize);
	  for (ii=0; ii<elsize; ii++)
	    connecth[ii]=(elemssd[j]->connect)[ii];

	  ptCoord=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level); 

	  // get val of solution at nodes of the element
	  Vector<Double> valFnc; 
	  valFnc.Resize(elsize);
	  Integer in;
	  for (in=0; in<elsize; in++)
	    valFnc[in]=sol_[connecth[in]-1];

	  Vector<Double> vecRHS; /* calculation of elem matrix:
				    int_{Elem}(ShapeFnc * Sol) */
	  ///!!! it is not clear what it is component
	    Integer component=0;
	    rhs->CalcElemVector4InterpolatedFnc(ptCoord,component,valFnc,vecRHS);
	    ptalgsys_->AddElementRHS(vecRHS.get(),connecth.get(),elsize,as_sysid_);

	  delete rhs;
	  delete [] ptCoord;
	}
    }
}

void Acou2dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  SetupMatrices(level,ptBCs);
  SetBCs(ptBCs,level,update,0);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
}

void Acou2dPDE::SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::SolveStepTrans" << std::endl;
#endif

  // restore initial mesh instead of coarsening
//   if (GridIsRefined_) {
//     ptgrid_->RestoreInitialMesh();
//     GridIsRefined_=FALSE;
//   }

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
  // calculation of error map
  CalcErrorMap();
}

void Acou2dPDE::SolveStepTransNewMesh(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::SolveStepTransNewMesh" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=TRUE;

  Double * ptsol;

  Integer update,job;

  update = 0;
  job = 1;
  SetupMatrices(level);
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);

  if (kstep) {
    ptalgsys_->ResetRHS(as_sysid_);
    ComputeRHS(lasttimecalc_);
  }

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

void Acou2dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::WriteResultsInFile" << std::endl;
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

  if (WriteErrorMap_) {   
     OutFile_->WriteDataOnCell(errorMap_,laststepcalc_,lasttimecalc_,"error_ep"); 
  }
  if (WriteMarkedElements_) {
    OutFile_->WriteDataOnCell(markingElems_,laststepcalc_,lasttimecalc_,"marking_elems"); 
  };
  
}

void Acou2dPDE :: CalcParameters(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::CalcParameters" << std::endl;
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

void Acou2dPDE :: CalcDerSol()
{
#ifdef TRACE
  (*trace) << " entering  Acou2dPDE :: CalcDerSol() " << std::endl;
#endif

  sol_der2_=(sol_ - sol_old_)*a0_ - (sol_der1_old_)*a2_ - sol_der2_old_*a3_;
  sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;
}

void Acou2dPDE::SaveSolAsPrevStep()
{
#ifdef TRACE
  (*trace) << " entering  Acou2dPDE::SaveSolAsPrevStep() " << std::endl;
#endif

  sol_old_=sol_;
  sol_der1_old_=sol_der1_;
  sol_der2_old_=sol_der2_;  
}

Double Acou2dPDE::CalcEnergyNorm()
{
 Double help1, help2;
 help1=ptalgsys_->CalcEnergyNorm(0,0,2,sol_.get());
 help2=ptalgsys_->CalcEnergyNorm(0,0,5,sol_.get());
 
 return sqrt(help1+help2);
}

void Acou2dPDE::CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff, Vector<Double> & coeffdamp)
{
  if (!MatFile_) Error("You didn't specialize material file. Check your config-file.");

  coeffmass.Resize(subdoms_.size());
  coeffstiff.Resize(subdoms_.size());
  if (!without_absBCs_) coeffdamp.Resize(subdoms_.size());

  Integer i,matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      // read density and compress with material number matnum
      Double density, compress;
      MatFile_->ReadDensityAndCompressity(density,compress,matnum,"fluid");

      coeffmass[i]  = density*density/compress;
      coeffstiff[i] = density;
      if (!without_absBCs_)
      coeffdamp[i]  = density/((sqrt(compress/density)));
    }
}

void Acou2dPDE::RestoreSol()
{
#ifdef TRACE
  (*trace) << " entering Acou2dPDE::RestoreSol() " << std::endl;
#endif

  Vector<Double> help;
 
  help=sol_old_; 
  ptgrid_->ProlongSol(sol_old_,help,0);
  sol_old_=help;

  help=sol_der1_old_;
  ptgrid_->ProlongSol(sol_der1_old_,help,0);
  sol_der1_old_=help;
  help=sol_der2_old_;
  ptgrid_->ProlongSol(sol_der2_old_,help,0);
  sol_der2_old_=help;

  Integer sizesol=sol_old_.size();
  sol_.Resize(sizesol);
  sol_der1_.Resize(sizesol);
  sol_der2_.Resize(sizesol);
}

void Acou2dPDE::RecoverySol(Vector<Double> & result)
{
#ifdef TRACE
  (*trace) << " entering Acou2dPDE::RecoverySol() " << std::endl;
#endif
  
  Integer update = 1;
  Integer job    = 1;
  Integer level  = 0;

  Double matrix_factor[4];
  matrix_factor[0] = 0.0;
  matrix_factor[1] = 0.0;
  matrix_factor[2] = 0.0;
  matrix_factor[3] = 1.0;

  ptalgsys_->ResetRHS(as_sysid_);
  ptalgsys_->ResetMatrix(0,0,1);
  
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);
  ComputeRHS4RecoverySol();

  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  Double * ptsol;
  ptsol = ptalgsys_->GetSolution(as_sysid_);

  Integer i;
  for (i=0; i<ptgrid_->GetMaxnumnodes(level); i++)
    result[i]=ptsol[i];   
}

void Acou2dPDE::CalcErrorMap()
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::CalcErrorMap" << std::endl;
#endif

  if (ptError_) ConstructorError();

  // loop over elements
  Vector<Double> result;
  Integer level=0;
  std::vector<Elem*> elemssd;
  Integer isd,j;
  std::vector<Elem*> *neighbours;
  //vector, where we store global number of nodes for which result is obtained in fnc RecoveryProcedure4ElemsPatch
  std::vector<Integer> locationsInResult;

  //vector, where we store SPR gradient
  Vector<Double> * SPRgrad=new Vector<Double>[2];
  //vector, where we store number of summered value of SPR at node. this vector is needed for averaging procedure of SPRgrad
  Vector<Integer> * numberAvergVals=new Vector<Integer>[2];

  Integer idm;
  for (idm=0;idm<2;idm++) {
    SPRgrad[idm].Resize(sol_.size());
    numberAvergVals[idm].Resize(sol_.size());
  }
  
  // computation of SPR gradient
  Integer icmp;
  for (icmp=0; icmp<2; icmp++) { // loop over components of gradient
  for (isd=0; isd<subdoms_.size(); isd++) // loop over all subdomains
    {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
       
      for (j=0; j< elemssd.size(); j++) // loop over elements of subdomain
	{
	  //  std::vector<Elem*> * neighbours;
	  // in: j - noOfElem, suddoms_[isd] - name of subdomain;
	  // out: pt to vector with neughbors-elements
	  neighbours=ptgrid_->GetNeighboursOfElem(j,subdoms_[isd]);

	  // add element in list of neighbors to get a full patch of elements
	  neighbours->push_back(elemssd[j]);

	  // in: elemssd[j] - list with elements of this subdomain
	  ptError_->RecoveryProcedure4ElemsPatch(*neighbours,ptgrid_,sol_,icmp,result,locationsInResult);

	  // form arrays for averaging values of SPRgrad over nodes
	  Integer ires;
	  for (ires=0;ires<result.size();ires++) {
	    SPRgrad[icmp][locationsInResult[ires]-1]+=result[ires];
	    numberAvergVals[icmp][locationsInResult[ires]-1]++;
	  }
	 
	} // loop over subdomains
    } // loop over component of gradient

  //average values of SPRgrad over nodes
  Integer igr;
  for (igr=0; igr<sol_.size();igr++)
    SPRgrad[icmp][igr]/=numberAvergVals[icmp][igr];

  } // end of loop over components of gradient

  // delete 
  delete [] numberAvergVals;				 

  // calculation of error for each element

  //resize arrays
  Integer maxelem=ptgrid_->GetMaxnumElem(level,subdoms_);
  errorMap_.Resize(maxelem);
  gradSPRElemL2norm_.Resize(maxelem);
  

    Integer counterElems=0;
    Integer iel;
    for (isd=0; isd<subdoms_.size(); isd++) // loop over all subdomains
      {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
       
       // loop over elements of subdomain
      for (iel=0; iel< elemssd.size(); iel++, counterElems++)
	{
	  Double errEl;
	  Double normGradSPR;
	  
	  CalcErrorForElem(elemssd[iel],SPRgrad, errEl, normGradSPR);
	    
	  gradSPRElemL2norm_[counterElems]=normGradSPR;
	  errorMap_[counterElems]=errEl;

	} // end of loop over elems in subdomain
    } // end of loop over subdomains

#ifdef TRACE
  (*trace) << "leaving Acou2dPDE::CalcErrorMap" << std::endl;
#endif
}

void Acou2dPDE::CalcErrorForElem(const Elem* elem, const Vector<Double>* SPRgrad,
 Double & error, Double & normGradSPR)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::CalcErrorForElem" << std::endl;
#endif

  Integer level=0;

  BaseElem * ptElem=elem->ptElem;

  Vector<Integer> connectVec=elem->connect;
  Integer nnodes=connectVec.size();
  Point<2> * ptCoord=new Point<2>[nnodes];
  ptgrid_->GetCoordNodesElem(connectVec,ptCoord,level);
  Double valsol;
	  
  // get value of shape fnc for linear interpolation 
  Vector<Double> shFncAtIP[4];
  Integer ish;
  for (ish=0; ish<nnodes; ish++)
    shFncAtIP[ish]=ptElem->GetShFncAtIP(ish+1);
	   
  // get integr. points
  Integer noIntPnts=ptElem->GetNumIntPoints();
  Vector<Double> * intWeights=ptElem->GetIntWeights();

  error=0;
  normGradSPR=0;
  
  //do a loop over int. points in element
  Integer iIP;
  for (iIP=0; iIP<noIntPnts; iIP++)
    {   
      Double diffGradL2norm=0;
      Double SPRGradL2norm=0;
   
      // calc of Jacobian 
      Jacobian<2> J;
      Vector<Double> JinvX,JinvY;
      ptElem->CalcJacobian(J,iIP,ptCoord, TRUE);  
      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);

      // calculation of FEM gradient at intgr. point iIP in local coord.
      Double gradFEM[2];
      gradFEM[0]=0;
      gradFEM[1]=0;
      Integer ish;
      for (ish=0; ish<nnodes; ish++)
	{  
	  valsol=sol_[connectVec[ish]-1];

	  Vector<Double> grad;
	  ptElem->GetGradientShFnc(grad,ish+1,iIP);     
	  
	  Vector<Double> glgrad(2);	  
	  glgrad[0]=JinvX*grad;
	  glgrad[1]=JinvY*grad;

	  gradFEM[0]+=valsol*glgrad[0];
	  gradFEM[1]+=valsol*glgrad[1];		
	}	 

      // calculation of SPR gradient at intgr. point iIP in local coord. 
      Double gradSPR[2];
      gradSPR[0]=0;
      gradSPR[1]=0;
      for (ish=0; ish<nnodes; ish++) {		
	gradSPR[0]+=SPRgrad[0][connectVec[ish]-1]*shFncAtIP[ish][iIP];
	gradSPR[1]+=SPRgrad[1][connectVec[ish]-1]*shFncAtIP[ish][iIP];
      }

      // calculation of L2 norm for difference of SPR and FEM gradient,
      // and L2 norm of SPR gradient
      diffGradL2norm=sqrt(sqr(gradSPR[0]-gradFEM[0])+sqr(gradSPR[1]-gradFEM[1]));
      SPRGradL2norm=sqrt(sqr(gradSPR[0])+sqr(gradSPR[1]));

      if (intWeights) {
	error+=diffGradL2norm*J.detJ*(*intWeights)[iIP];
	normGradSPR+=SPRGradL2norm*(*intWeights)[iIP]*J.detJ;
	}
	else {
	error+=diffGradL2norm*J.detJ;
	normGradSPR+=SPRGradL2norm*J.detJ;	
	}

    } // end of loop over integration points

#ifdef TRACE
  (*trace) << "leaving Acou2dPDE::CalcErrorForElem" << std::endl;
#endif
}

Boolean Acou2dPDE::TestError()
{
#ifdef TRACE
  (*trace) << "entering Acou2d::TestError" << std::endl;
#endif

  CalcErrorMap();
 
  normError_=0;
  Integer itm;
  for (itm=0; itm<gradSPRElemL2norm_.size(); itm++) {
    normError_+=gradSPRElemL2norm_[itm];
  }

  Double sumErrorMap=0;
  Integer i;
  for (i=0; i<errorMap_.size();i++) {
    sumErrorMap+=errorMap_[i];
  }

  Double totalError=sumErrorMap/normError_;
  
  relativeErrorMap_=errorMap_/normError_;
 
  conf->get("error_tolerance",errorTol_,"SpaceAdaptivity");

  if (InfoPrint)
  (*infofile) << " Total error: " << totalError << " diffError " << sumErrorMap << " norm " << normError_ << std::endl;

  if (totalError>errorTol_) return TRUE;
  else return FALSE;
}

void Acou2dPDE::RefineMesh()
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::RefineMesh" << std::endl;
#endif

   Integer level=0;

   // get max num elements for the domain,on which we have the equation
   Integer numElems=ptgrid_->GetMaxnumElem(level,subdoms_);

   // get pointer to array with elements
   std::vector<Elem*> elemssd;
   ptgrid_->GetElemSD(elemssd,subdoms_[0],level);

   // choose error tolerance for each of element
   // according to area of element
   Double sumArea=0;
   Vector<Double> areaElems(numElems);
   Integer iem;
   for (iem=0; iem<numElems; iem++) {
        areaElems[iem]=ptgrid_->CalcAreaElem(elemssd[iem]);
	sumArea+= areaElems[iem];
   }

   Vector<Double> errorLocalTol(numElems);
   for (iem=0; iem<numElems; iem++) 
     errorLocalTol[iem]=errorTol_*areaElems[iem]/sumArea;

  Double errLocTol=errorTol_/numElems; 

  relativeErrorMap_=errorMap_/normError_;
  markingElems_.Resize(numElems);

 //  if (InfoPrint)
//   (*infofile) << " relative error map: " << endl << relativeErrorMap_ << endl << " tolerance: " << errLocTol << endl;

  for (iem=0; iem<numElems; iem++)
    if (relativeErrorMap_[iem]>errorLocalTol[iem]) { 
      elemssd[iem]->refinementFlag=TRUE;
      markingElems_[iem]=1.;
    }
    else { elemssd[iem]->refinementFlag=FALSE;}

  ptgrid_->Refine();
  GridIsRefined_=TRUE;
}

//In this fnc we delete old pointer to Error-object & create new
void Acou2dPDE::ConstructorError()
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::ConstructorError" << std::endl;
#endif

  if (ptError_) delete ptError_;
  
  ptError_=new SpaceErrorEstimator<2>();
  ptError_->Init(this);

}

TimeErrorEstimator * Acou2dPDE::CreatePtTimeError()
{
  return ptTimeError_=new AcousticTimeErrorEstimator(this);
}

Acou2dPDE::~Acou2dPDE()
{
  if (ptTimeError_) delete ptTimeError_; 
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !! old stuff, which is not used !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void Acou2dPDE :: CalculationDerivativesSol(const Boolean Recalc)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::CalculationDerivativesSol" << std::endl;
#endif
  
  if (Recalc) std::cout << laststepcalc_ << std::endl;

  if (!Recalc)
  { sol_der2_old_=sol_der2_;
    sol_der1_old_=sol_der1_; 
  }

  if (!Recalc)
  sol_der2_=(sol_-sol_old_)*a0_-sol_der1_old_*a2_-sol_der2_old_*a3_;
  else
  sol_der2_=(sol_-s_oldold_)*a0_-sol_der1_old_*a2_-sol_der2_old_*a3_;

  sol_der1_=sol_der1_old_+sol_der2_old_*a6_+sol_der2_*a7_;


  if (!Recalc) { s_oldold_=sol_old_;}
  sol_old_=sol_;
}

void Acou2dPDE::CalcThirdDerivateFromEquation(Vector<Double> & result)
{
#ifdef TRACE
  (*trace) << "entering Acou2dPDE::CalcThirdDerivateFromEquation" << std::endl;
#endif

  Double * ptres;

  Double mat_factor[4];
  Integer update,job;

  mat_factor[0] = 0.0;
  mat_factor[1] = 0.0;
  mat_factor[2] = 0.0;
  mat_factor[3] = 1.0;

  job    = 1;

  ptalgsys_->ResetRHS(as_sysid_);
  ptalgsys_->ResetMatrix(0,0,1);
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,mat_factor);

  Integer matrix_id;
  Vector<Double> coeffMass;

  matrix_id = 2;
  Vector<Double> help=-sol_der1_;

  ptalgsys_->UpdateRHS(as_sysid_,as_sysid_,matrix_id, help.get());

  update=1; // use boundary condition at this time step

  //  SetBCs(ptBCs,level,update,lasttimecalc_);
  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);
  ptres = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Integer level=0; //   !!!!!!!!!!!!!!

  Vector<Double> transres(ptgrid_->GetMaxnumnodes(level), ptres);
  result=transres;

  (*infofile) << " result " << transres << std::endl;
}

} // end of namespace

// void Acou2dPDE::SetupMatrices(const Integer level)
// {
// #ifdef TRACE
//   (*trace) << "entering Acou2dPDE::SetupMatrices" << std::endl;
// #endif
 
//   Vector<Double> coeffm, coeffst;
//   CalcCoeff(coeffm, coeffst);

//  Integer i;
//  for (i=0; i<subdoms_.size(); i++)
// {
//   PutElemMatInAlgSys putelmatalgsys(ptalgsys_,ptgrid_,coeffm[i],coeffst[i],as_sysid_,level);

//   ptgrid_->forEachElemSd(putelmatalgsys,subdoms_[i]);
// }
// }

/*
void PutElemMatInAlgSys::operator()(Elem t)
{
  Matrix<Double> elemmat;

  BaseForm<2> * bilinear_stiff = new LaplaceInt<2>(t.ptElem,1);
  BaseForm<2> * bilinear_mass  = new MassInt<2>(t.ptElem,1);

  Point<2> * ptCoord=new Point<2>[t.connect.size()];
  ptgrid_->GetCoordNodesElem(t.connect,ptCoord,level_);

  // stiffness part
  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);

  elemmat*=coeffs_;

#ifdef DEBUG
      (*debug) << "Stiffnessmatrix, ElementNumber  "  << std::endl;

      (*debug) << elemmat << std::endl;
#endif

  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), t.connect.get(), t.connect.size(), sysid_, sysid_, matrix_stiff_);

      // mass part
  bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

  elemmat*=coeffm_;

#ifdef DEBUG
      (*debug) << "Massmatrix, ElementNumber  "  << std::endl;

      (*debug) << elemmat << std::endl;
#endif

  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), t.connect.get(), t.connect.size(), sysid_, sysid_,matrix_mass_);

  delete bilinear_stiff;
  delete bilinear_mass;
  delete [] ptCoord;
}
*/


