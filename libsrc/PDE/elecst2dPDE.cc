#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elecst2dPDE.hh"
#include "outUnverg.hh"
#include "forms_header.hh"
#include "linsystem.hh"
#include "spaceerror.hh" 

namespace CoupledField
{

Elecst2dPDE::Elecst2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, 
			 FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::Electst2dPDE " << std::endl;
#endif

  dofspernode_=1;
  ptgrid_=aptgrid;

  size_=ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);

  conf->getsubdompde(subdoms_,"Elecst2d");
  ReadBCs("Elecst2d");

  errorTol_ = 0;
  NeedMassMatrix_=TRUE;

  WriteErrorMap_=conf->get_option("write_error_map");

  SolverCFS_=FALSE;
  
  ptError_=NULL;
  solGrad_=NULL;

  WriteErrorMap_=conf->get_option("write_error_map");
  AbsValueElectricField_=conf->get_option("calc_electr_field");

}

void Elecst2dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, Integer &numeqcoarse, Double &coarsealpha )
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Elecst2d"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Elecst2d"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Elecst2d"); // max number of iterations
  conf->get("solvertype",solvertype,"Elecst2d"); // Richardson or CG
  conf->get("precondtype", precondtype, "Elecst2d"); //ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Elecst2d"); // number of equation for coarsing
  conf->get("coarsealpha",coarsealpha,"Elecst2d"); // coarsing parameter for AMG

} 

void Elecst2dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SetMatrixFactors" << std::endl;
#endif
  
  matrix_factor_[0] = 1.0; 
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 0.0;
}

void Elecst2dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SpecifyMatrices" << std::endl;
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
  matrixsystype[2] = 0;   // memory for the damping matrix
  matrixsystype[3] = 0;   // memory for the convection matrix
  if (NeedMassMatrix_) 
    matrixsystype[4] = 5;
  else matrixsystype[4] = 0;   // memory for the mass matrix

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3
                 // VOLUME = 4

  numdofpernode  = 1;
  numdirichlets = 1;
  numconstraints = 0;
}

void Elecst2dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Point<2> * ptCoord;

  BaseElem * ptElem;

  Integer matrix_stiff = 2;
  Integer matrix_mass = 5;

  if (InfoPrint)
    (*infofile) << " ------------------------- Element matrices --------------- " << std::endl;

  Vector<Double> coeffst;
  
  CalcCoeff(coeffst);  

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  Vector<Integer> connecth;
 
	  ptElem=elemssd[j]->ptElem;

	  BaseForm<2> * bilinear_stiff = new LaplaceInt<2>(ptElem,1);
	  BaseForm<2> * bilinear_mass = NULL;

	  if (NeedMassMatrix_)
	    bilinear_mass = new MassInt<2>(ptElem,1);

	  connecth=elemssd[j]->connect;

	  ptCoord=new Point<2>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];
	  
#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;

	  (*debug) << elemmat << std::endl;
#endif
	  if (InfoPrint)
	    (*infofile) << elemmat << std::endl;

	 
	  if (SolverCFS_) {
	  Integer elsize=connecth.size();
	  Integer irow,icln;
	  Integer ii, iii;
	    for (ii=0; ii<elsize; ii++)
	      for (iii=0; iii<elsize; iii++)
		{
		  irow=connecth[ii]-1;
		  icln=connecth[iii]-1;
		  sysmat_.Add(irow,icln,elemmat[ii][iii]);
		} 
	  }


	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_stiff);

	  //mass
	  if (NeedMassMatrix_){
	
	    bilinear_mass->CalcElemMatrix(ptCoord, elemmat);
	    
	    ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_mass);	    
	  }

	  delete bilinear_stiff;
	  if (bilinear_mass)   delete bilinear_mass;
	  delete [] ptCoord;
	}  
      
    }

#ifdef TRACE
  (*trace) << "Leaving Elecst2dPDE::SetupMatrices" << std::endl;
#endif
}

void Elecst2dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SetBCs" << std::endl;
#endif


  Double BigConst;
    if (SolverCFS_) {
  Integer matsize=sysmat_.getSize(),i;
  Double max=sysmat_(0,0);
  for (i=0;i<matsize;i++)
    if (sysmat_(i,i)>max) max=sysmat_(i,i);
 

   BigConst=(10e+10)*max;
    }
 
  Integer node;
  Double val, valueTF;

  //system matrix: id = 1
  Integer matrix_id = 1;

  Integer j=0;
  std::list<Integer> nodes;

  if (InfoPrint)
    (*infofile) << " ---------------- Dirichle boundary condition -------------" << std::endl;

  Integer i;
  for (i=0; i<bcs_hd_.size(); i++)
    {  
      nodes=ptBCs->GetNodesLevel(bcs_hd_[i]);
  
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  val=0; 
          if (update==1)
            {

	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
            }
          else
            {

	       if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,
					as_sysid_, matrix_id);
	      if (SolverCFS_) {
	      vecrhs_[node-1]+=BigConst*val;
 	      sysmat_(node-1,node-1)+=BigConst;
	      }
	     
            }
	}
    }

  for (i=0; i<bcs_id_.size(); i++)
    {
      nodes=ptBCs->GetNodesLevel(bcs_id_[i]);

      val=val_id_[i];

      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
          if (update==1)
            {	     

	       if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

              ptalgsys_->SetBCDirichletUpdate(j+1, node, val, dofspernode_, as_sysid_, as_sysid_, matrix_id);
            }
          else
            {
	      
	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;
	      
	      ptalgsys_->SetBCDirichlet(j+1, node, val, dofspernode_, as_sysid_,as_sysid_, matrix_id);
	   
		if (SolverCFS_) {
		vecrhs_[node-1]+=BigConst*val;
		sysmat_(node-1,node-1)+=BigConst;
		}
	
            }
	}
    }
  
#ifdef TRACE
  (*trace) << " leaving Elecst2d::SetBCs " << std::endl;
#endif 
}

void Elecst2dPDE::ComputeRHS()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::ComputeRHS" << std::endl;
#endif
  Vector<Double> elemvec;

  Integer level=0;
  // read fnc from conf-file for RHS

  std::vector<Elem*> elemssd;
  // loop over elements 
  Integer isd;
  for (isd=0; isd<subdoms_.size(); isd++) // loop over all subdomains
    {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);

      Integer j;
      for (j=0; j< elemssd.size(); j++) // loop over elements of subdomain
	{
	  BaseElem * ptElem=elemssd[j]->ptElem;

	  BaseForm<2> * linearRHS=new LinearForm<2>(ptElem);
	  Vector<Integer> connectVec=elemssd[j]->connect;
	  Integer nnodes=connectVec.size();
	  Point<2> * ptCoord=new Point<2>[nnodes];
	  ptgrid_->GetCoordNodesElem(connectVec,ptCoord,level);

	  linearRHS->CalcElemVector(ptCoord,elemvec);
 
	  elemvec*=-1; // multiply by -1, since the sysmat we assembly without -1

	  delete linearRHS;
	  delete [] ptCoord;
	  
	  ptalgsys_->AddElementRHS(elemvec.get(),connectVec.get(),connectVec.size(),as_sysid_);

	  if (SolverCFS_) {
	    Integer nd,sz;
	    sz=connectVec.size();
	    Integer ind;
	    for (ind=0; ind<sz; ind++)
	      vecrhs_[connectVec[ind]-1]+=elemvec[ind];
	  }
	  
	} // end of loop over elements
    } // end of loop over subdomains 
}

void Elecst2dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  if (SolverCFS_) {
  Integer size=ptgrid_->GetMaxnumnodes(level);
  sysmat_.Resize(size,size);
  vecrhs_.Resize(size);
  }
 
  SetupMatrices(level);
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);

  ComputeRHS();
  SetBCs(ptBCs,level,update,0);

  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);

  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);
  sol_=transsol;

  if (SolverCFS_) {
  Double tol=1e-20;
  LinSystem<Double,Matrix<Double> > LinSys(tol);
 
  LinSys.SetSysMatrix(sysmat_);
  LinSys.SetRHS(vecrhs_);

  Vector<Double> result;
  LinSys.CG(1000,Jacobi);
  LinSys.GetSolution(result); 
  
   Integer is;
   for (is=0; is<sol_.size(); is++) {
     Point<2> XY;
     ptgrid_->GetCoordinateNode(is,0,XY);
	std::cout << sol_[is] << " " << result[is] << " " << XY[0]*200 << std::endl;
   }

  } // end of fnc for solver CFS++ - LinAlg

  // check
   Integer is;
  for (is=0; is<sol_.size(); is++) {
     Point<2> XY;
     ptgrid_->GetCoordinateNode(is,0,XY);

  std::cout << sol_[is] << " exact: " << 200*XY[0] << std::endl;

//     //  cout << sol_[is] << " exact: " << (1-XY.x*XY.x-XY.y*XY.y)*0.25 << " x:" << XY.x << " y:" << XY.y << std::endl;


//   // sol_[is]= XY.x*XY.y*(1-XY.x)*(1-XY.y)*(1/tan(0.8*((XY.x+XY.y)/sqrt(2)-20)));
  }

 if  (AbsValueElectricField_) CalcElectricField();

 if (WriteErrorMap_) {
   ptError_=new SpaceErrorEstimator<2>();
   ptError_->Init(this);
   ptError_->CalcErrorMap(&sol_,subdoms_,ptgrid_,errorMap_,gradSPRElemL2norm_);
 }

}

void Elecst2dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::WriteResultsInFile" << std::endl;
#endif
 
  Integer step=0;
  Double time=0;
  if (OutFile_->IsGMV()) 
    OutFile_->WriteSolution(sol_,step,time,"electric_potential");
  else 
    OutFile_->WriteSolution(sol_,step,time,"electric potential");

  if (WriteErrorMap_) {   
     OutFile_->WriteDataOnCell(errorMap_,step,time,"error_ep"); 
  }

  if (AbsValueElectricField_) {
      OutFile_->WriteDataOnCell(absValueElectricField_,step,time,"elecfield");
  }
}

Double Elecst2dPDE::CalcEnergyNorm()
{
  Double help1;
  help1=ptalgsys_->CalcEnergyNorm(0,0,2,sol_.get());
 
  return sqrt(help1);
}

void Elecst2dPDE::CalcCoeff(Vector<Double> & coeff)
{
#ifdef TRACE
  (*trace) << " entering Elecst2dPDE::CalcCoeff " <<std::endl;
#endif

  if (!MatFile_) Error("You didn't specialize material file.");

  coeff.Resize(subdoms_.size());

  Integer i, matnum;
  for (i=0; i<subdoms_.size(); i++)
    {
      conf->get(subdoms_[i],matnum,"list_subdomains");

      Double dielectr;
      MatFile_->ReadDielectricTerms(dielectr,matnum); 

      coeff[i]=dielectr;    
    }

}

Boolean Elecst2dPDE::TestError()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::TestError" << std::endl;
#endif

  if (!ptError_) ConstructorError();

  // calculation of error map
  ptError_->CalcErrorMap(&sol_,subdoms_,ptgrid_,errorMap_,gradSPRElemL2norm_);
 
  // calculation of ||grad_SPR||_L2
  normError_=0;
  Integer itm;
  for (itm=0; itm<gradSPRElemL2norm_.size(); itm++) {
    normError_+=gradSPRElemL2norm_[itm];
  }

  // calculation of ||grad_SPR-grad_FEM||_L2
  Double sumErrorMap=0;
  Integer i;
  for (i=0; i<errorMap_.size();i++) {
    sumErrorMap+=errorMap_[i];
  }

    std::cout << errorMap_ << std::endl;
 std::cout << sumErrorMap << " " << normError_ << std::endl;

  Double totalError=sumErrorMap/normError_;
  
  relativeErrorMap_=errorMap_/normError_;

  if (!errorTol_)
    conf->get("error_tolerance",errorTol_,"SpaceAdaptivity");


  if (InfoPrint) {
  (*infofile) << " total Error: " << totalError << " tolerance:" << errorTol_ << std::endl;
  }

  std::cerr << " Total error: " << totalError << " tolerance " << errorTol_ << " diffError " << sumErrorMap << " norm " << normError_ << std::endl;

  if (totalError>errorTol_) return TRUE;
  else return FALSE;
}

//In this fnc we delete old pointer to Error-object & create new
void Elecst2dPDE::ConstructorError()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::ConstructorError" << std::endl;
#endif

  if (ptError_) delete ptError_;
  
  ptError_=new SpaceErrorEstimator<2>();
  ptError_->Init(this);
}

void Elecst2dPDE::RefineMesh()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::RefineMesh" << std::endl;
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

   // calculation of error tolerance for each element
   Vector<Double> errorLocalTol(numElems);
   for (iem=0; iem<numElems; iem++) 
     errorLocalTol[iem]=errorTol_*areaElems[iem]/sumArea;

   // form array from 1 and 0 for marking elements
   markingElems_.Resize(numElems);

  // mark element for refiniment
  for (iem=0; iem<numElems; iem++)
    if (relativeErrorMap_[iem]>errorLocalTol[iem]) { 
      elemssd[iem]->refinementFlag=TRUE;
      markingElems_[iem]=1.;
    }
    else { elemssd[iem]->refinementFlag=FALSE;}

  // do refinement
  ptgrid_->Refine();
}

void Elecst2dPDE::CalcGradSolElemIP(const Elem * element, const Integer level, Vector<Double>*gradElIP)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::CalcGradSolElemIP" << std::endl;
#endif

  if (solGrad_) delete [] solGrad_;
  
  solGrad_=new Vector<Double>[2];

  // 
  Vector<Double> grad(2),glgrad(2);
  // for transformation of element to standard
  Jacobian<2> J;
  Vector<Double> JinvX, JinvY;

  BaseElem * ptElem=element->ptElem;
  Vector<Integer> connectVec=element->connect;
  Integer nnodes=connectVec.size();
  Integer nIntgPnts=ptElem->GetNumIntPoints();
  Point<2> * ptCoord=new Point<2>[nnodes];
  ptgrid_->GetCoordNodesElem(connectVec,ptCoord,level);
  Double valsol;

 
  gradElIP[0].Resize(nIntgPnts);
  gradElIP[1].Resize(nIntgPnts);

  //do a loop over integration points in element
  Integer iIP;
  for (iIP=0; iIP<nIntgPnts; iIP++)
    {
	      
      ptElem->CalcJacobian(J,iIP,ptCoord,TRUE);
      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);
	      
      // do a loop over shape functions	  
      Integer ish;
      for (ish=0; ish<nnodes; ish++)
	{  
	  valsol=sol_[connectVec[ish]-1];
		  
	  ptElem->GetGradientShFnc(grad,ish+1,iIP);

	  glgrad[0]=JinvX*grad;
	  glgrad[1]=JinvY*grad;

	  gradElIP[0][iIP]+=valsol*glgrad[0];
	  gradElIP[1][iIP]+=valsol*glgrad[1];		
	}	 

    } // end of loop over integration points of the element

  delete [] ptCoord;
#ifdef TRACE
  (*trace) << "leaving Elecst2dPDE::CalcGradSolElemIP" << std::endl;
#endif 
}

void Elecst2dPDE::CalcElectricField()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::CalcElectricField" << std::endl;
#endif
   Integer level=0;
  Integer maxelem=ptgrid_->GetMaxnumElem(level,subdoms_);
  absValueElectricField_.Resize(maxelem);

 
  std::vector<Elem*> elemssd;
  Integer counterElems=0;
  // loop over all elements
  Integer isd;
   for (isd=0; isd<subdoms_.size(); isd++) // loop over all subdomains
      {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
       
       // loop over elements of subdomain
      Integer iel;
      for (iel=0; iel< elemssd.size(); iel++,counterElems++)
	{
	  // 
	  Vector<Double> grad(2),glgrad(2);

	  //for transformation of element to standard
	  Jacobian<2> J;
	  Vector<Double> JinvX, JinvY;

	  BaseElem * ptElem=elemssd[iel]->ptElem;
	  Vector<Integer> connectVec=elemssd[iel]->connect;
	  Integer nnodes=connectVec.size();
	  Point<2> * ptCoord=new Point<2>[nnodes];
	  ptgrid_->GetCoordNodesElem(connectVec,ptCoord,level);
	  Double valsol;

	  Vector<Double> gradElectricFieldAtCenterElem;
	  gradElectricFieldAtCenterElem.Resize(2);
	      
	  ptElem->CalcJacobianAtCenter(J,ptCoord,TRUE);
	  J.GetJinvX(JinvX);
	  J.GetJinvY(JinvY);
	      
	  // do a loop over shape functions	  
	  Integer ish;
	  for (ish=0; ish<nnodes; ish++)
	    {  
	      valsol=sol_[connectVec[ish]-1];
		  
	      ptElem->GetGradientShFncAtCenter(grad,ish+1);

	      glgrad[0]=JinvX*grad;
	      glgrad[1]=JinvY*grad;

	      gradElectricFieldAtCenterElem[0]+=valsol*glgrad[0];
	      gradElectricFieldAtCenterElem[1]+=valsol*glgrad[1];		
	    }	 

	  absValueElectricField_[counterElems]=gradElectricFieldAtCenterElem.normL2();

	  delete [] ptCoord;
	} // end of loop over elements of subdomains
      } // end of loop over subdomain
}

 //! write additional info (marked elements, relative error) to files with mesh
void Elecst2dPDE::PrintMeshesInfo(WriteResults * ptMeshes)
{
   
  ptMeshes->WriteDataOnCell(relativeErrorMap_,0,0,"relative_error");
  ptMeshes->WriteDataOnCell(markingElems_,0,0,"marked_elems");
 
}

Elecst2dPDE::~Elecst2dPDE()
{
#ifdef TRACE
  (*trace) << " entering Elecst2dPDE::~Elecst2dPDE " << std::endl;
#endif
  if (ptError_) delete ptError_ ;

  if (solGrad_) delete [] solGrad_;
}

} // end of namespace
