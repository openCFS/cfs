#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "elecst3dPDE.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>
 
namespace CoupledField
{

Elecst3dPDE::Elecst3dPDE(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut)
:BasePDE(ptalgsys,ptMaterial,aptFileType,aptOut,aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::Electst3dPDE " << std::endl;
#endif

  dofspernode_=1;
  ptgrid_=aptgrid;

  size_=ptgrid_->GetMaxnumnodes(0);
  sol_.Resize(size_);
  sol_.Init(0);

  conf->getsubdompde(subdoms_,"Elecst3d");
  ReadBCs("Elecst3d");

  WriteErrorMap_=conf->get_option("write_error_map");

  calcElecField_=FALSE;

  ptError_=NULL;
}

void Elecst3dPDE::SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, Integer &numeqcoarse, Double &coarsealpha )
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SpecifySolver" << std::endl;
#endif

  conf->get("eps",eps,"Elecst2d"); // relative accuracy in the precond. energy
  conf->get("dampiter",dampiter,"Elecst2d"); // damping parameter for Jacobi, SSOR
  conf->get("maxnumit",maxnumit,"Elecst2d"); // max number of iterations
  conf->get("solvertype",solvertype,"Elecst2d"); // Richardson or CG
  conf->get("precondtype", precondtype, "Elecst2d"); //ID or MG
  conf->get("numeqcoarse",numeqcoarse,"Elecst2d"); // number of equation for coarsing
  conf->get("coarsealpha",coarsealpha,"Elecst2d"); // coarsing parameter for AMG

} 


void Elecst3dPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 1.0;
  matrix_factor_[1] = 0.0;
  matrix_factor_[2] = 0.0;
  matrix_factor_[3] = 0.0;
}

void Elecst3dPDE::SpecifyMatrices(Integer &matrixtype, Integer * matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets,
Integer &numconstraints)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SpecifyMatrices" << std::endl;
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
  matrixsystype[4] = 0;   // memory for the mass matrix

  graphtype = 1; // NOGRAPH = 0
                 // NODE   = 1
                 // EDGE   = 2
                 // FACE   = 3
                 // VOLUME = 4

  numdofpernode  = 1;
  numdirichlets = 1;
  numconstraints = 0;
}

void Elecst3dPDE::SetupMatrices(const Integer level, BCs * ptBCs)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;
  Point<3> * ptCoord;

  BaseElem * ptElem;

  Integer matrix_stiff=2;

  Vector<Double> coeffst;
  CalcCoeff(coeffst);  

  Vector<Integer> connecth;
  std::vector<Elem*> elemssd;

  Integer i, j;
  for (i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm<3> * bilinear_stiff = new LaplaceInt<3>(ptElem,1);

	  connecth=elemssd[j]->connect;

	  ptCoord=new Point<3>[connecth.size()];
	  ptgrid_->GetCoordNodesElem(connecth,ptCoord,level);

	  // stiffness part
	  bilinear_stiff->CalcElemMatrix(ptCoord, elemmat);
	  elemmat*=coeffst[i];

	  if (InfoPrint)
	    (*infofile) << elemmat << std::endl;

#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;

	  (*debug) << elemmat << std::endl;
#endif

	  ptalgsys_->PutElemMatAlgSys(elemmat.getinarray(), connecth.get(), connecth.size(), as_sysid_, as_sysid_, matrix_stiff);

	  delete bilinear_stiff;
	  delete [] ptCoord;
	}     
    }

#ifdef TRACE
  (*trace) << "Leaving Elecst3dPDE::SetupMatrices" << std::endl;
#endif
}

void Elecst3dPDE::SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, valueTF;

  //system matrix: id = 1
  Integer matrix_id = 1;

  Integer i,j=0;
  std::list<Integer> nodes;

  for (i=0; i<bcs_hd_.size(); i++)
    {  
      nodes=ptBCs->GetNodesLevel(bcs_hd_[i], level);
  
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
            }
	}
    }

  for (i=0; i<bcs_id_.size(); i++)
    {
      nodes=ptBCs->GetNodesLevel(bcs_id_[i], level);

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
            }
	}
    }
  
#ifdef TRACE
  (*trace) << " leaving Elecst3d::ReadBCs " << std::endl;
#endif
 
}

void Elecst3dPDE::SolveStepStatic(BCs * ptBCs, Integer level)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  SetupMatrices(level);
  ptalgsys_->ComputeSysMatrix(as_sysid_,as_sysid_,matrix_factor_);

  SetBCs(ptBCs,level,update,0);

  ptalgsys_->ComputePrecond(job,as_sysid_);
  ptalgsys_->SolveAlgSys(as_sysid_);

  ptsol = ptalgsys_->GetSolution(as_sysid_);

  // save solution
  Vector<Double> transsol(ptgrid_->GetMaxnumnodes(level), ptsol);
  sol_=transsol;

// check
  Integer is;
  for (is=0; is<sol_.size(); is++) {
    Point<3> XY;
    ptgrid_->GetCoordinateNode(is,0,XY);
  }

 if (WriteErrorMap_) {
   ptError_=new SpaceErrorEstimator<3>();
   ptError_->Init(this);
   ptError_->CalcErrorMap(&sol_,subdoms_,ptgrid_,errorMap_,gradSPRElemL2norm_);
 }
  
  if (calcElecField_)
   CalcElectricField();

}

void Elecst3dPDE::CalcElectricField()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::CalcElectricField" << std::endl;
#endif

  Integer level=0;
  Integer maxelem=ptgrid_->GetMaxnumElem(level,subdoms_);

  Integer icp;
  for (icp=0; icp<3; icp++)
    elecFieldAtCenterElem_[icp].Resize(maxelem);
 
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
	  Vector<Double> grad(3),glgrad(3);

	  //for transformation of element to standard
	  Jacobian<3> J;
	  Vector<Double> JinvX, JinvY, JinvZ;

	  BaseElem * ptElem=elemssd[iel]->ptElem;
	  Vector<Integer> connectVec=elemssd[iel]->connect;
	  Integer nnodes=connectVec.size();
	  Point<3> * ptCoord=new Point<3>[nnodes];
	  ptgrid_->GetCoordNodesElem(connectVec,ptCoord,level);
	  Double valsol;

	  Vector<Double> gradElectricFieldAtCenterElem;
	  gradElectricFieldAtCenterElem.Resize(3);
	      
	  ptElem->CalcJacobianAtCenter(J,ptCoord,TRUE);
	  J.GetJinvX(JinvX);
	  J.GetJinvY(JinvY);
	  J.GetJinvZ(JinvZ);
	      
	  // do a loop over shape functions	  
	  Integer ish;
	  for (ish=0; ish<nnodes; ish++)
	    {  
	      valsol=sol_[connectVec[ish]-1];
		  
	      ptElem->GetGradientShFncAtCenter(grad,ish+1);

	      glgrad[0]=JinvX*grad;
	      glgrad[1]=JinvY*grad;
	      glgrad[2]=JinvZ*grad;

	      gradElectricFieldAtCenterElem[0]+=valsol*glgrad[0];
	      gradElectricFieldAtCenterElem[1]+=valsol*glgrad[1];
	      gradElectricFieldAtCenterElem[2]+=valsol*glgrad[2];
	    }	 

	  for (icp=0; icp<3; icp++)
	    elecFieldAtCenterElem_[icp][counterElems]=gradElectricFieldAtCenterElem[icp];
	   
	  //  absValueElectricField_[icp][counterElems]=gradElectricFieldAtCenterElem.normL2();

	  delete [] ptCoord;
	} // end of loop over elements of subdomains

      // test
   //   Integer ivr;
   //   for (ivr=0; ivr<counterElems; ivr++) {
//	for (icp=0; icp<3; icp++)
//	  std::cout << elecFieldAtCenterElem_[icp][ivr] << " ";
//      std::cout << std::endl;
//      }

      } // end of loop over subdomain
}

void Elecst3dPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;

  // write results
  if (OutFile_->IsGMV())
    OutFile_->WriteNodeSolution(sol_,step,time,"electric_potential");
  else
    OutFile_->WriteNodeSolution(sol_,step,time,"electric potential");

  // write error map
   if (WriteErrorMap_) {   
     OutFile_->WriteElemSolution(errorMap_,step,time,"error_ep"); 
  }

   // write electric field
   if (calcElecField_) 
       OutFile_->WriteElemSolution(elecFieldAtCenterElem_,step,time,"elec_field",3);

}

Double Elecst3dPDE::CalcEnergyNorm()
{
  Double help1;
  help1=ptalgsys_->CalcEnergyNorm(0,0,2,sol_.get());
 
  return sqrt(help1);
}

void Elecst3dPDE::CalcCoeff(Vector<Double> & coeff)
{
  if (!MatFile_) Error("You didn't specialize material file. Use option -m");

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

void Elecst3dPDE::CalcErrorMap()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::CalcErrorMap" << std::endl;
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
  Vector<Double> * SPRgrad=new Vector<Double>[3];
  //vector, where we store number of summered value of SPR at node. this vector is needed for averaging procedure of SPRgrad
  Vector<Integer> * numberAvergVals=new Vector<Integer>[3];

  Integer idm;
  for (idm=0;idm<3;idm++) {
    SPRgrad[idm].Resize(sol_.size());
    numberAvergVals[idm].Resize(sol_.size());
  }
  
  // computation of SPR gradient
  Integer icmp;
  for (icmp=0; icmp<3; icmp++) { // loop over components of gradient
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
  
//  std::cout << " SPR grad: " << std::endl;
//  Integer igr;
//  for (igr=0; igr<sol_.size();igr++)
//    std::cout << SPRgrad[0][igr] << " " <<  SPRgrad[1][igr] << " " <<  SPRgrad[2][igr] << std::endl;

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

	  std::cout << " error for elem: " << errEl << std::endl;

	} // end of loop over elems in subdomain
    } // end of loop over subdomains

#ifdef TRACE
  (*trace) << "leaving Elecst3dPDE::CalcErrorMap" << std::endl;
#endif
}

//In this fnc we delete old pointer to Error-object & create new
void Elecst3dPDE::ConstructorError()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::ConstructorError" << std::endl;
#endif

  if (ptError_) delete ptError_;
  
  ptError_=new SpaceErrorEstimator<3>();
  ptError_->Init(this);

}

void Elecst3dPDE::CalcErrorForElem(const Elem* elem, const Vector<Double>* SPRgrad,
 Double & error, Double & normGradSPR)
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::CalcErrorForElem" << std::endl;
#endif

  Integer level=0;

  BaseElem * ptElem=elem->ptElem;

  Vector<Integer> connectVec=elem->connect;
  Integer nnodes=connectVec.size();
  Point<3> * ptCoord=new Point<3>[nnodes];
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
      Jacobian<3> J;
      Vector<Double> JinvX,JinvY,JinvZ;
      ptElem->CalcJacobian(J,iIP,ptCoord, TRUE);  
      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);
      J.GetJinvZ(JinvZ);

      // calculation of FEM gradient at intgr. point iIP in local coord.
      Double gradFEM[3];
      gradFEM[0]=0;
      gradFEM[1]=0;
      gradFEM[2]=0;
      Integer ish;
      for (ish=0; ish<nnodes; ish++)
	{  
	  valsol=sol_[connectVec[ish]-1];

	  Vector<Double> grad;
	  ptElem->GetGradientShFnc(grad,ish+1,iIP);     
	  
	  Vector<Double> glgrad(3);	  
	  glgrad[0]=JinvX*grad;
	  glgrad[1]=JinvY*grad;
	  glgrad[2]=JinvZ*grad;

	  gradFEM[0]+=valsol*glgrad[0];
	  gradFEM[1]+=valsol*glgrad[1];	
	  gradFEM[2]+=valsol*glgrad[2];	
	}	 

      // calculation of SPR gradient at intgr. point iIP in local coord. 
      Double gradSPR[3];
      gradSPR[0]=0;
      gradSPR[1]=0;
      gradSPR[2]=0;

      for (ish=0; ish<nnodes; ish++) {		
	gradSPR[0]+=SPRgrad[0][connectVec[ish]-1]*shFncAtIP[ish][iIP];
	gradSPR[1]+=SPRgrad[1][connectVec[ish]-1]*shFncAtIP[ish][iIP];
	gradSPR[2]+=SPRgrad[2][connectVec[ish]-1]*shFncAtIP[ish][iIP];
      }

      std::cout << gradSPR[0] << " " << gradFEM[0] << 
	" y: " << gradSPR[1] << " " << gradFEM[1] << 
	" z: " <<  gradSPR[2] << " " << gradFEM[2] << std::endl;


      // calculation of L2 norm for difference of SPR and FEM gradient,
      // and L2 norm of SPR gradient
      diffGradL2norm=sqrt(sqr(gradSPR[0]-gradFEM[0])+
			  sqr(gradSPR[1]-gradFEM[1])+sqr(gradSPR[2]-gradFEM[2]));
      SPRGradL2norm=sqrt(sqr(gradSPR[0])+sqr(gradSPR[1])+sqr(gradSPR[2]));

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
  (*trace) << "leaving Elecst3dPDE::CalcErrorForElem" << std::endl;
#endif
}

Boolean Elecst3dPDE::TestError()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::TestError" << std::endl;
#endif

  if (!ptError_) ConstructorError();

    ptError_->CalcErrorMap(&sol_,subdoms_,ptgrid_,errorMap_,gradSPRElemL2norm_);
 
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

  if (InfoPrint) {
  (*infofile) << " total Error: " << totalError << " tolerance Error:" << errorTol_ << std::endl;
  }

  std::cerr << " Total error: " << totalError << " diffError " << sumErrorMap << " norm " << normError_ << std::endl;

  if (totalError>errorTol_) return TRUE;
  else return FALSE;
}

void Elecst3dPDE::RefineMesh()
{
#ifdef TRACE
  (*trace) << "entering Elecst3dPDE::RefineMesh" << std::endl;
#endif

   Integer level=0;

   // get max num elements for the domain,on which we have the equation
   Integer numElems=ptgrid_->GetMaxnumElem(level,subdoms_);

   // get pointer to array with elements
   std::vector<Elem*> elemssd;
   ptgrid_->GetElemSD(elemssd,subdoms_[0],level);

//    // choose error tolerance for each of element
//    // according to area of element
//    Double sumArea=0;
//    Vector<Double> areaElems(numElems);
//    Integer iem;
//    for (iem=0; iem<numElems; iem++) {
//         areaElems[iem]=ptgrid_->CalcAreaElem(elemssd[iem]);
// 	sumArea+= areaElems[iem];
//    }

//    Vector<Double> errorLocalTol(numElems);
//    for (iem=0; iem<numElems; iem++) 
//      errorLocalTol[iem]=errorTol_*areaElems[iem]/sumArea;

  Double errLocTol=errorTol_/numElems; 

  relativeErrorMap_=errorMap_/normError_;
  markingElems_.Resize(numElems);

  if (InfoPrint)
  (*infofile) << " relative error map: " << std::endl << relativeErrorMap_ << std::endl << " tolerance: " << errLocTol << std::endl;

 //  Integer i;
//   for (i=0; i< relativeErrorMap_.size(); i++) {
//     cout <<  relativeErrorMap_[i] << " " << errorLocalTol[i] << " " << errLocTol << endl;
//   }

  Integer iem;
  for (iem=0; iem<numElems; iem++)
    if (relativeErrorMap_[iem]>errLocTol) { 
      elemssd[iem]->refinementFlag=TRUE;
      markingElems_[iem]=1.;
    }
    else { elemssd[iem]->refinementFlag=FALSE;}

  ptgrid_->RefineUniform();
}

 //! write additional info (marked elements, relative error) to files with mesh
void Elecst3dPDE::PrintMeshesInfo(WriteResults * ptMeshes)
{
   
  ptMeshes->WriteElemSolution(relativeErrorMap_,0,0,"relative_error");
  ptMeshes->WriteElemSolution(markingElems_,0,0,"marked_elems");
 
}

Elecst3dPDE::~Elecst3dPDE()
{
 ;
}

} // end of namespace


