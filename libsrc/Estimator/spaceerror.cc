#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "grid.hh"
#include "spaceerror.hh"
 
#include "interface_piles.hh"
#include "forms_header.hh"
#include "linsystem.hh"

namespace CoupledField
{

template<Integer Dim>
SpaceErrorEstimator<Dim>::SpaceErrorEstimator(Grid*aptGrid)
{
  ptGrid_=aptGrid;
}

template<>
void SpaceErrorEstimator<2>::CalcErrorMap(const Vector<Double> * sol,
 std::vector<std::string> &subdoms, Grid * ptgrid,
 Vector<Double> & errorMap, Vector<Double> & gradSPRElemL2norm)
{
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimation::CalcErrorMap" << std::endl;
#endif
  
  // loop over elements
  Vector<Double> result;
  Integer level=0;
  std::vector<Elem*> elemssd;
  Integer isd,j;
  std::vector<Elem*> *neighbours;
  //vector, where we store global number of nodes for which result is obtained in fnc RecoveryProcedure4ElemsPatch
  std::vector<Integer> locationsInResult;

  Integer dim;
  dim=2;
 
  //vector, where we store SPR gradient
  Vector<Double> * SPRgrad=new Vector<Double>[dim];
  //vector, where we store number of summered value of SPR at node. this vector is needed for averaging procedure of SPRgrad
  Vector<Integer> * numberAvergVals=new Vector<Integer>[dim];

  Integer idm;
  for (idm=0;idm<dim;idm++) {
    SPRgrad[idm].Resize(sol->size());
    numberAvergVals[idm].Resize(sol->size());
  }
  
  // computation of SPR gradient
  Integer icmp;
  for (icmp=0; icmp<dim; icmp++) { // loop over components of gradient
  for (isd=0; isd<subdoms.size(); isd++) // loop over all subdomains
    {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid->GetElemSD(elemssd,subdoms[isd],level);
       
      for (j=0; j< elemssd.size(); j++) // loop over elements of subdomain
	{
	  //  std::vector<Elem*> * neighbours;
	  // in: j - noOfElem, suddoms_[isd] - name of subdomain;
	  // out: pt to vector with neughbors-elements
	  neighbours=ptgrid->GetNeighboursOfElem(j,subdoms[isd]);

	  // add element in list of neighbors to get a full patch of elements
	  neighbours->push_back(elemssd[j]);

	  // in: elemssd[j] - list with elements of this subdomain
	  RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,*sol,icmp,result,locationsInResult);

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
  for (igr=0; igr<sol->size();igr++)
    SPRgrad[icmp][igr]/=numberAvergVals[icmp][igr];

  } // end of loop over components of gradient

  // delete 
  delete [] numberAvergVals;				 
  
 //  cout << " SPR grad: " << std::endl;
//   Integer igr;
//   for (igr=0; igr<sol_.size();igr++)
//     cout << SPRgrad[0][igr] << " " <<  SPRgrad[1][igr] <<std::endl;

  // calculation of error for each element

  //resize arrays
  Integer maxelem=ptgrid->GetMaxnumElem(level,subdoms);
  errorMap.Resize(maxelem);
  gradSPRElemL2norm.Resize(maxelem);
  

    Integer counterElems=0;
    Integer iel;
    for (isd=0; isd<subdoms.size(); isd++) // loop over all subdomains
      {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid->GetElemSD(elemssd,subdoms[isd],level);
       
       // loop over elements of subdomain
      for (iel=0; iel< elemssd.size(); iel++, counterElems++)
	{
	  Double errEl;
	  Double normGradSPR;
	  
	  CalcErrorForElem(elemssd[iel],SPRgrad, errEl, normGradSPR, sol, ptgrid);
	    
	  gradSPRElemL2norm[counterElems]=normGradSPR;
	  errorMap[counterElems]=errEl;

	} // end of loop over elems in subdomain
    } // end of loop over subdomains

 
#ifdef TRACE
  (*trace) << "leaving SpaceErrorEstimator::CalcErrorMap" << std::endl;
#endif    
}

template<>
void SpaceErrorEstimator<3>::CalcErrorMap(const Vector<Double> * sol,
 std::vector<std::string> &subdoms, Grid * ptgrid,
 Vector<Double> & errorMap, Vector<Double> & gradSPRElemL2norm)
{
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimation::CalcErrorMap" << std::endl;
#endif

      // loop over elements
  Vector<Double> result;
  Integer level=0;
  std::vector<Elem*> elemssd;
  Integer isd,j;
  std::vector<Elem*> *neighbours;
  //vector, where we store global number of nodes for which result is obtained in fnc RecoveryProcedure4ElemsPatch
  std::vector<Integer> locationsInResult;

  Integer dim;
  dim=3;
 
  //vector, where we store SPR gradient
  Vector<Double> * SPRgrad=new Vector<Double>[dim];
  //vector, where we store number of summered value of SPR at node. this vector is needed for averaging procedure of SPRgrad
  Vector<Integer> * numberAvergVals=new Vector<Integer>[dim];

  Integer idm;
  for (idm=0;idm<dim;idm++) {
    SPRgrad[idm].Resize(sol->size());
    numberAvergVals[idm].Resize(sol->size());
  }
  
  // computation of SPR gradient
  Integer icmp;
  for (icmp=0; icmp<dim; icmp++) { // loop over components of gradient
  for (isd=0; isd<subdoms.size(); isd++) // loop over all subdomains
    {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid->GetElemSD(elemssd,subdoms[isd],level);
       
      for (j=0; j< elemssd.size(); j++) // loop over elements of subdomain
	{
	  //  std::vector<Elem*> * neighbours;
	  // in: j - noOfElem, suddoms_[isd] - name of subdomain;
	  // out: pt to vector with neughbors-elements
	  neighbours=ptgrid->GetNeighboursOfElem(j,subdoms[isd]);

	  // add element in list of neighbors to get a full patch of elements
	  neighbours->push_back(elemssd[j]);

	  // in: elemssd[j] - list with elements of this subdomain
	  RecoveryProcedure4ElemsPatch(*neighbours,ptgrid,*sol,icmp,result,locationsInResult);

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
  for (igr=0; igr<sol->size();igr++)
    SPRgrad[icmp][igr]/=numberAvergVals[icmp][igr];

  } // end of loop over components of gradient

  // delete 
  delete [] numberAvergVals;				 
  
  // calculation of error for each element

  //resize arrays
  Integer maxelem=ptgrid->GetMaxnumElem(level,subdoms);
  errorMap.Resize(maxelem);
  gradSPRElemL2norm.Resize(maxelem);
  

    Integer counterElems=0;
    Integer iel;
    for (isd=0; isd<subdoms.size(); isd++) // loop over all subdomains
      {      
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid->GetElemSD(elemssd,subdoms[isd],level);
       
       // loop over elements of subdomain
      for (iel=0; iel< elemssd.size(); iel++, counterElems++)
	{
	  Double errEl;
	  Double normGradSPR;
	  
	  CalcErrorForElem(elemssd[iel],SPRgrad, errEl, normGradSPR, sol, ptgrid);
	    
	  gradSPRElemL2norm[counterElems]=normGradSPR;
	  errorMap[counterElems]=errEl;

	} // end of loop over elems in subdomain
    } // end of loop over subdomains

 

#ifdef TRACE
  (*trace) << "leaving SpaceErrorEstimator::CalcErrorMap" << std::endl;
#endif     
}

template<Integer dim>
void SpaceErrorEstimator<dim>::CalcErrorForElem(const Elem* elem,
 const Vector<Double>* SPRgrad,  Double & error, Double & normGradSPR,
 const Vector<Double>* sol, Grid * ptgrid)
{
#ifdef TRACE
  (*trace) << "entering SpaceErrorEstimator::CalcErrorForElem" << std::endl;
#endif

  Integer level=0;

  BaseElem * ptElem=elem->ptElem;

  Vector<Integer> connectVec=elem->connect;
  Integer nnodes=connectVec.size();
  Point<dim> * ptCoord=new Point<dim>[nnodes];
  ptgrid->GetCoordNodesElem(connectVec,ptCoord,level);
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
      Jacobian<dim> J;
      Vector<Double> JinvX,JinvY,JinvZ;
      ptElem->CalcJacobian(J,iIP,ptCoord, TRUE);  
      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);
      if (dim==3) J.GetJinvZ(JinvZ);

      // calculation of FEM gradient at intgr. point iIP in local coord.
      Double gradFEM[dim];
      Integer idm;
      for (idm=0; idm<dim; idm++)
      gradFEM[idm]=0;
     
      Integer ish;
      for (ish=0; ish<nnodes; ish++)
	{  
	  valsol=(*sol)[connectVec[ish]-1];

	  Vector<Double> grad;
	  ptElem->GetGradientShFnc(grad,ish+1,iIP);     
	  
	  Vector<Double> glgrad(dim);	  
	  glgrad[0]=JinvX*grad;
	  glgrad[1]=JinvY*grad;
	  if (dim==3) glgrad[2]=JinvZ*grad;

	  for (idm=0; idm<dim; idm++)
	    gradFEM[idm]+=valsol*glgrad[idm];	
	}	 

      // calculation of SPR gradient at intgr. point iIP in local coord. 
      Double gradSPR[dim];
      for (idm=0; idm<dim; idm++)
	gradSPR[idm]=0;
    
      for (ish=0; ish<nnodes; ish++) {	
	for (idm=0; idm<dim; idm++)
	gradSPR[idm]+=SPRgrad[idm][connectVec[ish]-1]*shFncAtIP[ish][iIP];
      }

     std::cout << gradSPR[0] << " " << gradFEM[0] << " y: " << gradSPR[1] << " " << gradFEM[1] <<std::endl;

      // calculation of L2 norm for difference of SPR and FEM gradient,
      // and L2 norm of SPR gradient
     diffGradL2norm=0;
     SPRGradL2norm=0; 
     for (idm=0; idm<dim; idm++) {
      diffGradL2norm+=sqr(gradSPR[idm]-gradFEM[idm]);
      SPRGradL2norm+=sqr(gradSPR[idm]);
     }

     diffGradL2norm=sqrt(diffGradL2norm);
     SPRGradL2norm=sqrt(SPRGradL2norm);

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
  (*trace) << "leaving SpaceErrorEstimator::CalcErrorForElem" << std::endl;
#endif
}

// template<>
// void SpaceErrorEstimator<Point3D>::CalcErrorForElem(const Elem* elem, const Vector<Double>* SPRgrad,  Double & error, Double & normGradSPR,
//  const Vector<Double>* sol, Grid * ptgrid)
// {
// #ifdef TRACE
//   (*trace) << "entering SpaceErrorEstimator::CalcErrorForElem" << std::endl;
// #endif

//    Integer level=0;

//   BaseElem * ptElem=elem->ptElem;

//   Vector<Integer> connectVec=elem->connect;
//   Integer nnodes=connectVec.size();
//   Point3D * ptCoord=new Point3D[nnodes];
//   ptgrid->GetCoordNodesElem(connectVec,ptCoord,level);
//   Double valsol;
	  
//   // get value of shape fnc for linear interpolation 
//   Vector<Double> shFncAtIP[4];
//   Integer ish;
//   for (ish=0; ish<nnodes; ish++)
//     shFncAtIP[ish]=ptElem->GetShFncAtIP(ish+1);
	   
//   // get integr. points
//   Integer noIntPnts=ptElem->GetNumIntPoints();
//   Vector<Double> * intWeights=ptElem->GetIntWeights();

//   error=0;
//   normGradSPR=0;
  
//   //do a loop over int. points in element
//   Integer iIP;
//   for (iIP=0; iIP<noIntPnts; iIP++)
//     {   
//       Double diffGradL2norm=0;
//       Double SPRGradL2norm=0;
   
//       // calc of Jacobian 
//       Jacobian<Point3D> J;
//       Vector<Double> JinvX,JinvY,JinvZ;
//       ptElem->CalcJacobian(J,iIP,ptCoord, TRUE);  
//       J.GetJinvX(JinvX);
//       J.GetJinvY(JinvY);
//       J.GetJinvZ(JinvZ);

//       // calculation of FEM gradient at intgr. point iIP in local coord.
//       Double gradFEM[3];
//       gradFEM[0]=0;
//       gradFEM[1]=0;
//       gradFEM[2]=0;
//       Integer ish;
//       for (ish=0; ish<nnodes; ish++)
// 	{  
// 	  valsol=(*sol)[connectVec[ish]-1];

// 	  Vector<Double> grad;
// 	  ptElem->GetGradientShFnc(grad,ish+1,iIP);     
	  
// 	  Vector<Double> glgrad(3);	  
// 	  glgrad[0]=JinvX*grad;
// 	  glgrad[1]=JinvY*grad;
// 	  glgrad[2]=JinvZ*grad;

// 	  gradFEM[0]+=valsol*glgrad[0];
// 	  gradFEM[1]+=valsol*glgrad[1];	
// 	  gradFEM[2]+=valsol*glgrad[2];	
// 	}	 

//       // calculation of SPR gradient at intgr. point iIP in local coord. 
//       Double gradSPR[3];
//       gradSPR[0]=0;
//       gradSPR[1]=0;
//       gradSPR[2]=0;

//       for (ish=0; ish<nnodes; ish++) {		
// 	gradSPR[0]+=SPRgrad[0][connectVec[ish]-1]*shFncAtIP[ish][iIP];
// 	gradSPR[1]+=SPRgrad[1][connectVec[ish]-1]*shFncAtIP[ish][iIP];
// 	gradSPR[2]+=SPRgrad[2][connectVec[ish]-1]*shFncAtIP[ish][iIP];
//       }

//       std::cout << gradSPR[0] << " " << gradFEM[0] << 
// 	" y: " << gradSPR[1] << " " << gradFEM[1] << 
// 	" z: " <<  gradSPR[2] << " " << gradFEM[2] << std::endl;


//       // calculation of L2 norm for difference of SPR and FEM gradient,
//       // and L2 norm of SPR gradient
//       diffGradL2norm=sqrt(sqr(gradSPR[0]-gradFEM[0])+
// 			  sqr(gradSPR[1]-gradFEM[1])+sqr(gradSPR[2]-gradFEM[2]));
//       SPRGradL2norm=sqrt(sqr(gradSPR[0])+sqr(gradSPR[1])+sqr(gradSPR[2]));

//       if (intWeights) {
// 	error+=diffGradL2norm*J.detJ*(*intWeights)[iIP];
// 	normGradSPR+=SPRGradL2norm*(*intWeights)[iIP]*J.detJ;
// 	}
// 	else {
// 	error+=diffGradL2norm*J.detJ;
// 	normGradSPR+=SPRGradL2norm*J.detJ;	
// 	}

//     } // end of loop over integration points

// #ifdef TRACE
//   (*trace) << "leaving Elecst3dPDE::CalcErrorForElem" << std::endl;
// #endif
  
// }

template <Integer Dim>
void SpaceErrorEstimator<Dim>::ComputeRHS4RecoverySol(Grid * ptgrid,  const std::vector<std::string> & subdoms, const Integer as_sysid, AbstractAlgebraicSys * ptalgsys, const Vector<Double> & sol)
  {
#ifdef TRACE
    (*trace) << "entering Elecst2dPDE::ComputeRHS4RecoverySol" << std::endl;
#endif

    Integer level=0;
    BaseElem * ptEl;
    Vector<Integer> connecth;
    std::vector<Elem*> elemssd;

    Point<2> * ptCoord;

    Integer isd,j;
    for (isd=0; isd<subdoms.size(); isd++) // loop over all subdomains
      {
	// get vector of Elem of subdomain with color: subdoms_[isd]
	ptgrid->GetElemSD(elemssd,subdoms[isd],level);

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
	    ptgrid->GetCoordNodesElem(connecth,ptCoord,level); 

	    Vector<Double> valFnc;  // get val of sol at nodes of the element
	    valFnc.Resize(elsize);
	    Integer in;
	    for (in=0; in<elsize; in++)
	      valFnc[in]=sol[connecth[in]-1];

	    Vector<Double> vecRHS; /* calculation of elem matrix:
				      int_{element}(ShapeFnc * Sol) */
	    //!!!! component
	    Integer component=0;
	    rhs->CalcElemVector4InterpolatedFnc(ptCoord,component,valFnc,vecRHS);
	    ptalgsys->AddElementRHS(vecRHS.get(),connecth.get(),elsize,as_sysid);

	    delete rhs;
	    delete [] ptCoord;
	  }
      } 
  }

// template <Integer Dim>
// void SpaceErrorEstimator<Dim>::CalcError4EachElem(const Vector<Double>&sol, const Vector<Double>& sol_SPR, Grid * ptgrid, const std::vector<std::string> & subdoms, Vector<Double> & err4elems)
// {
// #ifdef TRACE
//   (*trace) << "entering SpaceErrorEstimator::CalcError4EachElem " << std::endl;
// #endif
  
//    Integer level=0;
//    BaseElem * ptEl;
//    Vector<Integer> connection;
//    std::vector<Elem*> elemssd;

//    Point2D * ptCoord;

//    Integer isd,j;
//    for (isd=0; isd<subdoms.size(); isd++) // loop over all subdomains
//      {
//        // get vector of Elem of subdomain with color: subdoms[isd]
//        ptgrid->GetElemSD(elemssd,subdoms[isd],level);
       
//        err4elems.Resize(elemssd.size());
//        for (j=0; j< elemssd.size(); j++) // loop over elements of subdomain
// 	 {
// 	   connection=elemssd[j]->connect;

// 	   Double result=0, aux;
// 	   Integer ine;	  
// 	   for (ine=0; ine<connection.size(); ine++) {
// 	     aux=sol[connection[ine]-1]-sol_SPR[connection[ine]-1];
// 	     result+=aux*aux;
// 	   }

// 	     err4elems[j]=sqrt(result);	   
// 	 }
//      }
// }

template <Integer Dim>
  void SpaceErrorEstimator<Dim>::RecoveryProcedure4Elem(Elem * aptElem, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result)
  {
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimator::RecoveryProcedure4Elem " << std::endl;
#endif
    Vector<Integer> vec_connect; // auxilialary vector, store connection of element
    Matrix<Double> elemmat; // sys matrix
    Vector<Double> vecRHS; // rhs vector
    Integer level=0;

    // system matrix
    BaseElem * ptBaseElem=aptElem->ptElem;
    BaseForm<2> * bilinear_mass  = new MassInt<2>(ptBaseElem,1);

    vec_connect=aptElem->connect;
    Integer elsize=vec_connect.size();
    Point<2> * ptCoord=new Point<2>[elsize];
 
    ptgrid->GetCoordNodesElem(vec_connect,ptCoord,level); 
    bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

    // rhs
    BaseForm<2> * rhs=new LinearForm<2>(ptBaseElem);
    Vector<Double> valFnc;  // get val of sol at nodes of the element
    valFnc.Resize(elsize);
    Integer in;
    for (in=0; in<elsize; in++)
      valFnc[in]=sol[vec_connect[in]-1];
    rhs->CalcElemVector4InterpolatedFnc(ptCoord,aComponent,valFnc,vecRHS);

    Double tol=1e-6;
    LinSystem<Double,Matrix<Double> > LinSys(tol);
    LinSys.SetSysMatrix(elemmat);
    LinSys.SetRHS(vecRHS);
    LinSys.CG(100,Jacobi);
    LinSys.GetSolution(result);

    delete bilinear_mass;
    delete rhs;
    delete [] ptCoord;
  }

template <Integer Dim>
void SpaceErrorEstimator<Dim>::RecoveryProcedure4ElemsPatch(const std::vector<Elem*> &Elems, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result,std::vector<Integer>&locations)
{
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimator::RecoveryProcedure4ElemsPatch " << std::endl;
#endif
    
    Vector<Integer> vec_connect; // auxilialary vector, stores connection of element
    Matrix<Double> elemmat; // sys matrix for element
    Vector<Double> vecRHS; // rhs vector for element
    Integer counter=0; // counter nodes in element-patch
    Matrix<Double> SysMatrix; // system matrix for elements-patch 
    Vector<Double> RHSVector; // rhs vector for elements-patch
    Integer level=0;
  
    // define dimension of system matrix and rhs 
    Integer nnodes=CalcNumberOfNodesInPatch(Elems);

    SysMatrix.Resize(nnodes,nnodes);
    RHSVector.Resize(nnodes);
// in this vector we store global numbers of nodes that occure in the patch.
// number of position is ,accordingly, local number in patch
    std::vector<Integer>  mapGlbNumNodes2LocInPatch;
    mapGlbNumNodes2LocInPatch.resize(nnodes);

    Integer iel;
    for(iel=0; iel<Elems.size(); iel++) {   // loop over elements in patch 

    // element system matrix
    BaseElem * ptBaseElem=Elems[iel]->ptElem;
    BaseForm<2> * bilinear_mass  = new MassInt<2>(ptBaseElem,1);

    vec_connect=Elems[iel]->connect;
    Integer elsize=vec_connect.size();
    Point<2> * ptCoord=new Point<2>[elsize];
  
    ptgrid->GetCoordNodesElem(vec_connect,ptCoord,level); 
    bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

    // element rhs
    BaseForm<2> * rhs=new LinearForm<2>(ptBaseElem);
    Vector<Double> valFnc;  // get val of sol at nodes of the element
    valFnc.Resize(elsize);
    Integer in;
    for (in=0; in<elsize; in++)
      valFnc[in]=sol[vec_connect[in]-1];
    rhs->CalcElemVector4InterpolatedFnc(ptCoord,aComponent,valFnc,vecRHS);

    // do renumbering of nodes for this element for assembling procedure
    Boolean UseOldNumber;
   
    std::vector<Integer> assembleConnect(vec_connect.size());
    Integer ivc,imp;
    // do a loop over connection vector 
    for (ivc=0; ivc<vec_connect.size(); ivc++) {    
      UseOldNumber=FALSE;

      // loop over vector with global nodes for previous elements
      for (imp=0; imp<mapGlbNumNodes2LocInPatch.size(); imp++) {
	// check that this node is not new
	if (mapGlbNumNodes2LocInPatch[imp] == vec_connect[ivc]) {
	  // if node is new, then we use previous local number in assembling
	  assembleConnect[ivc] = imp;      
	  UseOldNumber=TRUE;
	}
      }

      if (!UseOldNumber) {
	assembleConnect[ivc]=counter;
	mapGlbNumNodes2LocInPatch[counter]=vec_connect[ivc];
	counter++;
      }
    }

    // assemble sysmatrix
    Integer irow,icln;
    Integer ii, iii;
    for (ii=0; ii<elsize; ii++)
      for (iii=0; iii<elsize; iii++)
	{
	  irow=assembleConnect[ii];
	  icln=assembleConnect[iii];
	  SysMatrix.Add(irow,icln,elemmat[ii][iii]);
	} 
   
    // assemble rhs
    for (ii=0; ii<vecRHS.size();ii++) {
      irow=assembleConnect[ii];
      RHSVector[irow]+=vecRHS[ii];
    }    

    delete bilinear_mass;
    delete rhs;
    delete [] ptCoord;
    
    } // end of loop over elements in patch

    // solve lin system
    Double tol=1e-8;
    LinSystem<Double,Matrix<Double> > LinSys(tol);
    LinSys.SetSysMatrix(SysMatrix);
    LinSys.SetRHS(RHSVector);

    //   LinSys.BiCGSTAB(100,Jacobi);
    LinSys.CG(100,Jacobi);
    LinSys.GetSolution(result);

    locations=mapGlbNumNodes2LocInPatch;
  }  // end of fnc: RecoveryProcedure4ElemsPatch

template <Integer Dim>
void SpaceErrorEstimator<Dim>::RecoveryProcedure4ElemsPatch3d(const std::vector<Elem*> &Elems, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result,std::vector<Integer>&locations)
{
#ifdef TRACE
    (*trace) << "entering SpaceErrorEstimator::RecoveryProcedure4ElemsPatch3d " << std::endl;
#endif
    
    Vector<Integer> vec_connect; // auxilialary vector, stores connection of element
    Matrix<Double> elemmat; // sys matrix for element
    Vector<Double> vecRHS; // rhs vector for element
    Integer counter=0; // counter nodes in element-patch
    Matrix<Double> SysMatrix; // system matrix for elements-patch 
    Vector<Double> RHSVector; // rhs vector for elements-patch
    Integer level=0;
  
    // define dimension of system matrix and rhs 
    Integer nnodes=CalcNumberOfNodesInPatch(Elems);

    SysMatrix.Resize(nnodes,nnodes);
    RHSVector.Resize(nnodes);
// in this vector we store global numbers of nodes that occure in the patch.
// number of position is ,accordingly, local number in patch
    std::vector<Integer>  mapGlbNumNodes2LocInPatch;
    mapGlbNumNodes2LocInPatch.resize(nnodes);

    Integer iel;
    for(iel=0; iel<Elems.size(); iel++) {   // loop over elements in patch 

    // element system matrix
    BaseElem * ptBaseElem=Elems[iel]->ptElem;
    BaseForm<3> * bilinear_mass  = new MassInt<3>(ptBaseElem,1);

    vec_connect=Elems[iel]->connect;
    Integer elsize=vec_connect.size();
    Point<3> * ptCoord=new Point<3>[elsize];
  
    ptgrid->GetCoordNodesElem(vec_connect,ptCoord,level); 
    bilinear_mass->CalcElemMatrix(ptCoord, elemmat);

    // element rhs
    BaseForm<3> * rhs=new LinearForm<3>(ptBaseElem);
    Vector<Double> valFnc;  // get val of sol at nodes of the element
    valFnc.Resize(elsize);
    Integer in;
    for (in=0; in<elsize; in++)
      valFnc[in]=sol[vec_connect[in]-1];
    rhs->CalcElemVector4InterpolatedFnc(ptCoord,aComponent,valFnc,vecRHS);

    // do renumbering of nodes for this element for assembling procedure
    Boolean UseOldNumber;
   
    std::vector<Integer> assembleConnect(vec_connect.size());
    Integer ivc,imp;
    // do a loop over connection vector 
    for (ivc=0; ivc<vec_connect.size(); ivc++) {    
      UseOldNumber=FALSE;

      // loop over vector with global nodes for previous elements
      for (imp=0; imp<mapGlbNumNodes2LocInPatch.size(); imp++) {
	// check that this node is not new
	if (mapGlbNumNodes2LocInPatch[imp] == vec_connect[ivc]) {
	  // if node is new, then we use previous local number in assembling
	  assembleConnect[ivc] = imp;      
	  UseOldNumber=TRUE;
	}
      }

      if (!UseOldNumber) {
	assembleConnect[ivc]=counter;
	mapGlbNumNodes2LocInPatch[counter]=vec_connect[ivc];
	counter++;
      }
    }

    // assemble sysmatrix
    Integer irow,icln;
    Integer ii, iii;
    for (ii=0; ii<elsize; ii++)
      for (iii=0; iii<elsize; iii++)
	{
	  irow=assembleConnect[ii];
	  icln=assembleConnect[iii];
	  SysMatrix.Add(irow,icln,elemmat[ii][iii]);
	} 
   
    // assemble rhs
    for (ii=0; ii<vecRHS.size();ii++) {
      irow=assembleConnect[ii];
      RHSVector[irow]+=vecRHS[ii];
    }    

    delete bilinear_mass;
    delete rhs;
    delete [] ptCoord;
    
    } // end of loop over elements in patch

    // solve lin system
    Double tol=1e-8;
    LinSystem<Double,Matrix<Double> > LinSys(tol);
    LinSys.SetSysMatrix(SysMatrix);
    LinSys.SetRHS(RHSVector);

    //   LinSys.BiCGSTAB(100,Jacobi);
    LinSys.CG(100,Jacobi);
    LinSys.GetSolution(result);

    locations=mapGlbNumNodes2LocInPatch;
  }  // end of fnc: RecoveryProcedure4ElemsPatch

template <Integer Dim>
Integer SpaceErrorEstimator<Dim>::CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dPDE::CalcNumberOfNodesInPatch" << std::endl;
#endif

  Integer iels,ivc,imp;
  Vector<Integer> vec_connect;
  Boolean NewNode;
  std::vector<Integer> globalNumbers; // store global numbers of nodes of element that occurs in patch

  for (iels=0; iels<patch.size(); iels++) // loop over elements in patch
    {
     
      vec_connect=patch[iels]->connect;
      
      for (ivc=0; ivc<vec_connect.size(); ivc++) { 
	 NewNode=TRUE;
	// loop over vector with global nodes for previous elements
	for (imp=0; imp<globalNumbers.size(); imp++) {
	  // check that this node is not new
	  if (globalNumbers[imp] == vec_connect[ivc]) { 
	    NewNode=FALSE;
	  }	 
	}

      if (NewNode) {
	globalNumbers.push_back(vec_connect[ivc]);
      }

      } // end of loop over nodes in element

    } // end of loop over elements in patch

  return globalNumbers.size();
     
}

template <Integer Dim>
void SpaceErrorEstimator<Dim>::KellyError(std::vector<std::string> &sbdoms, const Integer level) 
{
#ifdef TRACE
  (*trace) << "entering SpaceErrorEstimator::KellyError " << std::endl;
#endif

  Integer i,j;
  BaseElem * ptEl;
  Double a;

  for (i=0; i<sbdoms.size(); i++){
    std::vector<Elem*> elemssd;
    ptGrid_->GetElemSD(elemssd,sbdoms[i],level); 

    // loop over faces to calculate error
    // IntegralOverFace_KellyError();

    for (j=0; j< elemssd.size(); j++)
      {
	ptEl=elemssd[j]->ptElem;
	
	// loop over faces of Elem and sum calculated integrals

      }
  }
} // end of KellyError

template <Integer Dim>
Double SpaceErrorEstimator<Dim>::IntegralOverRegularFace_KellyError()
{
  Double result;
  Integer ipoint,ifunc;
  Integer noQPoints, noQFunc; //number of integration points and shape functions 

  // 1. get gradients of the function at the first cell at the quadratire points. we know value of solution and shape functions in the quadrature points.
  // 2. the same for the second cell - neigbor
  // 3. calculate the difference. 
  // 4. take the normal
  // 5. now we have vector, each component of it is value at the quadrature point.
  // 6. take a square of each component of vector and sum up with integration weightes
  // 7. also here should be check is this face at the boundary or not
 
  // void calcNormal2Line(std::vector<Double> & normal,const Point2D a, const Point2D b) --> in tools.cc

  return result;
}

} // end of namespace
