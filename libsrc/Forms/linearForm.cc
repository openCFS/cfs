#include <iostream>
#include <fstream>
#include <math.h>

#include "linearForm.hh"


namespace CoupledField
{

  LinearForm::LinearForm(BaseFE * aptelem) : BaseForm(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::LinearForm" << std::endl;
#endif
  }



  LinearForm ::~LinearForm()
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::~LinearForm" << std::endl;
#endif
  }


  void LinearForm::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & Result)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::CalcElemVector" << std::endl;
#endif
  }

  LinearEdgeInt::LinearEdgeInt(BaseFE * aptelem, Double aVal, Integer aDirection,
			       std::vector<Double> * aCoilMidPt) 
    : LinearForm(aptelem), val_(aVal), direction_(aDirection)
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::LinearEdgeInt" << std::endl;
#endif
    
    if(aCoilMidPt)
      coilMidPt_ = new std::vector<Double>(*aCoilMidPt);
  }


  LinearEdgeInt ::~LinearEdgeInt()
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::~LinearEdgeInt" << std::endl;
#endif
  }



  void LinearEdgeInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering LinearEdgeInt::CalcElemVector" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrEdges  = ptelem->GetNumEdges();
    const Integer dim      = ptelem->GetDim();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Integer xDir = 0;
    const Integer yDir = 1;
    const Integer zDir = 2;
    
    Double jacDet;  


    // derivation of shape functions after global coordinates 
    Matrix<Double> shapeEdge;
    std::vector<Double> partElemVec;
    partElemVec.resize(nrEdges);


    // set vector to desired size and set all elements to zero    
    elemVec.resize(nrEdges); 
    for (Integer i=0; i<nrEdges; i++)
      elemVec[i] = 0;
    
    std::vector<Double> currentVec(dim, 0);
    // is needed, if coil is curved
    std::vector<Double> globCoord;
    Double len;
    
    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
	// define direction of current
	switch(direction_)
	  {	
	  case 1:
	  case 2:
	  case 3:
	    // if current direction goes parallel to a coordinate axis
	    currentVec[direction_-1] = val_;
	    break;

	  case 4: // current in xy-plane ==> z=0
	    ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	    currentVec[xDir] = -(globCoord[yDir] - (*coilMidPt_)[yDir]);
	    currentVec[yDir] = globCoord[xDir] - (*coilMidPt_)[xDir];
	    currentVec[zDir] = 0;

	    len = L2Norm(currentVec);
	    currentVec *= val_ / len;
	    break;

	  case 5: // current in yz-plane ==> x=0
	    ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	    currentVec[xDir] = 0;
	    currentVec[yDir] = -(globCoord[zDir] - (*coilMidPt_)[zDir]);
	    currentVec[zDir] = globCoord[yDir] - (*coilMidPt_)[yDir];

	    len = L2Norm(currentVec);
	    currentVec *= val_ / len;
	    break;

	  case 6: // current in xz-plane ==> y=0
	    ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	    currentVec[xDir] = globCoord[zDir] - (*coilMidPt_)[zDir];
	    currentVec[yDir] = 0;
	    currentVec[zDir] = -(globCoord[xDir] - (*coilMidPt_)[xDir]);

	    len = L2Norm(currentVec);
	    currentVec *= val_ / len;

	    break;
	    
	  default:
	    std::string errMsg = "Selected current direction with number ";
	    errMsg += direction_ + " not supported!";
	    Error(errMsg.c_str(),__FILE__,__LINE__);
	  }
	
	ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord);

	partElemVec = shapeEdge * currentVec;

	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

	for(Integer i=0; i<partElemVec.size(); i++)
	  elemVec[i] += partElemVec[i] * intWeights[actIntPt-1] * jacDet;	
      }
  

#ifdef DEBUG 
    (*debug) << "CalcElemVector:  "  << std::endl
	     << partElemVec << std::endl
	     << "\n jacDet " << jacDet << std::endl;
#endif


#ifdef TRACE
    (*trace) << "leaving linearEdgeInt::CalcElementVector" << std::endl;
#endif

  }



  // ==================================================================
  // nLinMech
  // ==================================================================


  nLinMech_linFormInt::nLinMech_linFormInt(BaseFE * aptelem, MaterialData & mat) 
    : LinearForm(aptelem), matData_(mat)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::nLinMech_linFormInt" << std::endl;
#endif
  }


  
  nLinMech_linFormInt ::~nLinMech_linFormInt()
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::~nLinMech_linFormInt" << std::endl;
#endif
  }



  void nLinMech_linFormInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::CalcElemVector" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer nrDofs   = getNrDofs();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    std::vector<Double> piolaStressVec;    
    std::vector<Double> partElemVec;

    Matrix<Double> linBMat; 
    Matrix<Double> nonLinBMat; 
    Matrix<Double> dMat;
    Matrix<Double> transpSumB;    // we need transposed of the b-matrices


    // This, as friend defined bilinearform holds the necessary differential operators
    nLinMech3dInt_PiolaStress piolaStressBiform(ptelem, matData_);

    if (!elemDisp_.size_row() || !elemDisp_.size_col()) 
      Error("Undefined displacements! ",__FILE__,__LINE__);

    partElemVec.resize(nrNodes * nrDofs);


   elemVec.resize(nrNodes*nrDofs);
   elemVec *= 0;    // set elems to 0
   


    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
	piolaStressBiform.setActElemDispl(elemDisp_);
	piolaStressBiform.calcPiolaStressVec(piolaStressVec, actIntPt, ptCoord);
	piolaStressBiform.calcNonLinBMat(nonLinBMat, actIntPt, ptCoord);
	piolaStressBiform.calcLinBMat(linBMat, actIntPt, ptCoord);

	nonLinBMat += linBMat;

	nonLinBMat.Transpose(transpSumB);
	
	partElemVec = transpSumB * piolaStressVec;

	Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
	partElemVec *= jacDet * intWeights[actIntPt-1];
	
	elemVec +=  partElemVec;

      }
  
#ifdef TRACE
    (*trace) << "leaving nLinMech_linFormInt::CalcElementVector" << std::endl;
#endif

  }


  ////////////////////////////////////////////////////////////////////////////
  //Members of class LinearFormFlowNoise

LinearFlowNoiseInt::LinearFlowNoiseInt(BaseFE * aptelem) : LinearForm(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::LinearFlowNoiseInt" << std::endl;
#endif
  }

LinearFlowNoiseInt::~LinearFlowNoiseInt()
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::~LinearFlowNoiseInt" << std::endl;
#endif
  }

void LinearFlowNoiseInt::CalcElemVector4Dip(Matrix<Double>& ptCoord,const Vector<Integer> & connecth, Vector<Double> & Result, const std::vector<Double> gradN_x_P)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::CalcElemVector4FlowSrc" << std::endl;
#endif

    Integer l=ptelem->GetNumIntPoints();
    Integer n=ptelem->GetNumNodes();


       
    Integer i,ii;

    Double jacDet;    
//     Double density=1.0;
    Result.Resize(n);
    std::vector<double>  Sf;
    Double multiplier;
    // Dipole Source Term: (w,deltp)onbnd
    std::vector<Double> normal;
    normal.resize(gradN_x_P.size());	

    switch(normal.size())
      {
      case 2:
	{
	  calcNormal2Line_Mat(normal,ptCoord);
	  std::cout<<"normal :"<<normal[0]<<"  "<<normal[1]<<std::endl; 
	  break;
	}
      case 3:
	{
	  calcNormal2Surface_Mat(normal,ptCoord);
	  std::cout<<"normal :"<<normal[0]<<"  "<<normal[1]<<"  "<<normal[2]<<std::endl;
	  break;
	}
      default:
	std::string errMsg = "Incorrect dimension of normal. ";
	Error(errMsg.c_str(),__FILE__,__LINE__);
      }

    for (i=0; i<l; i++)
      {
	jacDet = ptelem->CalcJacobianDetAtIp(i+1, ptCoord);

	ptelem->GetShFncAtIp(Sf, i+1);

	for (ii=0; ii < n; ii++)
	  {
	    multiplier=ScalarMult(normal,gradN_x_P); // n * gradShFnc
	    // Result[ii]+= -length/2.0*multiplier*Sf[ii][i];
	    // 	  std::cout<<"multiplier :"<<multiplier<<std::endl; 
	    // 	  std::cout<<"jacdet :"<<jacDet<<std::endl; 
	    // 	  std::cout<<"Sf[ii] :"<<Sf[ii]<<std::endl; 

	    Result[ii]+= -jacDet*multiplier*Sf[ii];
	    // 	  std::cout<<"Result["<<ii<<"]="<< Result[ii]<<std::endl; 
	  }

      }

  } // end of method

  //Doing it for acoustic RHS (quadrupole)

void LinearFlowNoiseInt::CalcElemVector4Quad(Matrix<Double>& ptCoord,const Vector<Integer> & connecth,const Matrix<Double> & FlowData, Vector<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4Quad" << std::endl;
#endif

  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();
  //   Vector<Double> * help=new Vector<Double>[n];
  Matrix<Double> xiDx;       
  Double jacDet;
  Integer actInt,ii;
  //   Vector<Double> JinvX, JinvY, JinvZ;
  //   Jacobian<dim> J;
  Double density=1.0;
  Result.Resize(n);
  //   Vector<Double> * Sf=new Vector<Double> [n];
  std::vector<double>  Sf;

  Double dTxxdx1part,dTxxdx2part,dTxydx1part,dTxydx2part;
  Double dTyxdy1part,dTyxdy2part,dTyydy1part,dTyydy2part;


  
  //   for (i=0; i < n; i++)
  //     Sf[i]=ptelem->GetShFncAtIP(i+1);


  for (actInt=1; actInt<=l; actInt++)
    {

      ptelem->GetShFncAtIp(Sf, actInt);

      //          ptelem->CalcJacobian(J,i,ptCoord);
      //           J.GetJinvX(JinvX);
      //           J.GetJinvY(JinvY);
      // Resetting derivatives of tensor
      dTxxdx1part=0;
      dTxxdx2part=0;
      dTxydx1part=0;
      dTxydx2part=0;
      dTyxdy1part=0;
      dTyxdy2part=0;
      dTyydy1part=0;
      dTyydy2part=0;
      // This work only for quad1 elements since we have values
      // for ux and uy only at the corners!
      // Here we compute the derivatives of the Lighthill's tensor needed in the quadrupole term
      // at the ith integration point


      ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet);

      for (ii=0; ii<n; ii++)
	{
	  //               ptelem->GetGradientShFnc(help[ii],ii+1,i);
	      
	  //dTxxdx1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]); 
	  //dTxxdx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];
	  dTxxdx1part+=(FlowData[1][connecth[ii]-1]*Sf[ii]); 
	  dTxxdx2part+=xiDx[ii][0]*FlowData[1][connecth[ii]-1];

	  //dTxydx1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
	  //dTxydx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];
	  dTxydx1part+=(FlowData[2][connecth[ii]-1]*Sf[ii]);
	  dTxydx2part+=xiDx[ii][0]*FlowData[1][connecth[ii]-1];

	  //dTyxdy1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
	  //dTyxdy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
	  dTyxdy1part+=(FlowData[1][connecth[ii]-1]*Sf[ii]);
	  dTyxdy2part+=xiDx[ii][1]*FlowData[2][connecth[ii]-1];

	  //dTyydy1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
	  //dTyydy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
	  dTyydy1part+=(FlowData[2][connecth[ii]-1]*Sf[ii]);
	  dTyydy2part+=xiDx[ii][1]*FlowData[2][connecth[ii]-1];
	}


      // Setting up the RHS (dw/di,dTij/di) part
      // 	  for (ii=0; ii < n; ii++) Result[ii]+= density*((help[ii]*JinvX)*(2.*dTxxdx1part*dTxxdx2part+dTxydx1part*dTxydx2part)+(help[ii]*JinvY)*(dTyxdy1part*dTyxdy2part+2.*dTyydy1part*dTyydy2part))*J.detJ;

      for (ii=0; ii < n; ii++) Result[ii]+= density*(xiDx[ii][0]*(2.*dTxxdx1part*dTxxdx2part+dTxydx1part*dTxydx2part)+
						     xiDx[ii][1]*(dTyxdy1part*dTyxdy2part+2.*dTyydy1part*dTyydy2part))*jacDet;
    }




  //28.08, NOCH NICHT UMGESTELLT, IN ZUKUNFT SOLLTE ALLGEMEIN SEIN!!!!
//   if (dim==3)
//     {

//       Double dTxzdx1part, dTxzdx2part, dTyzdy1part, dTyzdy2part, dTzxdz1part, dTzxdz2part, dTzydz1part, dTzydz2part, dTzzdz1part, dTzzdz2part;
//       for (i=0; i<l; i++)
// 	{
// 	  ptelem->CalcJacobian(J,i,ptCoord);  
// 	  J.GetJinvX(JinvX);
// 	  J.GetJinvY(JinvY);
// 	  J.GetJinvZ(JinvZ);
// 	  // Resetting derivatives of tensor
// 	  dTxxdx1part=0;
// 	  dTxxdx2part=0;
// 	  dTxydx1part=0;
// 	  dTxydx2part=0;
// 	  dTxzdx1part=0;
// 	  dTxzdx2part=0;
	  
// 	  dTyxdy1part=0;
// 	  dTyxdy2part=0;
// 	  dTyydy1part=0;
// 	  dTyydy2part=0;
// 	  dTyzdy1part=0;
// 	  dTyzdy2part=0;

// 	  dTzxdz1part=0;
// 	  dTzxdz2part=0;
// 	  dTzydz1part=0;
// 	  dTzydz2part=0;
// 	  dTzzdz1part=0;
// 	  dTzzdz2part=0;

// 	  // This work only for trilinear hexahedral elements since we have values
// 	  // for ux uy and uz only at the corners!
// 	  // Here we compute the derivatives of the Lighthill's tensor needed in the quadrupole term
// 	  // at the ith integration point
// 	  for (ii=0; ii<n; ii++)
// 	    {
// 	      ptelem->GetGradientShFnc(help[ii],ii+1,i);

// 	      dTxxdx1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
// 	      dTxxdx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];
// 	      dTxydx1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
// 	      dTxydx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];
// 	      dTxzdx1part+=(FlowData[3][connecth[ii]-1]*Sf[ii][i]);
// 	      dTxzdx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];

// 	      dTyxdy1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
// 	      dTyxdy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
// 	      dTyydy1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
// 	      dTyydy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
// 	      dTyzdy1part+=(FlowData[3][connecth[ii]-1]*Sf[ii][i]);
// 	      dTyzdy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];

// 	      dTzxdz1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
// 	      dTzxdz2part+=(help[ii]*JinvZ)*FlowData[3][connecth[ii]-1];
// 	      dTzydz1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
// 	      dTzydz2part+=(help[ii]*JinvZ)*FlowData[3][connecth[ii]-1];
// 	      dTzzdz1part+=(FlowData[3][connecth[ii]-1]*Sf[ii][i]);
// 	      dTzzdz2part+=(help[ii]*JinvZ)*FlowData[3][connecth[ii]-1];
// 	    }

// 	  // Setting up the RHS (dw/di,dTij/di) part
// 	  for (ii=0; ii < n; ii++) Result[ii]+= density*((help[ii]*JinvX)*(2.*dTxxdx1part*dTxxdx2part+dTxydx1part*dTxydx2part+dTxzdx1part*dTxzdx2part)+(help[ii]*JinvY)*(dTyxdy1part*dTyxdy2part+2.*dTyydy1part*dTyydy2part+dTyzdy1part*dTyzdy2part)+(help[ii]*JinvZ)*(dTzxdz1part*dTzxdz2part+dTzydz1part*dTzydz2part+2.*dTzzdz1part*dTzzdz2part))*J.detJ;
// 	}   
//     } // end of if (dim==3)

} // end of method


} // end of namespace
