#include <iostream>
#include <fstream>
#include <math.h>

#include "forms_header.hh"

namespace CoupledField {

LinearForm::LinearForm(BaseFE * aptelem) : BaseForm(aptelem)
{
  ENTER_FCN( "LinearForm::LinearForm" );
}


LinearForm::LinearForm() : BaseForm()
{
  ENTER_FCN( "LinearForm::LinearForm" );
}




LinearForm::~LinearForm()
{
  ENTER_FCN( "LinearForm::~LinearForm" );
}


void LinearForm::CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & Result)
{
  ENTER_FCN( "LinearForm::CalcElemVector" );
}


// =============================================================================
// edge integration
// =============================================================================

LinearEdgeInt::LinearEdgeInt(BaseFE * aptelem, Double aVal, Integer aDirection,
			     Vector<Double> * aCoilMidPt) 
  : LinearForm(aptelem), val_(aVal), direction_(aDirection)
{
  ENTER_FCN( "LinearEdgeInt::LinearEdgeInt" );
  
  if(aCoilMidPt)
    coilMidPt_ = new Vector<Double>(*aCoilMidPt);
}

LinearEdgeInt::~LinearEdgeInt()
{
  ENTER_FCN( "LinearEdgeInt::~LinearEdgeInt" );
}



void LinearEdgeInt::CalcElemVector(Matrix<Double>& ptCoord, 
				   Vector<Double> & elemVec)
{
  ENTER_FCN( "LinearEdgeInt::CalcElemVector" );
  
  const Integer nrIntPts = ptelem->GetNumIntPoints();
  const Integer nrEdges  = ptelem->GetNumEdges();
  const Integer dim      = ptelem->GetDim();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  
  const Integer xDir = 0;
  const Integer yDir = 1;
  const Integer zDir = 2;
  
  Double jacDet;  
  
  
  // derivation of shape functions after global coordinates 
  Matrix<Double> shapeEdge;
  Vector<Double> partElemVec;
  partElemVec.Resize(nrEdges);
  
  
  // set vector to desired size and set all elements to zero    
  elemVec.Resize(nrEdges); 
  for (Integer i=0; i<nrEdges; i++)
    elemVec[i] = 0;
  
  Vector<Double> currentVec(dim, 0);
  // is needed, if coil is curved
  Vector<Double> globCoord;
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
	  
	  len = currentVec.NormL2();
	  currentVec *= val_ / len;
	  break;

	case 5: // current in yz-plane ==> x=0
	  ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	  currentVec[xDir] = 0;
	  currentVec[yDir] = -(globCoord[zDir] - (*coilMidPt_)[zDir]);
	  currentVec[zDir] = globCoord[yDir] - (*coilMidPt_)[yDir];

	  len = currentVec.NormL2();
	  currentVec *= val_ / len;
	  break;

	case 6: // current in xz-plane ==> y=0
	  ptelem->GetGlobalEdgeIndicesAtIP(globCoord, actIntPt, ptCoord);
	  currentVec[xDir] = globCoord[zDir] - (*coilMidPt_)[zDir];
	  currentVec[yDir] = 0;
	  currentVec[zDir] = -(globCoord[xDir] - (*coilMidPt_)[xDir]);

	  len = currentVec.NormL2();
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

      for(Integer i=0; i<partElemVec.GetSize(); i++)
	elemVec[i] += partElemVec[i] * intWeights[actIntPt-1] * jacDet;	
    }
  

#ifdef DEBUG 
  (*debug) << "CalcElemVector:  "  << std::endl
	   << partElemVec << std::endl
	   << "\n jacDet " << jacDet << std::endl;
#endif


}


// =============================================================================
// volume source integration
// =============================================================================

VolumeSrcInt::VolumeSrcInt(Double aVal, Boolean isaxi)
  : LinearForm(), val_(aVal)
{
  ENTER_FCN( "VolumeSrcIntInt::VolumeSrcInt" );
  isaxi_ = isaxi;
}


VolumeSrcInt::~VolumeSrcInt()
{
  ENTER_FCN( "VolumeSrcIntInt::~VolumeSrcInt" );
}


void VolumeSrcInt::CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & elemVec)
{
  ENTER_FCN( "VolumeSrcInt::CalcElemVector" );

  const Integer nrIntPts = ptelem->GetNumIntPoints();
  const Integer nrNodes  = ptelem->GetNumNodes();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  
  Vector<Double> shapeFnc;
  
  elemVec.Resize(nrNodes);
  elemVec *= 0;    // set elems to 0
  
  Double factor;
  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
    {  
      ptelem->GetShFncAtIp(shapeFnc,actIntPt);
      factor = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord) * intWeights[actIntPt-1] * val_;
      
      if (isaxi_)
	{
	  Vector<Double> CoordAtIP;
	  CoordAtIP = ptCoord * shapeFnc;
	  factor *=  2 * PI * CoordAtIP[0];
	}
      
      shapeFnc *= factor;
      elemVec += shapeFnc;
    }
}


// =============================================================================
// permanent magnet in 2D
// =============================================================================

MagPerm2DInt::MagPerm2DInt(Vector<Double> vecVal, Double rel, Boolean isaxi)
  : LinearForm(), perm_(vecVal), reluctivity_(rel)
{
  ENTER_FCN( "MagPerm2DInt::VolumeSrcInt" );
  isaxi_ = isaxi;
}


MagPerm2DInt::~MagPerm2DInt()
{
  ENTER_FCN( "MagPerm2DInt::~MagPerm2DInt" );
}


void MagPerm2DInt::CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & elemVec)
{
  ENTER_FCN( "MagPerm2DInt::CalcElemVector" );

  const Integer nrIntPts = ptelem->GetNumIntPoints();
  const Integer nrNodes  = ptelem->GetNumNodes();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  
  Vector<Double> ShpFncAtIp, helpVec, CoordAtIP;
  Matrix<Double> xiDx;

  elemVec.Resize(nrNodes);
  elemVec.Init(0);  
  helpVec.Resize(nrNodes);
  helpVec.Init(0);
  
  Double factor;
  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

    Double jacDet = 0;
    ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

    if (isaxi_)
      {
	ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	CoordAtIP = ptCoord * ShpFncAtIp;
	for (Integer i=0; i<nrNodes; i++)
	    xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
	
	jacDet *= 2 * PI * CoordAtIP[0];
      }

    helpVec = xiDx * perm_;
    helpVec *=  intWeights[actIntPt-1] * jacDet * reluctivity_;
    elemVec += helpVec;
  }
}


// ==================================================================
// nLinMagnetics
// ==================================================================

nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt(BaseFE * aptelem, MaterialData & matData, 
												   Boolean isaxi) 
  : LinearForm(aptelem), matData_(matData)
{
  ENTER_FCN( "nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt" );
  isaxi_ = isaxi;
}

nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt(ApproxData *nlinFnc, Double startVal, 
												   Boolean axi)
  : LinearForm()
{
  ENTER_FCN( "nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt" );
  
  isaxi_       = axi;
  startmatVal_ = startVal;
  nlinFnc_     = nlinFnc;
}
  
nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt(Double startVal, Boolean axi)
  : LinearForm()
{
  ENTER_FCN( "nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt" );
  
  isaxi_       = axi;
  startmatVal_ = startVal;
  nlinFnc_     = NULL;   //since this is a linear subdomain
}

nLinMagNode2D_linFormInt::~nLinMagNode2D_linFormInt()
{
  ENTER_FCN( "nLinMagNode2D_linFormInt ::~nLinMagNode2D_linFormInt" );  
}

void nLinMagNode2D_linFormInt::CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & elemVec)
{
  ENTER_FCN("nLinMagNode2D_linFormInt :: ~CalcElemVector" );
  
  const Integer nrNodes  = ptelem->GetNumNodes();
  
  BaseForm * curlcurl2D;
  if (nlinFnc_== NULL) 
	//define the linear element matrix
	curlcurl2D = new CurlCurlNode2DInt(startmatVal_, isaxi_);
  else {
	//define the nonlinear element matrix
	curlcurl2D = new nLinCurlCurlNode2DInt(nlinFnc_, startmatVal_, isaxi_);
	//important to set method to FixPoint, since we compute the RHS!!
	curlcurl2D->SetNonLinMethod("fixPoint");
	
	//set the element solution vector to the bilinearform
	curlcurl2D->SetActElemSol(magPotinMatrix_);
  }

  //set the element-pointer
  curlcurl2D->SetElemPtr(ptelem);
  
  Matrix<Double> elemmat;
  curlcurl2D->CalcElementMatrix(ptCoord, elemmat);
  
  elemVec.Resize(nrNodes);
  elemVec = -elemmat * magPot_;
}




// ==================================================================
// nLinMech
// ==================================================================


nLinMech_linFormInt::nLinMech_linFormInt(BaseFE * aptelem, 
					 MaterialData & matData, 
					 Boolean isaxi) 
  : LinearForm(aptelem), matData_(matData)
{
  ENTER_FCN( "nLinMech_linFormInt::nLinMech_linFormInt" );
  isaxi_ = isaxi;
}


nLinMech_linFormInt::nLinMech_linFormInt(MaterialData & matData, Boolean isaxi) 
  : LinearForm(), matData_(matData)
{
  ENTER_FCN( "nLinMech_linFormInt::nLinMech_linFormInt" );
  isaxi_ = isaxi;
}



nLinMech_linFormInt::~nLinMech_linFormInt()
{
  ENTER_FCN( "nLinMech_linFormInt::~nLinMech_linFormInt" );
}



void nLinMech_linFormInt::CalcElemVector(Matrix<Double>& ptCoord, 
					 Vector<Double> & elemVec)
{
  ENTER_FCN( "nLinMech_linFormInt::CalcElemVector" );

  const Integer nrIntPts = ptelem->GetNumIntPoints();
  const Integer nrNodes  = ptelem->GetNumNodes();
  // getNrDofs() would not work, because CalcElemVec is used for 2d & 3d !
  //    const Integer nrDofs   = getNrDofs();
  const Integer nrDofs   = ptelem->GetDim(); 
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  
  Vector<Double> piolaStressVec;    
  Vector<Double> partElemVec;
  
  Matrix<Double> linBMat; 
  Matrix<Double> nonLinBMat; 
  Matrix<Double> dMat;
  Matrix<Double> transpSumB;    // we need transposed of the b-matrices
  
  
  nLinMechInt_PiolaStress * stressBiformInt;
  
  // These, as friend defined bilinearforms holds the necessary differential operators
  if (ptelem->GetDim() == 2 && !isaxi_)
    stressBiformInt = new nLinMechPlaneStrainInt_PiolaStress(ptelem, matData_);
  else if (ptelem->GetDim() == 2 && isaxi_)
    stressBiformInt = new nLinMechAxiInt_PiolaStress(ptelem, matData_);
  else if (ptelem->GetDim() == 3)
    stressBiformInt = new nLinMech3dInt_PiolaStress (ptelem, matData_);
  else 
    Error("Wrong space dimension of elements! ",__FILE__,__LINE__);
  
  
  
  if (!elemDisp_.GetSizeRow() || !elemDisp_.GetSizeCol()) 
    Error("Undefined displacements! ",__FILE__,__LINE__);
  
  partElemVec.Resize(nrNodes * nrDofs);
  
  
  elemVec.Resize(nrNodes*nrDofs);
  elemVec *= 0;    // set elems to 0
  
  
  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
    {    
      stressBiformInt->setActElemDispl(elemDisp_);
      stressBiformInt->CalcStressVec(piolaStressVec, actIntPt, ptCoord);
      stressBiformInt->calcNonLinBMat(nonLinBMat, actIntPt, ptCoord);
      stressBiformInt->calcLinBMat(linBMat, actIntPt, ptCoord);
      
      nonLinBMat += linBMat;
      
      nonLinBMat.Transpose(transpSumB);
      
      partElemVec = transpSumB * piolaStressVec;
      
      Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
      
      if (isaxi_)
	{
	  Vector<Double> shpFncAtIp;
	  Vector<Double> coordAtIP;
	  ptelem->GetShFncAtIp(shpFncAtIp, actIntPt);
	  coordAtIP = ptCoord * shpFncAtIp;
	  jacDet *= 2 * PI * coordAtIP[0];
	}
      
      partElemVec *= jacDet * intWeights[actIntPt-1];
      
      // the negative sign is due the fact, that this vector has to be subtracted from the RHS!! 
      // (see Kaltenbacher  "Numerical Sim. of Mechatronic Sensors and Actuators" p. 55
      elemVec -=  partElemVec;
    }
  delete stressBiformInt;
}




// ==================================================================
//  recovery technique. calculation of element matrix for RHS
// ==================================================================

RHSForRecoveryProcedure::RHSForRecoveryProcedure(BaseFE * aptelem) 
  : LinearForm(aptelem)
{
  ENTER_FCN( "RHSForRecoveryProcedure::RHSForRecoveryProcedure" );
}
 
RHSForRecoveryProcedure::~RHSForRecoveryProcedure()
{
  ENTER_FCN( "RHSForRecoveryProcedure::~RHSForRecoveryProcedure" );
}

void  RHSForRecoveryProcedure::CalcElemVectorRHSForSPR(Matrix<Double>& ptCoord,
						       Vector<Double> & fncNodesElem,
						       const Integer aComponent,
						       Vector<Double> & elemVec)
{
  ENTER_FCN( "RHSForRecoveryProcedure::CalcElemVectorRHSForSPR" );
	     
  const Integer nrIntPnts = ptelem->GetNumIntPoints();
  const Integer nrNodes   = ptelem->GetNumNodes();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  
  Vector<double> shFnc;
  Matrix<Double>      dvShFnc ;
  Double              jacDet;  
  Double              valFnc=0;
  Integer             iIntPnts,iShFnc;
  
  elemVec.Resize(nrNodes);
  
  for (iIntPnts=0; iIntPnts < nrIntPnts; iIntPnts ++)
    {
      valFnc = 0;
      
      jacDet = ptelem->CalcJacobianDetAtIp(iIntPnts+1, ptCoord);
      
      ptelem->GetGlobDerivShFncAtIp(dvShFnc, iIntPnts+1, ptCoord, jacDet);
      
      ptelem->GetShFncAtIp(shFnc, iIntPnts+1);
      
      for (iShFnc = 0; iShFnc < nrNodes; iShFnc ++)
	valFnc += fncNodesElem[iShFnc]*dvShFnc[iShFnc][aComponent];
      
      for (iShFnc = 0; iShFnc < nrNodes; iShFnc++) 
	elemVec[iShFnc]+=jacDet*valFnc*shFnc[iShFnc]*intWeights[iIntPnts];
      
    } // loop over integration points
  
}

// ==================================================================
// prestress linearform
// ==================================================================

PreStressLinFormInt::PreStressLinFormInt(BaseFE * aptelem, 
					 MaterialData & mat, 
					 Double aPreStressVal, 
					 Directions stressDir) 
  :nLinMech_linFormInt(aptelem, mat), 
   preStressVal_(aPreStressVal), preStressDir_(stressDir)
  
{
  ENTER_FCN( "PreStressLinFormInt::PreStressLinFormInt" );
}



PreStressLinFormInt::~PreStressLinFormInt()
{
  ENTER_FCN( "PreStressLinFormInt::~PreStressLinFormInt" );
}




void PreStressLinFormInt::CalcElemVector(Matrix<Double>& ptCoord, 
					 Vector<Double> & elemVec)
{
  ENTER_FCN( "PreStressLinFormInt::CalcElemVector" );
  Error("PreStressLinFormInt::CalcElemVector: not working");

//   const Integer nrIntPts = ptelem->GetNumIntPoints();
//   const Integer nrNodes  = ptelem->GetNumNodes();
//   const Integer nrDofs   = getNrDofs();
//   const Vector<Double> & intWeights = ptelem->GetIntWeights();  
//   Vector<Double> piolaStressVec;    
//   Vector<Double> partElemVec;
  
//   Matrix<Double> linBMat; 
//   Matrix<Double> nonLinBMat; 
//   Matrix<Double> dMat;
//   Matrix<Double> transpSumB;    // we need transposed of the b-matrices
  

//   // This, as friend defined bilinearform holds the necessary differential operators
//   PreStressInt preStressBiformInt(ptelem, matData_, preStressVal_, preStressDir_);

//   if (!elemDisp_.GetSizeRow() || !elemDisp_.GetSizeCol()) 
//     Error("Undefined displacements! ",__FILE__,__LINE__);

//   partElemVec.Resize(nrNodes * nrDofs);


//   elemVec.Resize( nrNodes * nrDofs );
//   elemVec *= 0;    // set elems to 0
   


//   for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
//     {    
//       preStressBiformInt.setActElemDispl(elemDisp_);
//       preStressBiformInt.CalcStressVec(piolaStressVec, actIntPt, ptCoord);
//       preStressBiformInt.calcNonLinBMat(nonLinBMat, actIntPt, ptCoord);
//       preStressBiformInt.calcLinBMat(linBMat, actIntPt, ptCoord);

//       nonLinBMat += linBMat;

//       nonLinBMat.Transpose(transpSumB);
	
//       partElemVec = transpSumB * piolaStressVec;

//       Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
//       partElemVec *= jacDet * intWeights[actIntPt-1];
	
//       elemVec +=  partElemVec;
//     }
}


// ========================================================================
// pressureLinForm 
// ========================================================================


PressureLinForm::PressureLinForm(Double aVal, Boolean isaxi)
  : LinearForm(), multiplier_(aVal)
{
  ENTER_FCN( "PressureLinForm::PressureLinForm" );
  isaxi_ = isaxi;
}



PressureLinForm::~PressureLinForm()
{
  ENTER_FCN( "PressureLinForm::~PressureLinForm" );
}


void PressureLinForm::CalcElemVector(Matrix<Double>& ptCoord, 
				       Vector<Double> & elemVec)
{
  ENTER_FCN( "PressureLinForm::CalcElemVector" );

  const Integer nrIntPts = ptelem->GetNumIntPoints();
  const Integer nrNodes  = ptelem->GetNumNodes();
  const Integer dim      = ptCoord. GetSizeRow();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();
  Vector<Double> shapeFnc;
  Vector<Double> normalVec(dim);

  //comput normal vector
  if (dim == 2) 
    {
      Double dx  = ptCoord[0][1] - ptCoord[0][0];
      Double dy  = ptCoord[1][1] - ptCoord[1][0];
      Double len = sqrt(dx*dx + dy*dy);
      if (len <= 0.0)
	Error("length of normal vector is zero!");
      normalVec[0] = dy/len;
      normalVec[1] = -dx/len;
    }
  else 
    {
      //compute the two vectors in the plane
      Vector<Double> vec1(3), vec2(3);
      for (Integer i=0; i<dim; i++) {
	vec1[i] = ptCoord[i][1] - ptCoord[i][0];
	vec2[i] = ptCoord[i][nrNodes-1] - ptCoord[i][0];
      }
      //compute cross product
      normalVec[0] = vec1[1] * vec2[2] - vec1[2]*vec2[1];
      normalVec[1] = vec1[2] * vec2[0] - vec1[0]*vec2[2];
      normalVec[2] = vec1[0] * vec2[1] - vec1[1]*vec2[0];
      //normalize the length to 1
      Double length = normalVec.NormL2();
      normalVec /= length;
    }

  Integer zerodim = 0;
  if (dofzero_ >0) {
    //now ndDofs = dim + 1 for electric potential
    elemVec.Resize(nrNodes*(dim+1));
    zerodim = 1;
  }
  else
    elemVec.Resize(nrNodes*dim);

  elemVec.Init(0);    // set elems to 0

  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
    {
      ptelem->GetShFncAtIp(shapeFnc, actIntPt);
      Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
      Double factor = multiplier_ * intWeights[actIntPt-1] * jacDet;

      if (isaxi_)
	{
	  Vector<Double> CoordAtIP;
	  CoordAtIP = ptCoord * shapeFnc;
	  factor *=  2 * PI * CoordAtIP[0];
	}
      
      shapeFnc *= factor;      
     
      Vector<Double> helpVec;
      Integer idx;
      for (Integer i=0; i<dim; i++) 
	{
	  //multiply with corresponding component of normal vector
	  //to get the x-,y-,z-component
	  helpVec = shapeFnc * normalVec[i];
	  for (Integer j=0; j<helpVec.GetSize(); j++) {
	    idx = j*(dim+zerodim) +i;
	    elemVec[idx] -= helpVec[j];
	  }
	}
    }
} // end of method


//==================================================================
//Members of class LinearFormFlowNoise
//==================================================================

LinearFlowNoiseInt::LinearFlowNoiseInt(BaseFE * aptelem) : LinearForm(aptelem)
{
  ENTER_FCN( "LinearFlowNoiseInt::LinearFlowNoiseInt" );
}

LinearFlowNoiseInt::~LinearFlowNoiseInt()
{
  ENTER_FCN( "LinearFlowNoiseInt::~LinearFlowNoiseInt" );
}

void LinearFlowNoiseInt::CalcElemVector4Dip(Matrix<Double>& ptCoord,
					    const StdVector<Integer> & connecth, 
					    Vector<Double> & Result, 
					    const Vector<Double> gradN_x_P)
{
  ENTER_FCN( "LinearForm::CalcElemVector4FlowSrc" );

  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();
  
  Integer i,ii;
  
  Double jacDet;    
  //     Double density=1.0;
  Result.Resize(n);
  for (i=0;i<n;i++)
    Result[i]=0.0;
  
  Vector<Double>  Sf;
  Double multiplier;
  // Dipole Source Term: (w,deltp)onbnd
  Vector<Double> normal;
  
  normal.Resize(gradN_x_P.GetSize());

  normal*=0;
  
  switch(normal.GetSize())
    {
    case 2:
      {
	calcNormal2Line_Mat(normal,ptCoord);
	break;
      }
    case 3:
      {
	calcNormal2Surface_Mat(normal,ptCoord);
	break;
      }
    default:
      {
	std::string errMsg = "Incorrect dimension of normal. ";
	Error(errMsg.c_str(),__FILE__,__LINE__);
      }
    }
  
  for (i=0; i<l; i++)
    {
      jacDet = ptelem->CalcJacobianDetAtIp(i+1, ptCoord);
      
      ptelem->GetShFncAtIp(Sf, i+1);
      
      for (ii=0; ii < n; ii++)
	{
	  multiplier= normal * gradN_x_P; // n * gradShFnc
	  std::cout<<"multiplier: "<<multiplier<<std::endl;
	  std::cout<<"jacDet: "<<jacDet<<std::endl;	
	  std::cout<<"Sf["<<i<<"]: "<<Sf[i]<<std::endl;    
	  Result[ii]+= -jacDet*multiplier*Sf[ii];
	}
      std::cout<<"Result: "<<Result<<std::endl; 
    }
  
} // end of method


void LinearFlowNoiseInt::CalcElemVector4Quad(Matrix<Double>& ptCoord,
					     const StdVector<Integer> & connecth,
					     const Matrix<Double> & FlowData, 
					     Vector<Double> & Result)
{
  ENTER_FCN( "LinearFlowNoiseInt::CalcElemVector4Quad" );

  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();
  Matrix<Double> xiDx;      
  Matrix<Double> elVel;
  Matrix<Double> elVelAtIP;    
  Double jacDet;
  Integer actInt;
  Double density=1.0;
  Result.Resize(n);

  for (int i=0;i<n;i++)
    Result[i]=0.0;

  Vector<Double>  Sf;
  Vector<Double> partResult;
  
  int dimelem=ptCoord.GetSizeRow();

//     std::vector<Double>  VelAtIP;
//     Matrix<Double> VelDerAtIP;
//     Matrix<Double> VelDerFromDiag;
//     Matrix<Double> VelMatrixforMult;
//     VelMatrixforMult.Resize(dimelem); 

//     std::vector<Double> dTij_di;
//     VelAtIP.resize(dimelem);
//     VelDerAtIP.Resize(dimelem);  
//     dTij_di.resize(dimelem);
//     std::vector<Double> helpVect;

  
  Vector<Double>  VelAtIP;
  Matrix<Double> VelDerAtIP;
  Matrix<Double> VelDerFromDiag;
  Matrix<Double> VelMatrixforMult;
  VelMatrixforMult.Resize(dimelem); 

  Vector<Double> dTij_di;
  VelAtIP.Resize(dimelem);
  VelDerAtIP.Resize(dimelem);  
  dTij_di.Resize(dimelem);
  Vector<Double> helpVect;

//     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// //     // For using with Anne Le Duc source
//     std::vector<Double> ddTij_dxidxi;
//     ddTij_dxidxi.resize(n);
// 	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  GetQttiesOfElement(elVel, FlowData, connecth, dimelem);
  

// 	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 	// Just for testing with source term from Anne Le Duc (ddTij/dxidxj) in MixLayer example
//     for (actInt=1; actInt<=l; actInt++)
//       {
// 	ddTij_dxidxi[actInt-1]=0; // Set to zero index to be filled in.
// 	ptelem->GetShFncAtIp(Sf, actInt);
// 	for (int ctrIP=1; ctrIP<=l; ctrIP++)
//  	  {
// 	    // Interpolate to IP and fill in vector
//  	ddTij_dxidxi[actInt-1]+=(FlowData[0][connecth[ctrIP-1]-1])*Sf[actInt-1]; // In files from MixL first value is ddTij/dxidxj!!
//  	  }
//       }
// 	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    elVelAtIP.Resize(dimelem, l);
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt);
	//jacDet = ptelem->CalcJacobianDetAtIp(actInt, ptCoord);
	VelAtIP=elVel*Sf;
	
	
	for (int comp=0;comp<dimelem;comp++)
	  {
	    elVelAtIP[comp][actInt-1]=VelAtIP[comp]; //Filling the matrix of velocity components at IP
	  }
      }
    

    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet);
      
	// This work only for quad1 and trilinear hexahedrals elements since we have values
	// of the flow quantities only at the corners!
	// Here we compute the derivatives of the Lighthill's tensor needed in the quadrupole term
	// at the ith integration point
     
	//Implementation 26.09.03
	//Modified 15.06.04



	//Wrong mathematic derivation!! 
// 	VelDerAtIP.GetDiagInMatrix(VelDerFromDiag);	 
// 	VelDerFromDiag.ConvertToVec_AppendRows(helpVect);
// 	for (int k=0;k<(dimelem-1);k++)
// 	  VelDerFromDiag.add_col(helpVect,1);

// 	VelDerFromDiag.ScaleDiagElems(2.0);
//	dTij_di=(VelDerFromDiag*VelAtIP);

// Implementation of new derivation

	VelAtIP=elVel*Sf;
	

	VelDerAtIP=(elVelAtIP*xiDx); // Now this is computed more accurately using elVelAtIP, before was elVel
	//VelDerAtIP*=jacDet; 

	VelDerAtIP.GetDiagInMatrix(VelDerFromDiag);
	VelDerFromDiag.ConvertToVec_AppendRows(helpVect);
	dTij_di=(VelDerAtIP*VelAtIP);
	
	for (int k=0;k<(dimelem);k++)
	  for (int j=0;j<(dimelem);j++)
	    VelMatrixforMult[j][k]=VelAtIP[j];

	//std::vector<Double> tempVect;
	Vector<Double> tempVect;
	tempVect.Resize(dimelem);
	
	tempVect=VelMatrixforMult*helpVect;

	dTij_di+=tempVect;

// END Implementation of new derivation



	dTij_di*=density;

	partResult=xiDx*dTij_di;
	partResult*=jacDet;




// 	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 	// Just for testing with source term from Anne Le Duc (ddTij/dxidxj) in MixLayer example
// 	ddTij_dxidxi*=density;
// 	for (int ii=0;ii<n;ii++)
// 	  {
// 	    partResult[ii]=Sf[ii]*ddTij_dxidxi[ii];
// 	    partResult[ii]*=jacDet;
// 	  }
// 	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	Result+=partResult;
      }

  } // end of method
 

void LinearFlowNoiseInt::GetQttiesOfElement(Matrix<Double>& elVel,
					    const Matrix<Double>& FlowData,
					    const StdVector<Integer>& connecth, 
					    Integer matrixRow)
{
  ENTER_FCN( "LinearFlowNoiseInt::GetVecOfElement" );
 
  // displacement of element nodes
  elVel.Resize(matrixRow, connecth.GetSize());
   
  for(Integer actNode=0; actNode<connecth.GetSize(); actNode++)
    for (int dim=0; dim < matrixRow; dim++)
      elVel[dim][actNode] = FlowData[dim+1][connecth[actNode]-1]; // dim+1 because index 0 in FlowData is used for storing pressure
}


// =============================================================================
// nonlinear RHS for nonlinear acoustics
// =============================================================================

nLinKuznetsovRHSInt::nLinKuznetsovRHSInt(Double aVal, Boolean isaxi)
  : LinearForm(), val_(aVal)
{
  ENTER_FCN( "nLinKuznetsovRHSInt::nLinKuznetsovRHSInt" );
  isaxi_ = isaxi;
}


nLinKuznetsovRHSInt::~nLinKuznetsovRHSInt()
{
  ENTER_FCN( "nLinKuznetsovRHSInt::~nLinKuznetsovRHSInt" );
}


void nLinKuznetsovRHSInt::CalcElemVector(Matrix<Double>& ptCoord, Vector<Double> & elemVec)
{
  ENTER_FCN( "nLinKuznetsovRHSInt::CalcElemVector" );

  const Integer nrIntPts = ptelem->GetNumIntPoints();
  const Integer nrNodes  = ptelem->GetNumNodes();
  const Vector<Double> & intWeights = ptelem->GetIntWeights();  

  // derivation of shape functions after global coordinates 
  Matrix<Double> xiDx;
  Matrix<Double> xiDxTransp;

  Vector<Double> ShpFncAtIp;
  Vector<Double> CoordAtIp;
  Vector<Double> solGradAtIp;
  Vector<Double> solDeriv1GradAtIp;
  Double solDeriv1AtIp;
  Double solDeriv2AtIp;

  Double jacDet;
  
  elemVec.Resize(nrNodes);
  elemVec.Init(0.0);

  Double factor;
  for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

	jacDet = 0;
	ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);
	
	if (isaxi_) {
	  CoordAtIp = ptCoord * ShpFncAtIp;
	  for (Integer i=0; i<nrNodes; i++)
		xiDx[i][0] += ShpFncAtIp[i] / CoordAtIp[0];
	  
	  jacDet *= 2 * PI * CoordAtIp[0];
	}
	
	xiDx.Transpose(xiDxTransp);
	
	//compute gradient of solution and 1st derivative at integration point
	solGradAtIp       = xiDxTransp * sol_;
	solDeriv1GradAtIp = xiDxTransp * solderiv1_;
	
	//get 1st and 2nd derivartive of solution at integration point
	solDeriv1AtIp = solderiv1_*ShpFncAtIp;
	solDeriv2AtIp = solderiv2_*ShpFncAtIp;
	
	Double factor=0;
	for (Integer j=0; j<xiDx.GetSizeCol(); j++)
	  factor += solGradAtIp[j]*solDeriv1GradAtIp[j];
	factor *= factorN2_;
	
	factor += factorN1_ * solDeriv1AtIp * solDeriv2AtIp;
	
	factor *= jacDet;
	for (Integer i=0; i< nrNodes; i++)
	  elemVec[i] += ShpFncAtIp[i] * factor;
	
  }
  
  //  std::cerr << "RHS in linearForm:\n" << elemVec << std::endl;
}


} // end of namespace
