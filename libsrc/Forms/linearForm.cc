#include <iostream>
#include <fstream>
#include <math.h>

#include "forms_header.hh"

namespace CoupledField
{

  LinearForm::LinearForm(BaseFE * aptelem) : BaseForm(aptelem)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::LinearForm" << std::endl;
#endif
  }


  LinearForm::LinearForm() : BaseForm()
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






// =============================================================================
// edge integration
// =============================================================================


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




// =============================================================================
// volume source integration
// =============================================================================



  VolumeSrcInt::VolumeSrcInt(Double aVal, Boolean isaxi)
    : LinearForm(), val_(aVal)
  {
#ifdef TRACE
    (*trace) << "entering VolumeSrcIntInt::VolumeSrcInt" << std::endl;
#endif
    isaxi_ = isaxi;
  }


  VolumeSrcInt ::~VolumeSrcInt()
  {
#ifdef TRACE
    (*trace) << "entering VolumeSrcIntInt::~VolumeSrcInLinearEdgeInt" << std::endl;
#endif
  }


  void VolumeSrcInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering VolumeSrcInt::CalcElemVector" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrNodes  = ptelem->GetNumNodes();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    std::vector<Double> shapeFnc;

    elemVec.resize(nrNodes);
    elemVec *= 0;    // set elems to 0
   
    Double factor;
    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {  
	ptelem->GetShFncAtIp(shapeFnc,actIntPt);
	factor = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord) * intWeights[actIntPt-1] * val_;

	if (isaxi_)
	  {
	    std::vector<Double> CoordAtIP;
	    CoordAtIP = ptCoord * shapeFnc;
	    factor *=  2 * PI * CoordAtIP[0];
	  }

	shapeFnc *= factor;
	elemVec += shapeFnc;
      }

#ifdef TRACE
    (*trace) << "leaving VolumeSrcInt::CalcElementVector" << std::endl;
#endif

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

  nLinMagNode2D_linFormInt ::~ nLinMagNode2D_linFormInt()
  {
    ENTER_FCN( "nLinMagNode2D_linFormInt ::~nLinMagNode2D_linFormInt" );

  }

  void nLinMagNode2D_linFormInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
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

      //set the element solution vector to the bilinearform
      curlcurl2D->SetActElemSol(magPotinMatrix_);
    }

    //set the element-pointer
    curlcurl2D->SetElemPtr(ptelem);

    Matrix<Double> elemmat;
    curlcurl2D->CalcElementMatrix(ptCoord, elemmat);
    
    elemVec.resize(nrNodes);
    elemVec = -elemmat * magPot_;
  }




  // ==================================================================
  // nLinMech
  // ==================================================================


  nLinMech_linFormInt::nLinMech_linFormInt(BaseFE * aptelem, MaterialData & matData, Boolean isaxi) 
    : LinearForm(aptelem), matData_(matData)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::nLinMech_linFormInt" << std::endl;
#endif
    isaxi_ = isaxi;
  }


  nLinMech_linFormInt::nLinMech_linFormInt(MaterialData & matData, Boolean isaxi) 
    : LinearForm(), matData_(matData)
  {
#ifdef TRACE
    (*trace) << "entering nLinMech_linFormInt::nLinMech_linFormInt" << std::endl;
#endif
    isaxi_ = isaxi;
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
    //    const Integer nrDofs   = getNrDofs();
    const Integer nrDofs   = ptelem->GetDim(); // getNrDofs() would not work, because CalcElemVec is used for 2d & 3d !
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    std::vector<Double> piolaStressVec;    
    std::vector<Double> partElemVec;

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

    partElemVec.resize(nrNodes * nrDofs);


   elemVec.resize(nrNodes*nrDofs);
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
	    std::vector<Double> shpFncAtIp;
	    std::vector<Double> coordAtIP;
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
#ifdef TRACE
    (*trace) << "entering RHSForRecoveryProcedure::RHSForRecoveryProcedure" << std::endl;
#endif
  }
 
  RHSForRecoveryProcedure::~RHSForRecoveryProcedure()
  {
#ifdef TRACE
    (*trace) << "entering RHSForRecoveryProcedure::~RHSForRecoveryProcedure" << std::endl;
#endif
  }

void  RHSForRecoveryProcedure::CalcElemVectorRHSForSPR(Matrix<Double>& ptCoord,
					       Vector<Double> & fncNodesElem,
						       const Integer aComponent,
						       Vector<Double> & elemVec)
{
#ifdef TRACE
    (*trace) << "entering RHSForRecoveryProcedure::CalcElemVectorRHSForSPR()" << std::endl;
#endif

    const Integer nrIntPnts = ptelem->GetNumIntPoints();
    const Integer nrNodes   = ptelem->GetNumNodes();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
    std::vector<double> shFnc;
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

  PreStressLinFormInt::PreStressLinFormInt(BaseFE * aptelem, MaterialData & mat, Double aPreStressVal, Directions stressDir) 
    :nLinMech_linFormInt(aptelem, mat), preStressVal_(aPreStressVal), preStressDir_(stressDir)
    
  {
#ifdef TRACE
    (*trace) << "entering PreStressLinFormInt::PreStressLinFormInt" << std::endl;
#endif
  }


  
  PreStressLinFormInt ::~PreStressLinFormInt()
  {
#ifdef TRACE
    (*trace) << "entering PreStressLinFormInt::~PreStressLinFormInt" << std::endl;
#endif
  }




  void PreStressLinFormInt::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering PreStressLinFormInt::CalcElemVector" << std::endl;
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
    PreStressInt preStressBiformInt(ptelem, matData_, preStressVal_, preStressDir_);

    if (!elemDisp_.GetSizeRow() || !elemDisp_.GetSizeCol()) 
      Error("Undefined displacements! ",__FILE__,__LINE__);

    partElemVec.resize(nrNodes * nrDofs);


   elemVec.resize( nrNodes * nrDofs );
   elemVec *= 0;    // set elems to 0
   


    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
	preStressBiformInt.setActElemDispl(elemDisp_);
	preStressBiformInt.CalcStressVec(piolaStressVec, actIntPt, ptCoord);
	preStressBiformInt.calcNonLinBMat(nonLinBMat, actIntPt, ptCoord);
	preStressBiformInt.calcLinBMat(linBMat, actIntPt, ptCoord);

	nonLinBMat += linBMat;

	nonLinBMat.Transpose(transpSumB);
	
	partElemVec = transpSumB * piolaStressVec;

	Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
	partElemVec *= jacDet * intWeights[actIntPt-1];
	
	elemVec +=  partElemVec;
      }
  }


  // ========================================================================
  // SurfaceIntLinForm 
  // ========================================================================


  SurfaceIntLinForm::SurfaceIntLinForm(BaseFE * aptelem, Directions aPressureDir) 
    : LinearForm(aptelem), pressureDir_(aPressureDir)
  {
#ifdef TRACE
    (*trace) << "entering SurfaceIntLinForm::SurfaceIntLinForm" << std::endl;
#endif
  }



  SurfaceIntLinForm ::~SurfaceIntLinForm()
  {
#ifdef TRACE
    (*trace) << "entering SurfaceIntLinForm::~SurfaceIntLinForm" << std::endl;
#endif
  }


void SurfaceIntLinForm::CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & elemVec)
  {
#ifdef TRACE
    (*trace) << "entering SurfaceIntLinForm::CalcElemVector" << std::endl;
#endif

    const Integer nrIntPts = ptelem->GetNumIntPoints();
    const Integer nrNodes  = ptelem->GetNumNodes();
    const Integer nrDofs   = getNrDofs();
    const std::vector<Double> & intWeights = ptelem->GetIntWeights();
    std::vector<Double> shapeFnc;
    Integer pressureDof;
    

    std::vector<Double> partElemVec;
    partElemVec.resize(nrNodes * nrDofs );
    partElemVec *= 0; //initialize vec
    

    elemVec.resize(nrNodes*nrDofs);
    elemVec *= 0;    // set elems to 0

    switch(pressureDir_)
      {
      case X: pressureDof = 0; break;
      case Y: pressureDof = 1; break;
      case Z: pressureDof = 2; break;
      default: 	Error("SurfaceIntLinForm: pressure direction not implemented!" 
		      ,__FILE__,__LINE__);
      }
    
    for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {
	ptelem->GetShFncAtIp(shapeFnc, actIntPt);

	Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
	
	std::vector<Double> helpVec = multiplier_ * intWeights[actIntPt-1] * jacDet * shapeFnc;

	SetSubVector(partElemVec, helpVec, pressureDof*nrNodes);
	elemVec += partElemVec;
      }
  } // end of method





  

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

void LinearFlowNoiseInt::CalcElemVector4Dip(Matrix<Double>& ptCoord,const Vector<Integer> & connecth, std::vector<Double> & Result, const std::vector<Double> gradN_x_P)
  {
#ifdef TRACE
    (*trace) << "entering LinearForm::CalcElemVector4FlowSrc" << std::endl;
#endif

    Integer l=ptelem->GetNumIntPoints();
    Integer n=ptelem->GetNumNodes();

    Integer i,ii;

    Double jacDet;    
//     Double density=1.0;
    Result.resize(n);
    for (i=0;i<n;i++)
      Result[i]=0.0;
    
    std::vector<Double>  Sf;
    Double multiplier;
    // Dipole Source Term: (w,deltp)onbnd
    std::vector<Double> normal;

    normal.resize(gradN_x_P.size());

    normal*=0;
    
    switch(normal.size())
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
	    multiplier=ScalarMult(normal,gradN_x_P); // n * gradShFnc
	    std::cout<<"multiplier: "<<multiplier<<std::endl;
	    std::cout<<"jacDet: "<<jacDet<<std::endl;	
	    std::cout<<"Sf["<<i<<"]: "<<Sf[i]<<std::endl;    
	    Result[ii]+= -jacDet*multiplier*Sf[ii];
	  }
      std::cout<<"Result: "<<Result<<std::endl; 
      }
    
  } // end of method


void LinearFlowNoiseInt::CalcElemVector4Quad(Matrix<Double>& ptCoord,const Vector<Integer> & connecth,const Matrix<Double> & FlowData, std::vector<Double> & Result)
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4Quad" << std::endl;
#endif

    Integer l=ptelem->GetNumIntPoints();
    Integer n=ptelem->GetNumNodes();
    Matrix<Double> xiDx;      
    Matrix<Double> elVel;    
    Double jacDet;
    Integer actInt;
    Double density=1.0;
    Result.resize(n);
    for (int i=0;i<n;i++)
      Result[i]=0.0;
    std::vector<Double>  Sf;
    std::vector<Double> partResult;
  
    int dimelem=ptCoord.GetSizeRow();
  
    std::vector<Double>  NodalVelAtIP;
    Matrix<Double> VelDerAtIP;
    Matrix<Double> VelDerFromDiag;
    std::vector<Double> dTij_di;
    NodalVelAtIP.resize(dimelem);
    VelDerAtIP.Resize(dimelem);  
    dTij_di.resize(dimelem);
    std::vector<Double> helpVect;

    GetQttiesOfElement(elVel, FlowData, connecth, dimelem);
  
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet);
      
	// This work only for quad1 and trilinear hexahedrals elements since we have values
	// of the flow quantities only at the corners!
	// Here we compute the derivatives of the Lighthill's tensor needed in the quadrupole term
	// at the ith integration point
     
	//Implementation 26.09.03
	  
	NodalVelAtIP=elVel*Sf;

	VelDerAtIP=(elVel*xiDx);
	VelDerAtIP*=jacDet; 
	VelDerAtIP.GetDiagInMatrix(VelDerFromDiag);	 
	VelDerFromDiag.ConvertToVec_AppendRows(helpVect);
	for (int k=0;k<(dimelem-1);k++)
	  VelDerFromDiag.AddColumn(helpVect,1);

	VelDerFromDiag.ScaleDiagElems(2.0);

	dTij_di=(VelDerFromDiag*NodalVelAtIP);
	dTij_di*=density;

	partResult=xiDx*dTij_di;
	partResult*=jacDet;

	Result+=partResult;
      }

  } // end of method


void LinearFlowNoiseInt::GetQttiesOfElement(Matrix<Double>& elVel,const Matrix<Double>& FlowData,const Vector<Integer>& connecth, Integer matrixRow)
{
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::GetVecOfElement" << std::endl;
#endif

    // displacement of element nodes
    elVel.Resize(matrixRow, connecth.GetSize());

    for(Integer actNode=0; actNode<connecth.GetSize(); actNode++)
      for (int dim=0; dim < matrixRow; dim++)
	elVel[dim][actNode] = FlowData[dim+1][connecth[actNode]-1]; // dim+1 because index 0 in FlowData is used for storing pressure
}


} // end of namespace
