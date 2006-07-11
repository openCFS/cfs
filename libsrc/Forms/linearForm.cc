#include <iostream>
#include <fstream>
#include <math.h>

#include "forms_header.hh"
#include "Utils/coordSystem.hh"
#include "Domain/domain.hh"

namespace CoupledField {


  LinearForm::LinearForm() : BaseForm( NULL )
  {
    ENTER_FCN( "LinearForm::LinearForm" );
  }




  LinearForm::~LinearForm()
  {
    ENTER_FCN( "LinearForm::~LinearForm" );
  }


  void LinearForm::CalcElemVector( Vector<Double> & result,
                                   EntityIterator& ent ) {
    ENTER_FCN( "LinearForm::CalcElemVector" );
  }


  // ================================================================
  // edge integration
  // ================================================================

  LinearEdgeInt::LinearEdgeInt(Double aVal, 
                               UInt aDirection,
                               Vector<Double> * aCoilMidPt) 
    : val_(aVal), direction_(aDirection)
  {
    ENTER_FCN( "LinearEdgeInt::LinearEdgeInt" );
    name_ = "LinearEdgeInt";
  
    if(aCoilMidPt)
      coilMidPt_ = new Vector<Double>(*aCoilMidPt);
  }

  LinearEdgeInt::~LinearEdgeInt()
  {
    ENTER_FCN( "LinearEdgeInt::~LinearEdgeInt" );
  }



  void LinearEdgeInt::CalcElemVector( Vector<Double> & elemVec,
                                      EntityIterator& ent ) 
  {
    ENTER_FCN( "LinearEdgeInt::CalcElemVector" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrEdges  = ptelem->GetNumEdges();
    const UInt dim      = ptelem->GetDim();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const UInt xDir = 0;
    const UInt yDir = 1;
    const UInt zDir = 2;
  
    Double jacDet;  
  
  
    // derivation of shape functions after global coordinates 
    Matrix<Double> shapeEdge;
    Vector<Double> partElemVec;
    partElemVec.Resize(nrEdges);
    partElemVec.Init();
  
  
    // set vector to desired size and set all elements to zero    
    elemVec.Resize(nrEdges); 
    elemVec.Init();

    Vector<Double> currentVec(dim, 0);
    // is needed, if coil is curved
    Vector<Double> globCoord;
    Double len;
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
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
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord_);
            currentVec[xDir] = -(globCoord[yDir] - (*coilMidPt_)[yDir]);
            currentVec[yDir] = globCoord[xDir] - (*coilMidPt_)[xDir];
            currentVec[zDir] = 0;
          
            len = currentVec.NormL2();
            currentVec *= val_ / len;
            break;

          case 5: // current in yz-plane ==> x=0
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord_);
            currentVec[xDir] = 0;
            currentVec[yDir] = -(globCoord[zDir] - (*coilMidPt_)[zDir]);
            currentVec[zDir] = globCoord[yDir] - (*coilMidPt_)[yDir];

            len = currentVec.NormL2();
            currentVec *= val_ / len;
            break;

          case 6: // current in xz-plane ==> y=0
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord_);
            currentVec[xDir] = globCoord[zDir] - (*coilMidPt_)[zDir];
            currentVec[yDir] = 0;
            currentVec[zDir] = -(globCoord[xDir] - (*coilMidPt_)[xDir]);

            len = currentVec.NormL2();
            currentVec *= val_ / len;

            break;
            
          default:
            std::string errMsg = "Selected current direction with num. ";
            errMsg += direction_ + " not supported!";
            Error(errMsg.c_str(),__FILE__,__LINE__);
          }
        
        ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord_);

        partElemVec = shapeEdge * currentVec;

        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_);

        for(UInt i=0; i<partElemVec.GetSize(); i++)
          elemVec[i] += partElemVec[i] * intWeights[actIntPt-1] * jacDet; 
      }
  

#ifdef DEBUG 
    (*debug) << "CalcElemVector:  "  << std::endl
             << partElemVec << std::endl
             << "\n jacDet " << jacDet << std::endl;
#endif


  }


  // ====================================================================
  // volume source integration
  // ====================================================================

  VolumeSrcInt::VolumeSrcInt(Double aVal, bool isaxi, bool coordUpdate )
    : LinearForm(), 
      val_(aVal)
  {
    ENTER_FCN( "VolumeSrcIntInt::VolumeSrcInt" );
    name_ = "VolumeSrcInt";
    isaxi_ = isaxi;
    coordUpdate_ = coordUpdate;
  }


  VolumeSrcInt::~VolumeSrcInt()
  {
    ENTER_FCN( "VolumeSrcIntInt::~VolumeSrcInt" );
  }


  void VolumeSrcInt::CalcElemVector( Vector<Double> & elemVec,
                                     EntityIterator& ent ) 
  {
    ENTER_FCN( "VolumeSrcInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc;
  
    elemVec.Resize(nrNodes);
    elemVec.Init();
  
    Double factor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {  
        ptelem->GetShFncAtIp(shapeFnc,actIntPt);
        factor = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_) * 
          intWeights[actIntPt-1] * val_;
      
        if (isaxi_)
          {
            Vector<Double> CoordAtIP;
            CoordAtIP = ptCoord_ * shapeFnc;
            factor *=  2 * PI * CoordAtIP[0];
          }
      
        shapeFnc *= factor;
        elemVec += shapeFnc;
      }
  }


  // ===================================================================
  // permanent magnet in 2D
  // ===================================================================

  MagPerm2DInt::MagPerm2DInt(Vector<Double> vecVal, Double rel, 
                             bool isaxi, bool coordUpdate )
    : LinearForm(), 
      perm_(vecVal), 
      reluctivity_(rel)
  {
    ENTER_FCN( "MagPerm2DInt::VolumeSrcInt" );
    name_ = "MagPerm2DInt";
    isaxi_ = isaxi;
    coordUpdate_ = coordUpdate;
  }


  MagPerm2DInt::~MagPerm2DInt()
  {
    ENTER_FCN( "MagPerm2DInt::~MagPerm2DInt" );
  }


  void MagPerm2DInt::CalcElemVector( Vector<Double> & elemVec,
                                     EntityIterator& ent )
  {
    ENTER_FCN( "MagPerm2DInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> ShpFncAtIp, helpVec, CoordAtIP;
    Matrix<Double> xiDx;

    elemVec.Resize(nrNodes);
    elemVec.Init(0);  
    helpVec.Resize(nrNodes);
    helpVec.Init(0);
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = 0;
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);

      if (isaxi_)
        {
          ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
          CoordAtIP = ptCoord_ * ShpFncAtIp;
          for (UInt i=0; i<nrNodes; i++)
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

  nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt(ApproxData *nlinFnc, 
                                                     Double startVal, 
                                                     bool axi,
                                                     bool coordUpdate)
    : LinearForm()
  {
    ENTER_FCN( "nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt" );
    name_ = "nLinMagNode2D_linFormInt";
    isaxi_       = axi;
    coordUpdate_ = coordUpdate;
    startmatVal_ = startVal;
    nlinFnc_     = nlinFnc;
    isSolDependent_ = true;
  }
  
  nLinMagNode2D_linFormInt::~nLinMagNode2D_linFormInt()
  {
    ENTER_FCN( "nLinMagNode2D_linFormInt ::~nLinMagNode2D_linFormInt" );  
  }

  void nLinMagNode2D_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                                 EntityIterator& ent )
  {
    ENTER_FCN("nLinMagNode2D_linFormInt :: ~CalcElemVector" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    const UInt nrNodes  = ptelem->GetNumNodes();
  
    BaseForm * curlcurl2D;
    if (nlinFnc_== NULL) 
      //define the linear element matrix
      curlcurl2D = new CurlCurlNode2DInt(startmatVal_, isaxi_);
    else {
      //define the nonlinear element matrix
      curlcurl2D = new nLinCurlCurlNode2DInt(nlinFnc_, startmatVal_, isaxi_);
      //important to set method to FixPoint, since we compute the RHS!!
      curlcurl2D->SetNonLinMethod("fixPoint");
        
      //set the solution class to the operator
      curlcurl2D->SetSolution( * sol_ );
    }

    Matrix<Double> elemmat;
    curlcurl2D->CalcElementMatrix(elemmat, ent, ent);

    // Get element solution
    Vector<Double> magPot;
    sol_->GetElemSolution( magPot, ent.GetElem()->connect );

    elemVec.Resize(nrNodes);
    elemVec = -(elemmat * magPot);

    delete curlcurl2D;
  }




  // ==================================================================
  // nLinMech
  // ==================================================================



  nLinMech_linFormInt::nLinMech_linFormInt(BaseMaterial* matData,
                                           bool isaxi) 
    : LinearForm(), matData_(matData)
  {
    ENTER_FCN( "nLinMech_linFormInt::nLinMech_linFormInt" );
    name_ = "nLinMech_linFormInt";
    isaxi_ = isaxi;
    isSolDependent_ = true;
  }



  nLinMech_linFormInt::~nLinMech_linFormInt()
  {
    ENTER_FCN( "nLinMech_linFormInt::~nLinMech_linFormInt" );
  }



  void nLinMech_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent )
  {
    ENTER_FCN( "nLinMech_linFormInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    // getNrDofs() would not work, because CalcElemVec is used for 2d&3d!
    //    const UInt nrDofs   = getNrDofs();
    const UInt nrDofs   = ptelem->GetDim(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> piolaStressVec;    
    Vector<Double> partElemVec;
  
    Matrix<Double> linBMat; 
    Matrix<Double> nonLinBMat; 
    Matrix<Double> dMat;
    Matrix<Double> transpSumB;    // we need transposed of the b-matrices
  
  
    nLinMechInt_PiolaStress * stressBiformInt;
  
    // These, as friend defined bilinearforms holds the necessary 
    // differential operators
    if (ptelem->GetDim() == 2 && !isaxi_)
      stressBiformInt = new nLinMechPlaneStrainInt_PiolaStress(matData_);
    else if (ptelem->GetDim() == 2 && isaxi_)
      stressBiformInt = new nLinMechAxiInt_PiolaStress(matData_);
    else if (ptelem->GetDim() == 3)
      stressBiformInt = new nLinMech3dInt_PiolaStress (matData_);
    else {
      Error("Wrong space dimension of elements! ",__FILE__,__LINE__);
    }  
  
    // Get element solution
    Matrix<Double> elemDisp;
    sol_->GetElemSolutionAsMatrix( elemDisp, ent.GetElem()->connect );

    partElemVec.Resize(nrNodes * nrDofs);
    partElemVec.Init();
  
    elemVec.Resize(nrNodes*nrDofs);
    elemVec.Init();
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {    
      stressBiformInt->CalcStressVec(piolaStressVec, actIntPt,
				     ptelem, ptCoord_, elemDisp);
      stressBiformInt->calcNonLinBMat(nonLinBMat, actIntPt, 
				      ptelem, ptCoord_, elemDisp);
      stressBiformInt->calcLinBMat(linBMat, actIntPt, ptelem, ptCoord_);
      
      nonLinBMat += linBMat;
	
      nonLinBMat.Transpose(transpSumB);
	
      partElemVec = transpSumB * piolaStressVec;
	
      Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_);
	
      if (isaxi_) {
	Vector<Double> shpFncAtIp;
	Vector<Double> coordAtIP;
	ptelem->GetShFncAtIp(shpFncAtIp, actIntPt);
	coordAtIP = ptCoord_ * shpFncAtIp;
	jacDet *= 2 * PI * coordAtIP[0];
      }
      
      partElemVec *= jacDet * intWeights[actIntPt-1];
      
      // the negative sign is due the fact, that this vector has
      // to be subtracted from the RHS!! 
      // (see Kaltenbacher  "Numerical Sim. of Mechatronic Sensors 
      // and Actuators" p. 55
      elemVec -=  partElemVec;
    }

    delete stressBiformInt;
  }




  // ==================================================================
  //  recovery technique. calculation of element matrix for RHS
  // ==================================================================

  RHSForRecoveryProcedure::RHSForRecoveryProcedure() 
  {
    ENTER_FCN( "RHSForRecoveryProcedure::RHSForRecoveryProcedure" );
    name_ = "RHSForRecoveryProcedure";
  }
 
  RHSForRecoveryProcedure::~RHSForRecoveryProcedure()
  {
    ENTER_FCN( "RHSForRecoveryProcedure::~RHSForRecoveryProcedure" );
  }

  void  RHSForRecoveryProcedure::CalcElemVectorRHSForSPR(BaseFE * aptelem,
                                                         Matrix<Double>&
                                                         ptCoord,
                                                         Vector<Double> & fncNodesElem,
                                                         const UInt aComponent,
                                                         Vector<Double> & elemVec)
  {
    ENTER_FCN( "RHSForRecoveryProcedure::CalcElemVectorRHSForSPR" );
             
    const UInt nrIntPnts = ptelem->GetNumIntPoints();
    const UInt nrNodes   = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<double> shFnc;
    Matrix<Double>      dvShFnc ;
    Double              jacDet;  
    Double              valFnc=0;
    UInt             iIntPnts,iShFnc;
  
    elemVec.Resize(nrNodes);
    elemVec.Init();
  
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

  PreStressLinFormInt::PreStressLinFormInt(BaseMaterial* mat, 
                                           Double aPreStressVal, 
                                           Directions stressDir) 
    :nLinMech_linFormInt(mat), 
     preStressVal_(aPreStressVal), preStressDir_(stressDir)
  
  {
    ENTER_FCN( "PreStressLinFormInt::PreStressLinFormInt" );
    name_ = "PreStressLinFormInt";
  }



  PreStressLinFormInt::~PreStressLinFormInt()
  {
    ENTER_FCN( "PreStressLinFormInt::~PreStressLinFormInt" );
  }




  void PreStressLinFormInt::CalcElemVector( Vector<Double> & result,
                                            EntityIterator& ent )
  {
    ENTER_FCN( "PreStressLinFormInt::CalcElemVector" );
    Error( "PreStressLinFormInt::CalcElemVector: not working",
           __FILE__, __LINE__ );

    //   const UInt nrIntPts = ptelem->GetNumIntPoints();
    //   const UInt nrNodes  = ptelem->GetNumNodes();
    //   const UInt nrDofs   = getNrDofs();
    //   const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    //   Vector<Double> piolaStressVec;    
    //   Vector<Double> partElemVec;
  
    //   Matrix<Double> linBMat; 
    //   Matrix<Double> nonLinBMat; 
    //   Matrix<Double> dMat;
    //  // we need transposed of the b-matrices
    //   Matrix<Double> transpSumB; 
  

    //   // This, as friend defined bilinearform holds 
    //   //the necessary differential operators
    //   PreStressInt preStressBiformInt(ptelem, matData_, 
    //                  preStressVal_, preStressDir_);

    //   if (!elemDisp_.GetSizeRow() || !elemDisp_.GetSizeCol()) 
    //     Error("Undefined displacements! ",__FILE__,__LINE__);

    //   partElemVec.Resize(nrNodes * nrDofs);


    //   elemVec.Resize( nrNodes * nrDofs );
    //   elemVec *= 0;    // set elems to 0
   


    //   for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
    //     {    
    //       preStressBiformInt.setActElemDispl(elemDisp_);
    //       preStressBiformInt.CalcStressVec(piolaStressVec, actIntPt,
    //                                                ptCoord);
    //       preStressBiformInt.calcNonLinBMat(nonLinBMat, actIntPt, 
    //                                                ptCoord);
    //       preStressBiformInt.calcLinBMat(linBMat, actIntPt, ptCoord);

    //       nonLinBMat += linBMat;

    //       nonLinBMat.Transpose(transpSumB);
        
    //       partElemVec = transpSumB * piolaStressVec;

    //       Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, 
    //                                                ptCoord);
    //       partElemVec *= jacDet * intWeights[actIntPt-1];
        
    //       elemVec +=  partElemVec;
    //     }
  }



  //==================================================================
  //Members of class LinearFormFlowNoise
  //==================================================================

  LinearFlowNoiseInt::LinearFlowNoiseInt(BaseFE * aptelem) 
  {
    ENTER_FCN( "LinearFlowNoiseInt::LinearFlowNoiseInt" );
    name_ = "LinearFlowNoiseInt";
    ptelem = aptelem;
  }

  LinearFlowNoiseInt::~LinearFlowNoiseInt()
  {
    ENTER_FCN( "LinearFlowNoiseInt::~LinearFlowNoiseInt" );
  }

  void LinearFlowNoiseInt::CalcElemVector4Dip(Matrix<Double>& ptCoord,
                                              const StdVector<UInt> & connecth, 
                                              Vector<Double> & Result, 
                                              const Vector<Double> gradN_x_P)
  {
    ENTER_FCN( "LinearForm::CalcElemVector4FlowSrc" );

    UInt l=ptelem->GetNumIntPoints();
    UInt n=ptelem->GetNumNodes();
  
    UInt i,ii;
  
    Double jacDet;    
    //     Double density=1.0;
    Result.Resize(n);
    Result.Init();
  
    Vector<Double>  Sf;
    Double multiplier;
    // Dipole Source Term: (w,deltp)onbnd
    Vector<Double> normal;
  
    normal.Resize(gradN_x_P.GetSize());
    normal.Init();

  
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

  void LinearFlowNoiseInt::CalcElemVec_withdTijdi(const Matrix<Double>& ptCoord,
                                                  const Matrix<Double>& dTijdi,
                                                  Vector<Double> & Result)
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector_withdTijdi" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
  
    // Source term =  integral(grad[Sf].div[T])

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();

    Matrix<Double> xiDx;      
    Vector<Double>  Sf;

    Vector<Double>  dTijdiAtIP;
    
    Double jacDet;
    Integer actInt;
    Double density=1.204; // At 20C

    Result.Resize(n);
    Result.Init();

    Vector<Double> partResult;
    partResult.Resize(n);
    partResult.Init();
    
    // Loop over all integration points 
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet);

	// dTijdi at integration point: (dTijdx  dTijdy)^T  (2x1)
	dTijdiAtIP =  dTijdi * Sf;
        
        dTijdiAtIP *= density;
        
        // Multiplication with the derivatives of the shape functions
	partResult  = xiDx * dTijdiAtIP;        
        partResult *= jacDet;
        Result     += partResult;
      }

  } // end of method


  void LinearFlowNoiseInt::CalcElemVec_withVortexVel(const Matrix<Double>& ptCoord,
                                                     const Matrix<Double> & NodalVel,
                                                     Vector<Double> & Result)
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector_withVortexVel" << std::endl;
#endif

    //std::cout << "NodalVel:\n" << NodalVel << std::endl;

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
  
    // Source term =  integral(grad[Sf].div[T])

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dimelem = ptCoord.GetSizeRow();

    Matrix<Double> xiDx;      
    Vector<Double>  Sf;

    Vector<Double>  VelAtIP;
    Matrix<Double> VelDerAtIP;

    Double jacDet;
    Integer actInt;
    Double density=1.204; //at 20C

    Result.Resize(n);
    Result.Init();

    Vector<Double> partResult;
    partResult.Resize(n);
    partResult.Init();
    
    Vector<Double> helpVec;
    helpVec.Resize(dimelem);
    helpVec.Init();

    
    // Loop over all integration points 
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet);

	// velocity at integration point: (vx  vy)^T  (2x1)
	VelAtIP = NodalVel * Sf;

	//first derivative of velocity at integration point: (2x2)
	//  vx,x   vx,y
        //  vy,x   vy,y
	//
        VelDerAtIP = NodalVel * xiDx;

	helpVec[0] = 2.0 * VelAtIP[0] * VelDerAtIP[0][0] 
	  + VelAtIP[1] *  VelDerAtIP[0][1] 
	  + VelAtIP[0] *  VelDerAtIP[1][1]; 

	helpVec[1] = 2.0 * VelAtIP[1] * VelDerAtIP[1][1] 
	  + VelAtIP[0] *  VelDerAtIP[1][0] 
	  + VelAtIP[1] *  VelDerAtIP[0][0]; 
	  
	helpVec *= density;
        
        // Multiplication with the derivatives of the shape functions
	partResult  = xiDx * helpVec;        
        partResult *= jacDet;
        Result     += partResult;
      }
    //    std::cout<<"ElemVect computed with test velocities: "<<std::endl;
    //    std::cout<<Result<<std::endl;
  } // end of method

  void LinearFlowNoiseInt::CalcElemVector4Quad(Matrix<Double>& ptCoord,
                                               const StdVector<UInt> & connecth,
                                               const Matrix<Double> & FlowData, 
                                               Vector<Double> & Result)
  {
    ENTER_FCN( "LinearFlowNoiseInt::CalcElemVector4Quad" );

    // This is used if we get the interpolated velocity components with MpCCI 
    // from a fluid computation

    UInt l=ptelem->GetNumIntPoints();
    UInt n=ptelem->GetNumNodes();
    Matrix<Double> xiDx;      
    Matrix<Double> elVel;
    Matrix<Double> elVelAtIP;    
    Double jacDet;
    UInt actInt;
    Double density=1.0;
    Result.Resize(n);
    Result.Init();

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
    VelMatrixforMult.Resize(dimelem,true); 
    VelMatrixforMult.Init();

    Vector<Double> dTij_di;
    VelAtIP.Resize(dimelem);
    VelDerAtIP.Resize(dimelem,true);  
    VelDerAtIP.Init();  
    dTij_di.Resize(dimelem);
    dTij_di.Init();
    Vector<Double> helpVect;

    //     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // //     // For using with Anne Le Duc source
    //     std::vector<Double> ddTij_dxidxi;
    //     ddTij_dxidxi.resize(n);
    //      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    GetQttiesOfElement(elVel, FlowData, connecth, dimelem);
  

    //      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //      // Just for testing with source term from Anne Le Duc
    //      // (ddTij/dxidxj) in MixLayer example
    //     for (actInt=1; actInt<=l; actInt++)
    //       {
    //      // Set to zero index to be filled in.
    //      ddTij_dxidxi[actInt-1]=0; 
    //      ptelem->GetShFncAtIp(Sf, actInt);
    //      for (int ctrIP=1; ctrIP<=l; ctrIP++)
    //        {
    //          // Interpolate to IP and fill in vector
    //          // In files from MixL first value is ddTij/dxidxj!!
    //      ddTij_dxidxi[actInt-1]+=(FlowData[0][connecth[ctrIP-1]-1])*
    //      Sf[actInt-1]; 
    //        }
    //       }
    //      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    elVelAtIP.Resize(dimelem, l);
    elVelAtIP.Init();
    for (actInt=1; actInt<=l; actInt++)
      {
        ptelem->GetShFncAtIp(Sf, actInt);
        //jacDet = ptelem->CalcJacobianDetAtIp(actInt, ptCoord);
        VelAtIP=elVel*Sf;
        
        
        for (int comp=0;comp<dimelem;comp++)
          {
            //Filling the matrix of velocity components at IP
            elVelAtIP[comp][actInt-1]=VelAtIP[comp];
          }
      }
    

    for (actInt=1; actInt<=l; actInt++)
      {
        ptelem->GetShFncAtIp(Sf, actInt);
        ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet);
      
        // This work only for quad1 and trilinear hexahedrals 
        // elements since we have values
        // of the flow quantities only at the corners!
        // Here we compute the derivatives of the Lighthill's tensor 
        // needed in the quadrupole term
        // at the ith integration point
     
        //Implementation 26.09.03
        //Modified 15.06.04

        VelAtIP=elVel*Sf;
        

        //VelDerAtIP=(elVelAtIP*xiDx);
        VelDerAtIP=(elVel*xiDx);

        //VelDerAtIP*=jacDet; 

        VelDerAtIP.GetDiagInMatrix(VelDerFromDiag);
        VelDerFromDiag.ConvertToVec_AppendRows(helpVect);

        //B_vector of derivation
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




        //      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //      // Just for testing with source term from Anne Le Duc 
        //      //(ddTij/dxidxj) in MixLayer example
        //      ddTij_dxidxi*=density;
        //      for (int ii=0;ii<n;ii++)
        //        {
        //          partResult[ii]=Sf[ii]*ddTij_dxidxi[ii];
        //          partResult[ii]*=jacDet;
        //        }
        //      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        Result+=partResult;
      }

  } // end of method
 

  void LinearFlowNoiseInt::GetQttiesOfElement(Matrix<Double>& elVel,
                                              const Matrix<Double>& FlowData,
                                              const StdVector<UInt>& connecth, 
                                              UInt matrixRow)
  {
    ENTER_FCN( "LinearFlowNoiseInt::GetVecOfElement" );
 
    // displacement of element nodes
    elVel.Resize(matrixRow, connecth.GetSize());
   
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++)
      for (UInt dim=0; dim < matrixRow; dim++)
        // dim+1 because index 0 in FlowData is used for storing pressure
        elVel[dim][actNode] = FlowData[dim+1][connecth[actNode]-1]; 
  }


  // ===================================================================
  // nonlinear RHS for nonlinear acoustics
  // ===================================================================

  nLinKuznetsovRHSInt::nLinKuznetsovRHSInt(Double aVal, bool isaxi)
    : LinearForm(), val_(aVal)
  {
    ENTER_FCN( "nLinKuznetsovRHSInt::nLinKuznetsovRHSInt" );
    name_ = "nLinKuznetsovRHSInt";
    isaxi_ = isaxi;
  }


  nLinKuznetsovRHSInt::~nLinKuznetsovRHSInt()
  {
    ENTER_FCN( "nLinKuznetsovRHSInt::~nLinKuznetsovRHSInt" );
  }


  void nLinKuznetsovRHSInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent ) 
  {
    ENTER_FCN( "nLinKuznetsovRHSInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();

    // derivative of shape functions after global coordinates 
    Matrix<Double> xiDx, xiDxTransp;

    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIp;
    Vector<Double> solGradAtIp, solDeriv1GradAtIp;
    Double solDeriv1AtIp, solDeriv2AtIp;

    Double jacDet;
  
    elemVec.Resize(nrNodes);
    elemVec.Init(0.0);

    Double totalfactor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = 0;
      ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);
        
      if (isaxi_) {
        CoordAtIp = ptCoord_ * ShpFncAtIp;
        for (UInt i=0; i<nrNodes; i++)
          xiDx[i][0] += ShpFncAtIp[i] / CoordAtIp[0];
          
        jacDet *= 2 * PI * CoordAtIp[0];
      }
        
      xiDx.Transpose(xiDxTransp);

      //compute gradient of solution + 1st derivative at integration point
      solGradAtIp       = xiDxTransp * sol_;
      solDeriv1GradAtIp = xiDxTransp * solderiv1_;

      //get 1st and 2nd derivative of solution at integration point
      solDeriv1AtIp = solderiv1_*ShpFncAtIp;
      solDeriv2AtIp = solderiv2_*ShpFncAtIp;
        
      totalfactor=0;
      for (UInt j=0; j<xiDx.GetSizeCol(); j++)
        totalfactor += solGradAtIp[j]*solDeriv1GradAtIp[j];
      totalfactor *= factorN2_;
        
      totalfactor += factorN1_ * solDeriv1AtIp * solDeriv2AtIp;
        
      totalfactor *= jacDet;
      for (UInt i=0; i< nrNodes; i++)
        elemVec[i] += ShpFncAtIp[i] * totalfactor;
        
    }
    //  std::cerr << "RHS in linearForm:\n" << elemVec << std::endl;
  }

  nLinWesterveltRHSInt::nLinWesterveltRHSInt(Double aVal, bool isaxi)
    : LinearForm(), val_(aVal)
  {
    ENTER_FCN( "nLinWesterveltRHSInt::nLinWesterveltRHSInt" );
    name_ = "nLinWesterveltRHSInt";
    isaxi_ = isaxi;
  }


  nLinWesterveltRHSInt::~nLinWesterveltRHSInt()
  {
    ENTER_FCN( "nLinWesterveltRHSInt::~nLinWesterveltRHSInt" );
  }


  void nLinWesterveltRHSInt::CalcElemVector( Vector<Double> & elemVec,
                                             EntityIterator& ent )
  {
    ENTER_FCN( "nLinWesterveltRHSInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();

    Vector<Double> ShpFncAtIp;
    Vector<Double> CoordAtIp;
    Double solAtIp, solDeriv1AtIp, solDeriv2AtIp;

    Double jacDet;
  
    elemVec.Resize(nrNodes);
    elemVec.Init(0.0);

    Double totalfactor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_);
      ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
        
      if (isaxi_) {
        CoordAtIp = ptCoord_ * ShpFncAtIp;
        jacDet *= 2 * PI * CoordAtIp[0];
      }
        
      //get solution and derivatives at integration point
      solAtIp = sol_*ShpFncAtIp;
      solDeriv1AtIp = solderiv1_*ShpFncAtIp;
      solDeriv2AtIp = solderiv2_*ShpFncAtIp;
        
      totalfactor = jacDet*factor_*2.0*(solDeriv1AtIp * solDeriv1AtIp +
                                        solAtIp*solDeriv2AtIp);
        
      for (UInt i=0; i< nrNodes; i++)
        elemVec[i] += ShpFncAtIp[i] * totalfactor;
    }
  }


  // ====================================================================
  // electric polarization
  // ====================================================================

  ElecPolarizationInt::ElecPolarizationInt(bool isaxi, bool coordUpdate )
    : LinearForm()
  {
    ENTER_FCN( "ElecPolarizationInt::ElecPolarizationInt" );
    name_ = "ElecPolarizationInt";
    isaxi_ = isaxi;
    coordUpdate_ = coordUpdate;
  }


  ElecPolarizationInt::~ElecPolarizationInt()
  {
    ENTER_FCN( "ElecPolarizationInt::~ElecPolarizationInt" );
  }


  void ElecPolarizationInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent )
  {
    ENTER_FCN( "ElecPolarizationInt::CalcElemVector" );
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> ShpFncAtIp, CoordAtIP;
    Matrix<Double> xiDx;

    elemVec.Resize(nrNodes);
    elemVec.Init(0);

    Double factor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = 0;
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);

      if (isaxi_) {
        ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
        CoordAtIP = ptCoord_ * ShpFncAtIp;
        for (UInt i=0; i<nrNodes; i++)
          xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
      
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      factor = intWeights[actIntPt-1] * jacDet;
      for (UInt i=0; i<nrNodes; i++) {
        for (UInt j=0; j<D_.GetSize(); j++) {
          elemVec[i] += xiDx[i][j] * D_[j];
        }
      }
    }

  }


  // ==================================================================
  // piezoelectric polarization
  // ==================================================================

  PiezoPolarizationInt::PiezoPolarizationInt(UInt dir, UInt numdof, 
                                             bool isaxi)
    : LinearForm(), comp_(dir-1), numDofs_(numdof) 
  {
    ENTER_FCN( "PiezoPolarizationInt::PiezoPolarizationInt" );
    name_ = "PiezoPolarizationInt";
    isaxi_ = isaxi;
    Pval_  = 0.0;
  }


  PiezoPolarizationInt::~PiezoPolarizationInt()
  {
    ENTER_FCN( "PiezoPolarizationInt::~PiezoPolarizationInt" );
  }


  void PiezoPolarizationInt::CalcElemVector( Vector<Double> & elemVec,
                                             EntityIterator& ent ) 
    
  {
    ENTER_FCN( "PiezoPolarizationInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> ShpFncAtIp, partElemVec, CoordAtIP;
    Matrix<Double> xiDx;

    partElemVec.Resize(nrNodes);
    partElemVec.Init(0);  

    Double factor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = 0;
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, jacDet);

      if (isaxi_) {
        if ( comp_ == 0 ) {
          ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
          CoordAtIP = ptCoord_ * ShpFncAtIp;
          for (UInt i=0; i<nrNodes; i++)
            xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
        }
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      factor = intWeights[actIntPt-1] * jacDet * Pval_;
      for (UInt i=0; i<nrNodes; i++) {
        partElemVec[i] += xiDx[i][comp_] * factor;
      }
    }

    // std::cerr << "rhs=" << partElemVec << std::endl; 

    //compute element vector by correctly putting the components of
    //partlementVec into elemVec, since the values for the mechanical
    //degree of freedoms are zero
    elemVec.Resize(nrNodes*numDofs_);
    elemVec.Init(0);  
    UInt idx;
    for (UInt i=0; i<nrNodes; i++) {
      idx = numDofs_*i + numDofs_ - 1;
      elemVec[idx] = partElemVec[i];
    }

  }


  // =========================================================================
  // mechanic rhs volume integrator
  // =========================================================================
  
  
  MechVolForceInt::MechVolForceInt(UInt numDof, bool isaxi) {

    ENTER_FCN( "MechVolForceInt::MechVolForceInt" );

    name_ = "MechVolForceInt";
    isaxi_ = isaxi;
    numDofs_ = numDof;

  }
    
  MechVolForceInt::~MechVolForceInt() {

    ENTER_FCN( "MechVolForceInt::~MechVolForceInt" );

  }
    
  void MechVolForceInt::SetVolForceVector(StdVector<std::string> & volForce, 
                                          const CoordSystem * coordSys,
                                          bool isUnit, Double volume ) {

    ENTER_FCN( "MechVolForceInt::SetVolForceVector" );

    locForce_ = volForce;
    coordSys_ = coordSys;
    isUnitValue_ = isUnit;
    volume_ = volume;
    
  }
    
  void MechVolForceInt::CalcElemVector( Vector<Double> & elemVec,
                                        EntityIterator& ent ) {

    ENTER_FCN( "MechVolForceInt::CalcElemVector");

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc, CoordAtIP, globVec;

    // get global coordinate system and math parser
    MathParser * parser = domain->GetMathParser();

    // First, map force to global coordinate system
    Vector<Double> globMidPoint, locMidPoint;

    ptelem->GetCoordMidPoint(locMidPoint);
    ptelem->Local2GlobalCoord(globMidPoint, locMidPoint, ptCoord_);

    // Update variables for mathParser
    parser->SetCoordinates( mHandle_, *coordSys_, globMidPoint );

    // Now evaluate each entry 
    Vector<Double> locLoadVec( numDofs_ );
    for ( UInt i = 0; i < locForce_.GetSize(); i++ ) {
      parser->SetExpr( mHandle_, locForce_[i] );
      locLoadVec[i] = parser->Eval( mHandle_ );
    }

    // If load is not unit load, divide by volume
    if ( isUnitValue_ == false ) {
      locLoadVec /= volume_;
    }
    
    // Map local load vector to global one
    coordSys_->Local2GlobalVector(globVec, locLoadVec,  globMidPoint);

    // Then, calculate element vector
    elemVec.Resize(nrNodes * numDofs_);
    elemVec.Init(0.0);
    Double factor;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_);
      ptelem->GetShFncAtIp(shapeFnc,actIntPt);

      if (isaxi_) {
        CoordAtIP = ptCoord_ * shapeFnc;
        jacDet *= 2 * PI * CoordAtIP[0];
      }
          
      factor = intWeights[actIntPt-1] * jacDet;
      for (UInt iNode = 0; iNode < nrNodes; iNode++) {
        for (UInt iDof = 0; iDof < numDofs_; iDof++) {
          elemVec[iNode*numDofs_ + iDof] += 
            factor * shapeFnc[iNode] * globVec[iDof]; 
        }
      }
    }
  }



  // ==================================================================
  // linMech: add mechanical stress as force to RHS
  // ==================================================================

  AddStressRHSInt::AddStressRHSInt(BaseMaterial* matData,
				     Vector<Double>& stressVec,
				     SubTensorType type)
    : LinearForm(), matData_(matData)
  {
    ENTER_FCN( "AddStressRHSInt::AddStressRHSInt" );
    name_ = "AddStressRHSInt";
    if ( type == AXI ) 
      isaxi_ = true;

    subTensorType_ = type;
    if ( type == AXI ) {
      addStress_.Resize(4);
      addStress_.Init(0.0);
      addStress_[0] = stressVec[0];
      addStress_[1] = stressVec[1];
      addStress_[2] = stressVec[5];
      addStress_[3] = stressVec[2];
    }
    else if ( type == PLANE_STRAIN ) {
      addStress_.Resize(3);
      addStress_.Init(0.0);
      addStress_[0] = stressVec[0];
      addStress_[1] = stressVec[1];
      addStress_[2] = stressVec[5];

    }
    else if ( type == FULL ) {
      addStress_.Resize(6);
      addStress_.Init(0.0);
      addStress_[0] = stressVec[0];
      addStress_[1] = stressVec[1];
      addStress_[2] = stressVec[2];
      addStress_[3] = stressVec[3];      
      addStress_[4] = stressVec[4];
      addStress_[5] = stressVec[5];
    }
  }

  AddStressRHSInt::~AddStressRHSInt()
  {
    ENTER_FCN( "AddStressRHSInt::~AddStressRHSInt" );
  }

  void AddStressRHSInt::CalcElemVector( Vector<Double> & elemVec,
					 EntityIterator& ent )
  {
    ENTER_FCN( "AddStressRHSInt::CalcElemVector" );

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();

    const UInt nrDofs   = ptelem->GetDim(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> piolaStressVec;    
    Vector<Double> partElemVec;
  
    Matrix<Double> linBMat; 
    Matrix<Double> transpB;    // we need transposed of the b-matrices
  
    linElastInt* bilinearStiff = new linElastInt(matData_, subTensorType_);

    partElemVec.Resize(nrNodes * nrDofs);
    partElemVec.Init();
  
    elemVec.Resize(nrNodes*nrDofs);
    elemVec.Init(0.0);
  
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
        bilinearStiff->calcBMatOnly(linBMat, actIntPt, ptelem, ptCoord_);
      
	linBMat.Transpose(transpB);
      
        partElemVec = transpB * addStress_;
      
        Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_);
      
        if (isaxi_) {
	  Vector<Double> shpFncAtIp;
	  Vector<Double> coordAtIP;
	  ptelem->GetShFncAtIp(shpFncAtIp, actIntPt);
	  coordAtIP = ptCoord_ * shpFncAtIp;
	  jacDet *= 2 * PI * coordAtIP[0];
	}
      
        partElemVec *= jacDet * intWeights[actIntPt-1];
	elemVec +=  partElemVec;
      }

    //    std::cout << "\n  Forcfe:\n" <<  elemVec << std::endl;
    delete bilinearStiff;
  }

} // end of namespace
