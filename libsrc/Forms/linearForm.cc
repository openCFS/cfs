#include <iostream>
#include <fstream>
#include <math.h>

#include "forms_header.hh"
#include "Utils/coordSystem.hh"
#include "Domain/domain.hh"

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


  void LinearForm::CalcElemVector(Matrix<Double>& ptCoord, 
                                  Vector<Double> & Result)
  {
    ENTER_FCN( "LinearForm::CalcElemVector" );
  }


  // ================================================================
  // edge integration
  // ================================================================

  LinearEdgeInt::LinearEdgeInt(BaseFE * aptelem, Double aVal, 
                               UInt aDirection,
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
  
  
    // set vector to desired size and set all elements to zero    
    elemVec.Resize(nrEdges); 
    for (UInt i=0; i<nrEdges; i++)
      elemVec[i] = 0;
  
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
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord);
            currentVec[xDir] = -(globCoord[yDir] - (*coilMidPt_)[yDir]);
            currentVec[yDir] = globCoord[xDir] - (*coilMidPt_)[xDir];
            currentVec[zDir] = 0;
          
            len = currentVec.NormL2();
            currentVec *= val_ / len;
            break;

          case 5: // current in yz-plane ==> x=0
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord);
            currentVec[xDir] = 0;
            currentVec[yDir] = -(globCoord[zDir] - (*coilMidPt_)[zDir]);
            currentVec[zDir] = globCoord[yDir] - (*coilMidPt_)[yDir];

            len = currentVec.NormL2();
            currentVec *= val_ / len;
            break;

          case 6: // current in xz-plane ==> y=0
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord);
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
        
        ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord);

        partElemVec = shapeEdge * currentVec;

        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);

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


  void VolumeSrcInt::CalcElemVector(Matrix<Double>& ptCoord, 
                                    Vector<Double> & elemVec)
  {
    ENTER_FCN( "VolumeSrcInt::CalcElemVector" );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc;
  
    elemVec.Resize(nrNodes);
    elemVec *= 0;    // set elems to 0
  
    Double factor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {  
        ptelem->GetShFncAtIp(shapeFnc,actIntPt);
        factor = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord) * 
          intWeights[actIntPt-1] * val_;
      
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


  // ===================================================================
  // permanent magnet in 2D
  // ===================================================================

  MagPerm2DInt::MagPerm2DInt(Vector<Double> vecVal, Double rel, 
                             Boolean isaxi)
    : LinearForm(), perm_(vecVal), reluctivity_(rel)
  {
    ENTER_FCN( "MagPerm2DInt::VolumeSrcInt" );
    isaxi_ = isaxi;
  }


  MagPerm2DInt::~MagPerm2DInt()
  {
    ENTER_FCN( "MagPerm2DInt::~MagPerm2DInt" );
  }


  void MagPerm2DInt::CalcElemVector(Matrix<Double>& ptCoord, 
                                    Vector<Double> & elemVec)
  {
    ENTER_FCN( "MagPerm2DInt::CalcElemVector" );

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
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

      if (isaxi_)
        {
          ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
          CoordAtIP = ptCoord * ShpFncAtIp;
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

  nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt(BaseFE * aptelem, 
                                                     MaterialData & matData,
                                                     Boolean isaxi) 
    : LinearForm(aptelem), matData_(matData)
  {
    ENTER_FCN( "nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt" );
    isaxi_ = isaxi;
  }

  nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt(ApproxData *nlinFnc, 
                                                     Double startVal, 
                                                     Boolean axi)
    : LinearForm()
  {
    ENTER_FCN( "nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt" );
  
    isaxi_       = axi;
    startmatVal_ = startVal;
    nlinFnc_     = nlinFnc;
  }
  
  nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt(Double startVal, 
                                                     Boolean axi)
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

  void nLinMagNode2D_linFormInt::CalcElemVector(Matrix<Double>& ptCoord,
                                                Vector<Double> & elemVec)
  {
    ENTER_FCN("nLinMagNode2D_linFormInt :: ~CalcElemVector" );
  
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
        
      //set the element solution vector to the bilinearform
      curlcurl2D->SetActElemSol(magPotinMatrix_);
    }

    //set the element-pointer
    curlcurl2D->SetElemPtr(ptelem);
  
    Matrix<Double> elemmat;
    curlcurl2D->CalcElementMatrix(ptCoord, elemmat);
  
    elemVec.Resize(nrNodes);
    elemVec = -(elemmat * magPot_);

    delete curlcurl2D;
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


  nLinMech_linFormInt::nLinMech_linFormInt(MaterialData & matData,
                                           Boolean isaxi) 
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
      stressBiformInt = new nLinMechPlaneStrainInt_PiolaStress(ptelem,
                                                               matData_);
    else if (ptelem->GetDim() == 2 && isaxi_)
      stressBiformInt = new nLinMechAxiInt_PiolaStress(ptelem, matData_);
    else if (ptelem->GetDim() == 3)
      stressBiformInt = new nLinMech3dInt_PiolaStress (ptelem, matData_);
    else {
      Error("Wrong space dimension of elements! ",__FILE__,__LINE__);
    }  
  
  
    if (!elemDisp_.GetSizeRow() || !elemDisp_.GetSizeCol()) 
      Error("Undefined displacements! ",__FILE__,__LINE__);
  
    partElemVec.Resize(nrNodes * nrDofs);
  
  
    elemVec.Resize(nrNodes*nrDofs);
    elemVec *= 0;    // set elems to 0
  
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
        stressBiformInt->setActElemDispl(elemDisp_);
        stressBiformInt->CalcStressVec(piolaStressVec, actIntPt,ptCoord);
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

  RHSForRecoveryProcedure::RHSForRecoveryProcedure(BaseFE * aptelem) 
    : LinearForm(aptelem)
  {
    ENTER_FCN( "RHSForRecoveryProcedure::RHSForRecoveryProcedure" );
  }
 
  RHSForRecoveryProcedure::~RHSForRecoveryProcedure()
  {
    ENTER_FCN( "RHSForRecoveryProcedure::~RHSForRecoveryProcedure" );
  }

  void  RHSForRecoveryProcedure::CalcElemVectorRHSForSPR(Matrix<Double>&
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
    : LinearForm(aptelem)
  {
    ENTER_FCN( "LinearFlowNoiseInt::LinearFlowNoiseInt" );
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
    Integer dimelem = ptCoord.GetSizeRow();

    Matrix<Double> xiDx;      
    Vector<Double>  Sf;

    Vector<Double>  dTijdiAtIP;
    
    Double jacDet;
    Integer actInt;
    Double density=1.0;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
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
    VelDerAtIP.Resize(dimelem);  

    Double jacDet;
    Integer actInt;
    Double density=1.0;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
    Vector<Double> helpVec;
    helpVec.Resize(dimelem);

    
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

    for (UInt i=0;i<n;i++)
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


  void nLinKuznetsovRHSInt::CalcElemVector(Matrix<Double>& ptCoord,
                                           Vector<Double> & elemVec)
  {
    ENTER_FCN( "nLinKuznetsovRHSInt::CalcElemVector" );

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
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);
        
      if (isaxi_) {
        CoordAtIp = ptCoord * ShpFncAtIp;
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

  nLinWesterveltRHSInt::nLinWesterveltRHSInt(Double aVal, Boolean isaxi)
    : LinearForm(), val_(aVal)
  {
    ENTER_FCN( "nLinWesterveltRHSInt::nLinWesterveltRHSInt" );
    isaxi_ = isaxi;
  }


  nLinWesterveltRHSInt::~nLinWesterveltRHSInt()
  {
    ENTER_FCN( "nLinWesterveltRHSInt::~nLinWesterveltRHSInt" );
  }


  void nLinWesterveltRHSInt::CalcElemVector(Matrix<Double>& ptCoord,
                                            Vector<Double> & elemVec)
  {
    ENTER_FCN( "nLinWesterveltRHSInt::CalcElemVector" );

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

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
      ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
        
      if (isaxi_) {
        CoordAtIp = ptCoord * ShpFncAtIp;
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

  ElecPolarizationInt::ElecPolarizationInt(Boolean isaxi)
    : LinearForm() 
  {
    ENTER_FCN( "ElecPolarizationInt::ElecPolarizationInt" );
    isaxi_ = isaxi;
  }


  ElecPolarizationInt::~ElecPolarizationInt()
  {
    ENTER_FCN( "ElecPolarizationInt::~ElecPolarizationInt" );
  }


  void ElecPolarizationInt::CalcElemVector(Matrix<Double>& ptCoord, 
                                           Vector<Double> & elemVec)
  {
    ENTER_FCN( "ElecPolarizationInt::CalcElemVector" );

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
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

      if (isaxi_) {
        ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
        CoordAtIP = ptCoord * ShpFncAtIp;
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
                                             Boolean isaxi)
    : LinearForm(), comp_(dir-1), numDofs_(numdof) 
  {
    ENTER_FCN( "PiezoPolarizationInt::PiezoPolarizationInt" );
    isaxi_ = isaxi;
    Pval_  = 0.0;
  }


  PiezoPolarizationInt::~PiezoPolarizationInt()
  {
    ENTER_FCN( "PiezoPolarizationInt::~PiezoPolarizationInt" );
  }


  void PiezoPolarizationInt::CalcElemVector(Matrix<Double>& ptCoord,
                                            Vector<Double> & elemVec)
  {
    ENTER_FCN( "PiezoPolarizationInt::CalcElemVector" );

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
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord, jacDet);

      if (isaxi_) {
        if ( comp_ == 0 ) {
          ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt);
          CoordAtIP = ptCoord * ShpFncAtIp;
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
  
  
  MechVolForceInt::MechVolForceInt(UInt numDof, Boolean isaxi) {

    ENTER_FCN( "MechVolForceInt::MechVolForceInt" );
    
    isaxi_ = isaxi;
    numDofs_ = numDof;

    // Obtain mathParser handler
    mHandler_ = domain->GetMathParser()->GetNewHandler();

  }
    
  MechVolForceInt::~MechVolForceInt() {

    ENTER_FCN( "MechVolForceInt::~MechVolForceInt" );

  }
    
  void MechVolForceInt::SetVolForceVector(StdVector<std::string> & volForce, 
                                          const CoordSystem * coordSys,
                                          Boolean isUnit, Double volume ) {

    ENTER_FCN( "MechVolForceInt::SetVolForceVector" );

    locForce_ = volForce;
    coordSys_ = coordSys;
    isUnitValue_ = isUnit;
    volume_ = volume;
    
  }
    
  void MechVolForceInt::CalcElemVector(Matrix<Double>& ptCoord, 
                                       Vector<Double> & elemVec) {

    ENTER_FCN( "MechVolForceInt::CalcElemVector");

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt nrNodes  = ptelem->GetNumNodes();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc, CoordAtIP, globVec;

    // get global coordinate system and math parser
    MathParser * parser = domain->GetMathParser();

    // First, map force to global coordinate system
    Vector<Double> globMidPoint, locMidPoint;

    ptelem->GetCoordMidPoint(locMidPoint);
    ptelem->Local2GlobalCoord(globMidPoint, locMidPoint, ptCoord);

    // Update variables for mathParser
    parser->SetCoordinates( mHandler_, *coordSys_, globMidPoint );

    // Now evaluate each entry 
    Vector<Double> locLoadVec( numDofs_ );
    for ( UInt i = 0; i < locForce_.GetSize(); i++ ) {
      locLoadVec[i] = parser->Eval( mHandler_, locForce_[i] );
    }

    // If load is not unit load, divide by volume
    if ( isUnitValue_ == FALSE ) {
      locLoadVec /= volume_;
    }
    
    // Map local load vector to global one
    coordSys_->Local2GlobalVector(globVec, locLoadVec,  globMidPoint);

    // Then, calculate element vector
    elemVec.Resize(nrNodes * numDofs_);
    elemVec.Init(0.0);
    Double factor;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord);
      ptelem->GetShFncAtIp(shapeFnc,actIntPt);

      if (isaxi_) {
        CoordAtIP = ptCoord * shapeFnc;
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


} // end of namespace
