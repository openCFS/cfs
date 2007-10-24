// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <math.h>

#include "forms_header.hh"
#include "Utils/coordSystem.hh"
#include "Domain/domain.hh"
#include "Utils/SmoothSpline.hh"

namespace CoupledField {


  LinearForm::LinearForm( BaseMaterial* matData) : BaseForm( matData )
  {
  }




  LinearForm::~LinearForm()
  {
  }


  void LinearForm::CalcElemVector( Vector<Double> & result,
                                   EntityIterator& ent ) {
  }

  void LinearForm::CalcElemVector( Vector<Complex> & result,
                                   EntityIterator& ent ) {

    Vector<Double> helpVec;
    CalcElemVector( helpVec, ent);
    
    Complex amplitude = Complex( 1.0, 0.0 );
    result = helpVec * amplitude;
  }

  // ================================================================
  // edge integration
  // ================================================================

  LinearEdgeInt::LinearEdgeInt( const std::string&  aVal, 
                               UInt aDirection,
                               Vector<Double> * aCoilMidPt) 
    : direction_(aDirection)
  {
    name_ = "LinearEdgeInt";
  
    if(aCoilMidPt)
      coilMidPt_ = new Vector<Double>(*aCoilMidPt);
    mParser_->SetExpr( mHandle_, aVal );
  }

  LinearEdgeInt::~LinearEdgeInt()
  {
  }



  void LinearEdgeInt::CalcElemVector( Vector<Double> & elemVec,
                                      EntityIterator& ent ) 
  {
  
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
  
    Double val = mParser_->Eval( mHandle_ );
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {    
        // define direction of current
        switch(direction_)
          {       
          case 1:
          case 2:
          case 3:
            // if current direction goes parallel to a coordinate axis
            currentVec[direction_-1] = val; 
            break;
          
          case 4: // current in xy-plane ==> z=0
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord_);
            currentVec[xDir] = -(globCoord[yDir] - (*coilMidPt_)[yDir]);
            currentVec[yDir] = globCoord[xDir] - (*coilMidPt_)[xDir];
            currentVec[zDir] = 0;
          
            len = currentVec.NormL2();
            currentVec *= val / len;
            break;

          case 5: // current in yz-plane ==> x=0
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord_);
            currentVec[xDir] = 0;
            currentVec[yDir] = -(globCoord[zDir] - (*coilMidPt_)[zDir]);
            currentVec[zDir] = globCoord[yDir] - (*coilMidPt_)[yDir];

            len = currentVec.NormL2();
            currentVec *= val / len;
            break;

          case 6: // current in xz-plane ==> y=0
            ptelem->GetGlobalEdgeIndicesAtIP(globCoord,actIntPt,ptCoord_);
            currentVec[xDir] = globCoord[zDir] - (*coilMidPt_)[zDir];
            currentVec[yDir] = 0;
            currentVec[zDir] = -(globCoord[xDir] - (*coilMidPt_)[xDir]);

            len = currentVec.NormL2();
            currentVec *= val / len;

            break;
            
          default:
            std::string errMsg = "Selected current direction with num. ";
            errMsg += direction_ + " not supported!";
            Error(errMsg.c_str(),__FILE__,__LINE__);
          }
        
        ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord_);

        partElemVec = shapeEdge * currentVec;

        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem());

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

  VolumeSrcInt::VolumeSrcInt( const std::string& aVal, 
                              bool isaxi, bool coordUpdate )
    : LinearForm()
  {
    name_ = "VolumeSrcInt";
    isaxi_ = isaxi;
    coordUpdate_ = coordUpdate;
    mParser_->SetExpr( mHandle_, aVal );
  }


  VolumeSrcInt::~VolumeSrcInt()
  {
  }

  void VolumeSrcInt::SetFactor( const std::string& factor ) {
    mParser_->SetExpr( mHandle_, factor );
  }
  

  void VolumeSrcInt::CalcElemVector( Vector<Double> & elemVec,
                                     EntityIterator& ent ) 
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc;
  
    elemVec.Resize(numFncs);
    elemVec.Init();
  
    Double factor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
      {  
        ptelem->GetShFncAtIp(shapeFnc,actIntPt, ent.GetElem() );
        factor = 
          ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem()) * 
          intWeights[actIntPt-1] * mParser_->Eval( mHandle_ );
      
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
    name_ = "MagPerm2DInt";
    isaxi_ = isaxi;
    coordUpdate_ = coordUpdate;
  }


  MagPerm2DInt::~MagPerm2DInt()
  {
  }


  void MagPerm2DInt::CalcElemVector( Vector<Double> & elemVec,
                                     EntityIterator& ent )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> ShpFncAtIp, helpVec, CoordAtIP;
    Matrix<Double> xiDx;

    elemVec.Resize(numFncs);
    elemVec.Init(0);  
    helpVec.Resize(numFncs);
    helpVec.Init(0);
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = 0;
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, 
                                    ptCoord_, jacDet, ent.GetElem() );

      if (isaxi_)
        {
          ptelem->GetShFncAtIp(ShpFncAtIp,actIntPt, ent.GetElem());
          CoordAtIP = ptCoord_ * ShpFncAtIp;
          for (UInt i=0; i<numFncs; i++)
            xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
        
          jacDet *= 2 * PI * CoordAtIP[0];
        }

      helpVec = xiDx * perm_;
      helpVec *=  intWeights[actIntPt-1] * jacDet * reluctivity_;
      elemVec += helpVec;
    }
  }

  // ===================================================================
  // permanent magnet in 3D
  // ===================================================================

  MagPerm3DInt::MagPerm3DInt(Vector<Double> vecVal, Double rel, 
                             bool isaxi, bool coordUpdate )
    : LinearForm(), 
      perm_(vecVal), 
      reluctivity_(rel)
  {
    name_ = "MagPerm3DInt";
    isaxi_ = false;
    coordUpdate_ = coordUpdate;
    nrDofs_ = 3;

    curlOp_ = new CurlCurlNode3DInt( NULL, coordUpdate );
  }


  MagPerm3DInt::~MagPerm3DInt()
  {
  }


  void MagPerm3DInt::CalcElemVector( Vector<Double>& elemVec,
                                     EntityIterator& ent )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    // Extract pointer to reference element for curl-bilinear form
    curlOp_->ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    Vector<Double> CoordAtIP;
    Matrix<Double> bMatCurl, bMatDiv, bMatCurlTrans; 
    Vector<Double> helpVec;
    Double jacDet, factor;

    elemVec.Resize( numFncs * nrDofs_ );
    helpVec.Resize( numFncs * nrDofs_ );
    elemVec.Init(0);  
    helpVec.Init(0);
  
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      // Setup the B matrix for current integration point
      curlOp_->calcBMat( bMatCurl, bMatDiv, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent.GetElem() );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
	(*error) << "CurlCurlNode3DInt::CalcElementMatrix: Encountered "
		 << "negative Jacobian determinant!";
	Error( __FILE__, __LINE__ );
      }

      factor = intWeights[actIntPt-1] * jacDet * reluctivity_;
      bMatCurl.Transpose(bMatCurlTrans);  

      //      std::cout << "bMatT:\n" << bMatCurlTrans << std::endl;
      helpVec  = bMatCurlTrans * perm_;
      helpVec *= intWeights[actIntPt-1] * jacDet * reluctivity_;
      elemVec += helpVec;
    }
//      std::cout << "Mvec:\n" << elemVec << std::endl;
//      std::cout << "M:\n" << perm_ << std::endl;
//      std::cout << "rel:\n" << reluctivity_ << std::endl;

  }


  // ==================================================================
  // nLinMagnetics
  // ==================================================================

  nLinMagNode2D_linFormInt::nLinMagNode2D_linFormInt( BaseMaterial* matData, 
                                                      bool axi,
                                                      bool coordUpdate)
    : LinearForm( matData )
  {

    name_ = "nLinMagNode2D_linFormInt";
    isSolDependent_ = true;

    isaxi_       = axi;
    coordUpdate_ = coordUpdate;
    ptMaterial->GetScalar( startmatVal_, MAG_RELUCTIVITY,REAL);

    // need nonlinear BH curve approximation
    if (  ptMaterial->GetNonlinFileName() != "" )  
      ptMaterial->NeedApproxMatCurve( magBH );
  }
  
  nLinMagNode2D_linFormInt::~nLinMagNode2D_linFormInt()
  {
  }

  void nLinMagNode2D_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                                 EntityIterator& ent )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );

    // get pointer to nonlinear BH curve approximation
    nlinFnc_ = ptMaterial->GetNonlinFncBH();

    BaseForm * curlcurl2D;
    if (nlinFnc_== NULL) 
      //define the linear element matrix
      curlcurl2D = new CurlCurlNode2DInt( ptMaterial, isaxi_, coordUpdate_ );
    else {
      //define the nonlinear element matrix
      curlcurl2D = new nLinCurlCurlNode2DInt( ptMaterial, isaxi_, coordUpdate_ );
      //important to set method to FixPoint, since we compute the RHS!!
      curlcurl2D->SetNonLinMethod("fixPoint");
        
      //set the solution class to the operator
      curlcurl2D->SetSolution( * sol_ );
    }

    Matrix<Double> elemmat;
    curlcurl2D->CalcElementMatrix(elemmat, ent, ent);

    // Get element solution
    Vector<Double> magPot;
    sol_->GetElemSolution( magPot, ent  );

    elemVec.Resize(numFncs);
    elemVec = -(elemmat * magPot);

    delete curlcurl2D;
  }


  nLinMagNode3D_linFormInt::nLinMagNode3D_linFormInt( BaseMaterial*matData, 
                                                      bool coordUpdate)
    : LinearForm( matData )
  {
    name_ = "nLinMagNode3D_linFormInt";
    isSolDependent_ = true;

    isaxi_       = false;
    coordUpdate_ = coordUpdate;
    ptMaterial->GetScalar( startmatVal_, MAG_RELUCTIVITY,REAL);

    // need nonlinear BH curve approximation
    // need nonlinear BH curve approximation
    if (  ptMaterial->GetNonlinFileName() != "" )  
      ptMaterial->NeedApproxMatCurve( magBH );
  }
  
  nLinMagNode3D_linFormInt::~nLinMagNode3D_linFormInt()
  {
  }

  void nLinMagNode3D_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                                 EntityIterator& ent )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );

    // get pointer to nonlinear BH curve approximation
    nlinFnc_ = ptMaterial->GetNonlinFncBH();

    BaseForm * curlcurl3D;
    if (nlinFnc_== NULL) 
      //define the linear element matrix
      curlcurl3D = new CurlCurlNode3DInt( ptMaterial, coordUpdate_);
    else {
      //define the nonlinear element matrix
      curlcurl3D = new nLinCurlCurlNode3DInt( ptMaterial, coordUpdate_ );
      //important to set method to FixPoint, since we compute the RHS!!
      curlcurl3D->SetNonLinMethod("fixPoint");
        
      //set the solution class to the operator
      curlcurl3D->SetSolution( *sol_ );
    }

    // Get element solution
    Vector<Double> magPot;
    sol_->GetElemSolution( magPot, ent  );

    Matrix<Double> elemmat;
    curlcurl3D->CalcElementMatrix(elemmat, ent, ent);

    elemVec = -(elemmat * magPot);

    delete curlcurl3D;
  }


  // ==================================================================
  // nLinMech
  // ==================================================================



  nLinMech_linFormInt::nLinMech_linFormInt(BaseMaterial* matData,
                                           bool isaxi) 
    : LinearForm(), matData_(matData)
  {
    name_ = "nLinMech_linFormInt";
    isaxi_ = isaxi;
    isSolDependent_ = true;
  }



  nLinMech_linFormInt::~nLinMech_linFormInt()
  {
  }



  void nLinMech_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    
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
    
    stressBiformInt->SetActEntities( ent, ent );

    // Get element solution
    Matrix<Double> elemDisp;
    sol_->GetElemSolutionAsMatrix( elemDisp, ent  );

    partElemVec.Resize(numFncs * nrDofs);
    partElemVec.Init();
    
  
    elemVec.Resize(numFncs*nrDofs);
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
      
      Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, 
                                                  ptCoord_, ent.GetElem() );
      
      if (isaxi_)
        {
          Vector<Double> shpFncAtIp;
          Vector<Double> coordAtIP;
          ptelem->GetShFncAtIp(shpFncAtIp, actIntPt, ent.GetElem());
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
    name_ = "RHSForRecoveryProcedure";
  }
 
  RHSForRecoveryProcedure::~RHSForRecoveryProcedure()
  {
  }

  void  RHSForRecoveryProcedure::CalcElemVectorRHSForSPR(BaseFE * aptelem,
                                                         Matrix<Double>&
                                                         ptCoord,
                                                         Vector<Double> & fncNodesElem,
                                                         const UInt aComponent,
                                                         Vector<Double> & elemVec)
  {

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPnts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<double> shFnc;
    Matrix<Double>      dvShFnc ;
    Double              jacDet;  
    Double              valFnc=0;
    UInt             iIntPnts,iShFnc;
  
    elemVec.Resize(numFncs);
    elemVec.Init();
  
    for (iIntPnts=0; iIntPnts < nrIntPnts; iIntPnts ++)
      {
        valFnc = 0;
      
        jacDet = ptelem->CalcJacobianDetAtIp(iIntPnts+1, ptCoord, NULL );
      
        ptelem->GetGlobDerivShFncAtIp(dvShFnc, iIntPnts+1, 
                                      ptCoord, jacDet, NULL );
      
        ptelem->GetShFncAtIp(shFnc, iIntPnts+1, NULL  );
      
        for (iShFnc = 0; iShFnc < numFncs; iShFnc ++)
          valFnc += fncNodesElem[iShFnc]*dvShFnc[iShFnc][aComponent];
      
        for (iShFnc = 0; iShFnc < numFncs; iShFnc++) 
          elemVec[iShFnc]+=jacDet*valFnc*shFnc[iShFnc]*intWeights[iIntPnts];
      
      } // loop over integration points
  
  }

  // ==================================================================
  // prestress linearform
  // ==================================================================

  PreStressLinFormInt::PreStressLinFormInt(BaseMaterial* mat, 
                                           const std::string & aPreStressVal, 
                                           Directions stressDir) 
    :nLinMech_linFormInt(mat), 
     preStressDir_(stressDir)
  
  {
    name_ = "PreStressLinFormInt";
    mParser_->SetExpr( mHandle_, aPreStressVal );
  }



  PreStressLinFormInt::~PreStressLinFormInt()
  {
  }




  void PreStressLinFormInt::CalcElemVector( Vector<Double> & result,
                                            EntityIterator& ent )
  {
    Error( "PreStressLinFormInt::CalcElemVector: not working",
           __FILE__, __LINE__ );

    //   const UInt nrIntPts = ptelem->GetNumIntPoints();
    //   const UInt numFncs  = ptelem->GetNumNodes();
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

    //   partElemVec.Resize(numFncs * nrDofs);


    //   elemVec.Resize( numFncs * nrDofs );
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
    name_ = "LinearFlowNoiseInt";
    ptelem = aptelem;
  }

  LinearFlowNoiseInt::~LinearFlowNoiseInt()
  {
  }

  void LinearFlowNoiseInt::CalcElemVector4Dip(Matrix<Double>& ptCoord,
                                              const StdVector<UInt> & connecth, 
                                              Vector<Double> & Result, 
                                              const Vector<Double> gradN_x_P)
  {

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
        jacDet = ptelem->CalcJacobianDetAtIp(i+1, ptCoord, NULL );
      
        ptelem->GetShFncAtIp(Sf, i+1, NULL );
      
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
	ptelem->GetShFncAtIp(Sf, actInt, NULL );
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, 
                                      jacDet, NULL );

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
    Double density=1.0; //1.0 since all properties of material for vortex are 1.0

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
	ptelem->GetShFncAtIp(Sf, actInt, NULL  );
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, 
                                      jacDet, NULL );

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
        ptelem->GetShFncAtIp(Sf, actInt, NULL );
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
        ptelem->GetShFncAtIp(Sf, actInt, NULL );
        ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, 
                                      jacDet, NULL );
      
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
 
    // displacement of element nodes
    elVel.Resize(matrixRow, connecth.GetSize());
   
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++)
      for (UInt dim=0; dim < matrixRow; dim++)
        // dim+1 because index 0 in FlowData is used for storing pressure
        elVel[dim][actNode] = FlowData[dim+1][connecth[actNode]-1]; 
  }

  void LinearFlowNoiseInt::CalcElemVec4QuadwithVel(const Matrix<Double>& ptCoord,
                                                   const Matrix<Double> & NodalVel,
                                                   Vector<Double> & Result,
                                                   const Elem* elem)
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4Quad" << std::endl;
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
    Vector<double> intWeights = ptelem->GetIntWeights();

    Double jacDet;
    UInt actInt;
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
	ptelem->GetShFncAtIp(Sf, actInt, elem);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet, elem);

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
        if(n == 8)
        {
          helpVec[0] += VelAtIP[2] *  VelDerAtIP[0][2]
                        + VelAtIP[0] *  VelDerAtIP[2][2]; 
        }
        

	helpVec[1] = 2.0 * VelAtIP[1] * VelDerAtIP[1][1] 
	  + VelAtIP[0] *  VelDerAtIP[1][0] 
	  + VelAtIP[1] *  VelDerAtIP[0][0]; 
        if(n == 8)
        {
            helpVec[1] += VelAtIP[1] *  VelDerAtIP[2][2]
                          + VelAtIP[2] *  VelDerAtIP[1][2]; 
        }
        

        if(n == 8)
        {
            helpVec[2] = 2.0 * VelAtIP[2] * VelDerAtIP[2][2] 
                         + VelAtIP[0] *  VelDerAtIP[2][0] 
                         + VelAtIP[2] *  VelDerAtIP[0][0]
                         + VelAtIP[1] *  VelDerAtIP[2][1]
                         + VelAtIP[2] *  VelDerAtIP[1][1]; 
        }

	  
	helpVec *= density;
        
        // Multiplication with the derivatives of the shape functions
	partResult  = xiDx * helpVec;
        partResult *= jacDet * intWeights[actInt-1];
        //std::cout<<"JacDet= "<<jacDet<<std::endl;
        Result     += partResult;
      }

  } // end of method

//===================================================================================
//=========================== Combustion Noise ======================================
//===================================================================================


void LinearFlowNoiseInt::CalcElemVec4CombustionTij(const Matrix<Double>& ptCoord,
                                                   const Matrix<Double>& NodalVel,
                                                   const Vector<Double>& NodalRho,
                                                   Vector<Double>& Result,
                                                   const Elem* elem) {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4CombustionTij" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
    // T_ij = rho v_i v_j
    // Source term =  integral(grad[Sf].div[T])

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetSizeRow();

    Matrix<Double> xiDx;      
    Matrix<Double> xiDxTransp;      
    Vector<Double>  Sf;

    Double RhoAtIP;
    Vector<Double> RhoDerAtIP;
    Vector<Double> VelAtIP;
    Matrix<Double> VelDerAtIP;
    VelDerAtIP.Resize(dim);  

    Vector<double> intWeights = ptelem->GetIntWeights();

    Double jacDet;
    UInt actInt;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
    Vector<Double> helpVec;
    helpVec.Resize(dim);

    
    // Loop over all integration points 
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt, elem);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet, elem);

	// velocity at integration point: (vx  vy)^T  (2x1)
	VelAtIP = NodalVel * Sf;

	//density at integration point
	RhoAtIP = NodalRho * Sf;

	//first derivative of velocity at integration point: (2x2)
	// xiDx: 
	//  vx,x   vx,y
        //  vy,x   vy,y
	//
        VelDerAtIP = NodalVel * xiDx;

	//first order derivative of density at integration point
	xiDx.Transpose(xiDxTransp);
	RhoDerAtIP = xiDxTransp * NodalRho;

	ComputeDerivOfTij(RhoAtIP,RhoDerAtIP,VelAtIP, VelDerAtIP,dim, helpVec);
        
        // Multiplication with the derivatives of the shape functions
	partResult  = xiDx * helpVec;
        partResult *= jacDet * intWeights[actInt-1];
        Result     += partResult;
      }

  } // end of method


void LinearFlowNoiseInt::CalcElemVec4CombustionVec(const Matrix<Double>& ptCoord,
						   const Matrix<Double>& NodalVec,
                                                   Vector<Double>& Result,
                                                   const Elem* elem) {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4CombustionVec" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times  vec
    // vec = partial / partial t ( (1-rho) * velocity )
    // Source term =  integral(grad[Sf].vec)

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetSizeRow();

    Matrix<Double> xiDx;      
    Vector<Double>  Sf;
    Vector<Double> VecAtIP;

    Vector<double> intWeights = ptelem->GetIntWeights();

    Double jacDet;
    UInt actInt;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
    Vector<Double> helpVec;
    helpVec.Resize(dim);

    
    // Loop over all integration points 
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt, elem);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet, elem);

	// velocity at integration point: (vx  vy)^T  (2x1)
	VecAtIP = NodalVec * Sf;

        // Multiplication with the derivatives of the shape functions
	partResult  = xiDx * VecAtIP;
        partResult *= jacDet * intWeights[actInt-1];
        Result     += partResult;
      }

  } // end of method


void LinearFlowNoiseInt::CalcElemVec4CombustionScalar(const Matrix<Double>& ptCoord,
						      const Vector<Double>& NodalScalar,
						      Vector<Double>& Result,
                                                      const Elem* elem) {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4CombustionScalar" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the shape function times  the scalar value
    // scalar = partial / partial t ( 1/rho * partial_C ( rho*partial omega / partial t ) )
    // Source term =  integral(Sf * scalar)

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();

    Vector<Double>  Sf;
    Double ScalarAtIP;
    Vector<Double> intWeights = ptelem->GetIntWeights();

    Double jacDet;
    UInt actInt;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
    // Loop over all integration points 
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt, elem);
	jacDet = ptelem->CalcJacobianDetAtIp(actInt, ptCoord, elem);

	// velocity at integration point: (vx  vy)^T  (2x1)
	ScalarAtIP = NodalScalar * Sf;

        // Multiplication with the derivatives of the shape functions
        Sf *= ScalarAtIP * jacDet * intWeights[actInt-1];
        Result     += Sf;
      }

} // end of method


void LinearFlowNoiseInt::CalcElemVec4CombustionTijOnSurface(const Matrix<Double>& ptCoord,
							    const Matrix<Double>& NodalVel,
							    const Vector<Double>& NodalRho,
							    Vector<Double>& Result,
                                                            const Elem* elem) {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4CombustionTij" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
    // T_ij = rho v_i v_j
    // Source term =  integral(grad[Sf].div[T])

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetSizeRow();

    Matrix<Double> xiDx;      
    Matrix<Double> xiDxTransp;      
    Vector<Double>  Sf;

    Double RhoAtIP;
    Vector<Double> RhoDerAtIP;
    Vector<Double> VelAtIP;
    Matrix<Double> VelDerAtIP;
    VelDerAtIP.Resize(dim);  

    Vector<double> intWeights = ptelem->GetIntWeights();

    Double jacDet;
    UInt actInt;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;
    
    Vector<Double> helpVec;
    helpVec.Resize(dim);

    Vector<Double> nVec(dim);
    ComputeNormalVec( ptCoord, nVec);
    
    // Loop over all integration points 
    for (actInt=1; actInt<=l; actInt++)
      {
	ptelem->GetShFncAtIp(Sf, actInt, elem);
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet, elem);

	// velocity at integration point: (vx  vy)^T  (2x1)
	VelAtIP = NodalVel * Sf;

	//density at integration point
	RhoAtIP = NodalRho * Sf;

	//first derivative of velocity at integration point: (2x2)
	// xiDx: 
	//  vx,x   vx,y
        //  vy,x   vy,y
	//
        VelDerAtIP = NodalVel * xiDx;

	//first order derivative of density at integration point
	xiDx.Transpose(xiDxTransp);
	RhoDerAtIP = xiDxTransp * NodalRho;

	ComputeDerivOfTij(RhoAtIP,RhoDerAtIP,VelAtIP, VelDerAtIP,dim,helpVec);
        
	// scalar product with normal vector
	Double val = nVec * helpVec;

        // Multiplication with the derivatives of the shape functions
        Sf *= val * jacDet * intWeights[actInt-1];
        Result     += Sf;
      }

} // end of method


void LinearFlowNoiseInt::CalcElemVec4CombustionVectorOnSurface(const Matrix<Double>& ptCoord,
							       const Matrix<Double>& NodalVel,
							       Vector<Double>& Result,
                                                               const Elem* elem) {
#ifdef TRACE
  (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4CombustionVectorOnSurface" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the shape function times scalar product of the vector with he normal vector
    // vector = \partial / \partial t ( (1-rho) v )

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetSizeRow();

    Matrix<Double> xiDx;      
    Vector<Double>  Sf;
    Vector<Double> VelAtIP;
    Vector<double> intWeights = ptelem->GetIntWeights();

    Double jacDet;
    UInt actInt;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;
    
    Vector<Double> nVec(dim);
    ComputeNormalVec( ptCoord, nVec);
    
    // Loop over all integration points 
    for (actInt=1; actInt<=l; actInt++)
      {
	jacDet = ptelem->CalcJacobianDetAtIp(actInt, ptCoord, elem);
	ptelem->GetShFncAtIp(Sf, actInt, elem);

	// velocity at integration point: (vx  vy)^T  (2x1)
	VelAtIP = NodalVel * Sf;

	// scalar product with normal vector
	Double val = nVec * VelAtIP;

	Sf    *= val * jacDet * intWeights[actInt-1];
        Result += Sf; 
      }

} // end of method


void LinearFlowNoiseInt::ComputeDerivOfTij( const Double RhoAtIP,
					    const Vector<Double>& RhoDerAtIP,
					    const Vector<Double>& VelAtIP, 
					    const Matrix<Double>& VelDerAtIP,
					    const Integer dim,
					    Vector<Double>& helpVec) {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::ComputeDerivOfTij" << std::endl;
#endif

    // This functions computes derivative of T_ij
    // T_ij = rho v_i v_j

    helpVec[0] = VelAtIP[0]*VelAtIP[0]*RhoDerAtIP[0] 
      + 2.0 * RhoAtIP * VelAtIP[0] * VelDerAtIP[0][0] 
      + VelAtIP[0] * VelAtIP[1] * RhoDerAtIP[1] 
      + RhoAtIP * VelAtIP[1] *  VelDerAtIP[0][1]
      + RhoAtIP * VelAtIP[0] * VelDerAtIP[1][1] ;
    
    if ( dim == 3 ) {
      helpVec[0] += VelAtIP[0] * VelAtIP[2] * RhoDerAtIP[2]
	+ RhoAtIP * VelAtIP[0] * VelDerAtIP[2][2]
	+ RhoAtIP * VelAtIP[2] * VelDerAtIP[0][2];
    }
    
    helpVec[1] = RhoAtIP * VelAtIP[0] * VelDerAtIP[1][0] 
      + RhoAtIP * VelAtIP[1] *  VelDerAtIP[0][0] 
      + VelAtIP[0] * VelAtIP[1] *  RhoDerAtIP[0]
      + VelAtIP[1] * VelAtIP[1] * RhoDerAtIP[1]
      + 2 * RhoAtIP * VelAtIP[1] * VelDerAtIP[1][1];
    
    if( dim == 3 ) {
      helpVec[1] += VelAtIP[1] *  VelAtIP[2] * RhoDerAtIP[2] 
	+ RhoAtIP * VelAtIP[1] * VelDerAtIP[2][2]
	+ RhoAtIP * VelAtIP[2] * VelDerAtIP[1][2];
      
      helpVec[2] = RhoAtIP * VelAtIP[0] * VelDerAtIP[2][0] 
	+ RhoAtIP * VelAtIP[2] *  VelDerAtIP[0][0] 
	+ VelAtIP[0] * VelAtIP[2] *  RhoDerAtIP[0]
	+ RhoAtIP * VelAtIP[1] * VelDerAtIP[2][1]
	+ RhoAtIP * VelAtIP[2] * VelDerAtIP[1][1]
	+ VelAtIP[1] * VelAtIP[2] * RhoDerAtIP[1]
	+ VelAtIP[2] * VelAtIP[2] * RhoDerAtIP[2]
	+ 2 * RhoAtIP * VelAtIP[2] * VelDerAtIP[2][2];
    }
}


void LinearFlowNoiseInt::ComputeNormalVec( const Matrix<Double>& ptCoord,
					    Vector<Double>& nVec) {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::ComputeNormalVec" << std::endl;
#endif

    // compute the normal vector
    UInt dim = nVec.GetSize();
    Vector<Double> vec1(dim);
    Vector<Double> vec2(dim);

    if ( dim == 3) {
      UInt idx;
      if ( ptCoord.GetSizeCol() == 4) {
	idx = 3;
      }
      else if ( ptCoord.GetSizeCol() == 3) {
	idx = 2;
      }
      else {
	Error("Surface Element is no quadrilateral or triangle",__FILE__,__LINE__);
      }

      for ( UInt i=0; i<dim; i++ ) {
	vec1[i] = ptCoord[i][1]   - ptCoord[i][0];
	vec2[i] = ptCoord[i][idx] - ptCoord[i][0];
      }

      // cross product: vec1 x vec 2
      nVec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
      nVec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
      nVec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    }
    else {
      nVec[0] =  ptCoord[1][1] - ptCoord[1][0];
      nVec[1] = -ptCoord[0][1] + ptCoord[0][0];
    }

    Double invLength = 1.0 / nVec.NormL2();
    nVec *= invLength;
}



  // ===================================================================
  // nonlinear RHS for nonlinear acoustics
  // ===================================================================

  nLinKuznetsovRHSInt::nLinKuznetsovRHSInt( bool isaxi )
    : LinearForm()
  {
    name_ = "nLinKuznetsovRHSInt";
    isaxi_ = isaxi;

  }


  nLinKuznetsovRHSInt::~nLinKuznetsovRHSInt()
  {
  }


  void nLinKuznetsovRHSInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent ) 
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt spaceDim = ptelem->GetDim();  

    elemVec.Resize(numFncs);
    elemVec.Init(0.0);

    solGradAtIp_.Resize(numFncs);
    solDeriv1GradAtIp_.Resize(numFncs);

//     // calculate real factors
//     Double fac1,fac2;
//     mParser_->SetExpr( mHandle_, factorN1_);
//     fac1 = mParser_->Eval( mHandle_ );
//     mParser_->SetExpr( mHandle_, factorN2_);
//     fac2 = mParser_->Eval( mHandle_ );


    Double jacDet, totalfactor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = 0;
      ptelem->GetShFncAtIp(ShpFncAtIp_,actIntPt, ent.GetElem() );
      ptelem->GetGlobDerivShFncAtIp(xiDx_, actIntPt, ptCoord_, 
                                    jacDet, ent.GetElem() );
        
      if (isaxi_) {
        CoordAtIp_ = ptCoord_ * ShpFncAtIp_;
        for (UInt i=0; i<numFncs; i++)
          xiDx_[i][0] += ShpFncAtIp_[i] / CoordAtIp_[0];
          
        jacDet *= 2 * PI * CoordAtIp_[0];
      }
        

      //compute gradient of solution + 1st derivative at integration point
      solGradAtIp_.Init();
      solDeriv1GradAtIp_.Init();
      for ( UInt i=0; i<spaceDim; i++ ) {
        for ( UInt j=0; j<numFncs; j++ ) {
          solGradAtIp_[i] += xiDx_[j][i] * sol_[j];
          solDeriv1GradAtIp_[i] += xiDx_[j][i] * solderiv1_[j];
        }
      }
//       xiDx.Transpose(xiDxTransp);
//       solGradAtIp       = xiDxTransp * sol_;
//       solDeriv1GradAtIp = xiDxTransp * solderiv1_;


      //get 1st and 2nd derivative of solution at integration point
      solDeriv1AtIp_ = solderiv1_ * ShpFncAtIp_;
      solDeriv2AtIp_ = solderiv2_ * ShpFncAtIp_;
        
      totalfactor=0;
      for (UInt j=0; j<xiDx_.GetSizeCol(); j++)
        totalfactor += solGradAtIp_[j] * solDeriv1GradAtIp_[j];
      totalfactor *= factorN2_;
        
      totalfactor += factorN1_ * solDeriv1AtIp_ * solDeriv2AtIp_;
        
      totalfactor *= jacDet;
      for (UInt i=0; i< numFncs; i++)
        elemVec[i] += ShpFncAtIp_[i] * totalfactor;
        
    }
    //  std::cerr << "RHS in linearForm:\n" << elemVec << std::endl;
  }

  nLinWesterveltRHSInt::nLinWesterveltRHSInt( bool isaxi )
    : LinearForm()
  {
    name_ = "nLinWesterveltRHSInt";
    isaxi_ = isaxi;
  }


  nLinWesterveltRHSInt::~nLinWesterveltRHSInt()
  {
  }


  void nLinWesterveltRHSInt::CalcElemVector( Vector<Double> & elemVec,
                                             EntityIterator& ent )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
   
    elemVec.Resize(numFncs);
    elemVec.Init(0.0);

    Double jacDet, totalfactor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem() );
      ptelem->GetShFncAtIp(ShpFncAtIp_,actIntPt, ent.GetElem() );
        
      if (isaxi_) {
        CoordAtIp_ = ptCoord_ * ShpFncAtIp_;
        jacDet *= 2 * PI * CoordAtIp_[0];
      }
        
      //get solution and derivatives at integration point
      solAtIp_ = sol_ * ShpFncAtIp_;
      solDeriv1AtIp_ = solderiv1_ * ShpFncAtIp_;
      solDeriv2AtIp_ = solderiv2_ * ShpFncAtIp_;
        
      totalfactor = jacDet * 2.0 * factor_ * (solDeriv1AtIp_ * solDeriv1AtIp_ +
                                              solAtIp_ * solDeriv2AtIp_);
        
      for (UInt i=0; i< numFncs; i++)
        elemVec[i] += ShpFncAtIp_[i] * totalfactor;
    }
  }


  // ====================================================================
  // electric polarization
  // ====================================================================

  ElecPolarizationInt::ElecPolarizationInt(bool isaxi, bool coordUpdate )
    : LinearForm()
  {
    name_ = "ElecPolarizationInt";
    isaxi_ = isaxi;
    coordUpdate_ = coordUpdate;
  }


  ElecPolarizationInt::~ElecPolarizationInt()
  {
  }


  void ElecPolarizationInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent )
  {
    
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
   
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> ShpFncAtIp, CoordAtIP;
    Matrix<Double> xiDx;

    elemVec.Resize(numFncs);
    elemVec.Init(0);

    Double factor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = 0;
      ptelem->GetGlobDerivShFncAtIp(xiDx, actIntPt, ptCoord_, 
                                    jacDet, ent.GetElem() );

      if (isaxi_) {
        ptelem->GetShFncAtIp(ShpFncAtIp, actIntPt, ent.GetElem() );
        CoordAtIP = ptCoord_ * ShpFncAtIp;
        for (UInt i=0; i<numFncs; i++)
          xiDx[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
      
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      factor = intWeights[actIntPt-1] * jacDet;
      for (UInt i=0; i<numFncs; i++) {
        for (UInt j=0; j<D_.GetSize(); j++) {
          elemVec[i] += xiDx[i][j] * D_[j];
        }
      }
    }

  }


  // ==================================================================
  // piezoelectric polarization: 
  // ==================================================================

  // base class
  PiezoPolarizationRhsInt::PiezoPolarizationRhsInt( BaseMaterial* matDataElec,
						    SubTensorType type ) 
    : LinearForm() {

    matDataElec_   = matDataElec;
    subTensorType_ = type;

    isaxi_ = false;
    if ( subTensorType_ == AXI )
      isaxi_ = true;

    // get direction of polarization
    std::string str;
    matDataElec_->GetScalar(str, P_DIRECTION);
    Directions dir;
    String2Enum(str,dir);
    dirP_ = dir;

    piezoBilinearForm_ = new linPiezoCoupling(matDataElec, type); 


    //get basis vector for irreversibel strain
    ComputeNormalizedSirr(subTensorType_, dirP_, baseVecSirr_);

  }

  PiezoPolarizationRhsInt::~PiezoPolarizationRhsInt() {

    delete EfieldOp_;
    delete piezoBilinearForm_;
    delete matDataElec_;
  }

  void PiezoPolarizationRhsInt::ComputeNormalizedSirr(SubTensorType type, 
						      UInt dirP,
						      Vector<Double>& baseVecSirr) {


    //get basis vector for irreversibel strain
    if ( type == FULL ) {
      baseVecSirr.Resize(6);
      baseVecSirr.Init();
      for ( UInt i=0; i<3; i++ )
	baseVecSirr[i] = -0.5;
      baseVecSirr[ dirP ] = 1.0;
    }
    else if ( type ==  AXI ) {
      baseVecSirr.Resize(4);
      baseVecSirr.Init();
      for ( UInt i=0; i<4; i++ )
	baseVecSirr[i] = -0.5;
      baseVecSirr[2] = 0.0;  // no shear strain
      baseVecSirr[ dirP ] = 1.0;
    }
    else {
      baseVecSirr.Resize(3);
      baseVecSirr.Init();
      for ( UInt i=0; i<1; i++ )
	baseVecSirr[i] = -0.5;
      baseVecSirr[ dirP ] = 1.0;
    }
  }


  void PiezoPolarizationRhsInt::Set4NonLinMaterial(Grid* ptGrid, 
						   StdPDE* ptPDE,
						   shared_ptr<EqnMap> eqnMap,
						   shared_ptr<ResultInfo> result) {


    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE, 
                                             eqnMap, *sol_, 
                                             ELEC_POTENTIAL, 
                                             result, isaxi_, 
                                             coordUpdate_);
  }


  // ==================================================================
  // piezoelectric polarization: RHS of electric equation
  // ==================================================================

  PiezoPolarizationElecRhsInt::PiezoPolarizationElecRhsInt( BaseMaterial* matDataElec,
							    BaseMaterial* matDataPiezo,
							    BaseMaterial* matDataMech,
							    SubTensorType type) 
    : PiezoPolarizationRhsInt( matDataElec, type ) {


    name_  = "PiezoPolarisationElecRhsInt";

    matDataPiezo_ = matDataPiezo;
    matDataMech_  = matDataMech;
  }


  PiezoPolarizationElecRhsInt::~PiezoPolarizationElecRhsInt() {

    delete matDataPiezo_;
    delete matDataMech_; 
  }


  void PiezoPolarizationElecRhsInt::CalcElemVector( Vector<Double> & elemVec,
						    EntityIterator& ent ) 
    
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrNodes  = ptelem->GetNumNodes();      
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
   
    Double jacDet;
    Vector<Double> partElemVec, Sirr, partSirr, elecD;
    Vector<Double> ShpFncAtIp;
    Matrix<Double> bMat, bMatTrans;

    elemVec.Resize(numFncs);
    elemVec.Init(0);  

    //compute polarization
    Vector<Double> LCoord, Efield;

    // calc E-Field:
    //compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);

    EfieldOp_->CalcElemGradField( Efield, ent, LCoord, 1);
    UInt nrEl = ent.GetElem()->elemNum;
    Double actP = matDataElec_->ComputeScalarHystVal( nrEl, Efield ) ; //Efield[dirP_] );
    Vector<Double> elecP( Efield.GetSize() );
    elecP.Init();
    elecP[ dirP_ ] = actP; 

    // get coupling tensor
    Matrix<Double> ematMatrix;
    matDataPiezo_->GetTensor( ematMatrix, PIEZO_TENSOR, REAL, subTensorType_);

    // compute irreversibel strain
    Vector<Double> coeff;
    matDataMech_->GetVector( coeff, COEFF_STRAIN_IRREVERSIBLE, REAL);

    Double val = 1.5 * ( coeff[0] +
			 coeff[1] * actP +
			 coeff[2] * actP * actP +
			 coeff[3] * actP * actP * actP + 
			 coeff[4] * actP * actP * actP * actP );

    Sirr = baseVecSirr_ * val;

    partSirr = ematMatrix * Sirr;

//     std::cout << "Sirr:\n" << Sirr << std::endl;
//     std::cout << "Pirr:\n" << elecP << std::endl;
//     std::cout << "partD_Sirr:\n" << partSirr << std::endl;

    // P - [e] * Sirr
    elecD = elecP - partSirr;

    piezoBilinearForm_->ExtractElemInfo( ent );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      // Setup the B matrix for current integration point
      piezoBilinearForm_->calcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,ent.GetElem() );

      //shape functions at IP  
      ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent.GetElem() );

      if ( isaxi_ ) {
        Double aux = 0.0;
        
        for ( UInt i = 0; i < nrNodes; i++ ) {
          aux += ptCoord_[0][i] * ShpFncAtIp[i];
        } 
        jacDet *= 2.0 * PI * aux;
      }
      jacDet *= intWeights[actIntPt-1];

  
      bMat.Transpose(bMatTrans);
      partElemVec = bMatTrans * elecD;

      for ( UInt i=0; i<partElemVec.GetSize(); i++)
	partElemVec[i] *=  ShpFncAtIp[i] * jacDet;

      elemVec += partElemVec;
    }
    
    //    std::cout << "elemPvec:\n" << elemVec << std::endl;
  }


  // ==================================================================
  // piezoelectric polarization: RHS of mechanic equation
  // ==================================================================

  PiezoPolarizationMechRhsInt::PiezoPolarizationMechRhsInt( BaseMaterial* matDataElec,
							    BaseMaterial* matDataMech,
							    SubTensorType type) 
    : PiezoPolarizationRhsInt( matDataElec, type ) {


    name_  = "PiezoPolarisationMechRhsInt";
    matDataMech_  = matDataMech;

  }


  PiezoPolarizationMechRhsInt::~PiezoPolarizationMechRhsInt()
  {

    delete matDataMech_; 
  }


  void PiezoPolarizationMechRhsInt::CalcElemVector( Vector<Double> & elemVec,
						    EntityIterator& ent ) 
    
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrNodes  = ptelem->GetNumNodes();      
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    UInt dim = 2;
    if ( subTensorType_ == FULL ) 
      dim = 3;

    Double jacDet;
    Vector<Double> partElemVec, strainIrr, stressIrr;
    Vector<Double> ShpFncAtIp;
    Matrix<Double> bMat;

    elemVec.Resize( numFncs*dim );
    elemVec.Init(0);  

    //    compute polarization
    Vector<Double> LCoord, Efield;

    // compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);

    EfieldOp_->CalcElemGradField( Efield, ent, LCoord, 1);
    UInt nrEl = ent.GetElem()->elemNum;
    Double actP = matDataElec_->ComputeScalarHystVal( nrEl, Efield ) ; //[dirP_] );

    // get coupling tensor
    Matrix<Double> cmatMatrix;
    matDataMech_->GetTensor( cmatMatrix, MECH_STIFFNESS_TENSOR, REAL, subTensorType_);

    // compute irreversibel strain
    Vector<Double> coeff;
    matDataMech_->GetVector( coeff, COEFF_STRAIN_IRREVERSIBLE, REAL);

    Double val = 1.5 * ( coeff[0] +
			 coeff[1] * actP +
			 coeff[2] * actP * actP +
			 coeff[3] * actP * actP * actP + 
			 coeff[4] * actP * actP * actP * actP );

    strainIrr = baseVecSirr_ * val;
    stressIrr = cmatMatrix * strainIrr;

    piezoBilinearForm_->ExtractElemInfo( ent );

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      // Setup the B matrix for current integration point
      piezoBilinearForm_->calcAMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,ent.GetElem() );

      // shape functions at IP  
      ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent.GetElem() );

      if ( isaxi_ ) {
        Double aux = 0.0;
        for ( UInt i = 0; i < nrNodes; i++ ) {
          aux += ptCoord_[0][i] * ShpFncAtIp[i];
        }
        jacDet *= 2.0 * PI * aux;
      }

      jacDet *= intWeights[actIntPt-1];

      partElemVec = bMat * stressIrr;
      UInt idx = 0;
      for ( UInt i=0; i<numFncs; i++) {
	for ( UInt j=0; j<dim; j++ ) {
	  partElemVec[idx] *=  ShpFncAtIp[i] * jacDet;
	  idx += 1;
	}
      }
      elemVec -= partElemVec;
    }
    
//     std::cout << "elemSvec:\n" << elemVec << std::endl;
  }



  // =========================================================================
  // mechanic rhs volume integrator
  // =========================================================================
  
  
  MechVolForceInt::MechVolForceInt(UInt numDof, 
                                   const std::string& phase,
                                   bool isaxi) {


    name_ = "MechVolForceInt";
    isaxi_ = isaxi;
    numDofs_ = numDof;
    phase_ = phase;

  }
    
  MechVolForceInt::~MechVolForceInt() {


  }
    
  void MechVolForceInt::SetVolForceVector(StdVector<std::string> & volForce, 
                                          const CoordSystem * coordSys,
                                          bool isUnit, Double volume ) {


    locForce_ = volForce;
    coordSys_ = coordSys;
    isUnitValue_ = isUnit;
    volume_ = volume;
    
  }
    
  void MechVolForceInt::CalcElemVector( Vector<Double> & elemVec,
                                        EntityIterator& ent ) {


    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    // get global coordinate system and math parser
    MathParser * parser = domain->GetMathParser();

    // First, map force to global coordinate system
    Vector<Double> globMidPoint, locMidPoint;

    ptelem->GetCoordMidPoint(locMidPoint  );
    ptelem->Local2GlobalCoord(globMidPoint, locMidPoint, 
                              ptCoord_, ent.GetElem() );

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
    Vector<Double> globVec;
    coordSys_->Local2GlobalVector(globVec, locLoadVec,  globMidPoint);

    // Calculate vector
    CalcPartVector( elemVec, globVec, ent );
  }

  void MechVolForceInt::CalcElemVector( Vector<Complex> & elemVec,
                                        EntityIterator& ent ) {


    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    // get global coordinate system and math parser
    MathParser * parser = domain->GetMathParser();

    // First, map force to global coordinate system
    Vector<Double> globMidPoint, locMidPoint;

    ptelem->GetCoordMidPoint(locMidPoint);
    ptelem->Local2GlobalCoord( globMidPoint, locMidPoint, 
                               ptCoord_, ent.GetElem() );

    // Update variables for mathParser
    parser->SetCoordinates( mHandle_, *coordSys_, globMidPoint );

    // Now evaluate each entry 
    Vector<Complex> locLoadVec( numDofs_ );
    Double amplitude, phase;

    // -- phase --
    parser->SetExpr( mHandle_, phase_ );
    phase = parser->Eval( mHandle_ );

    for ( UInt i = 0; i < locForce_.GetSize(); i++ ) {
      parser->SetExpr( mHandle_, locForce_[i] );
      amplitude = parser->Eval( mHandle_ );
      locLoadVec[i] = Complex( amplitude * cos(phase/180*PI), 
                               amplitude * sin(phase/180*PI) );
    }

    // If load is not unit load, divide by volume
    if ( isUnitValue_ == false ) {
      locLoadVec /= volume_;
    }
    
    // Map local load vector to global one
    Vector<Complex> globVec;
    coordSys_->Local2GlobalVector(globVec, locLoadVec,  globMidPoint);

    // Calculate vector
    CalcPartVector( elemVec, globVec, ent );
  }

  template<class TYPE>
  void MechVolForceInt::CalcPartVector( Vector<TYPE>& elemVec, 
                                        Vector<TYPE>& loadVec,
                                        EntityIterator& ent ) {
    

    /// ----- part of interior function ---
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Vector<Double> shapeFnc, CoordAtIP;

    // Then, calculate element vector
    elemVec.Resize(numFncs * numDofs_);
    elemVec.Init(0.0);
    Double factor;

    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem() );
      ptelem->GetShFncAtIp(shapeFnc,actIntPt, ent.GetElem() );

      if (isaxi_) {
        CoordAtIP = ptCoord_ * shapeFnc;
        jacDet *= 2 * PI * CoordAtIP[0];
      }
          
      factor = intWeights[actIntPt-1] * jacDet;
      for (UInt iFnc = 0; iFnc < numFncs; iFnc++) {
        for (UInt iDof = 0; iDof < numDofs_; iDof++) {
          elemVec[iFnc*numDofs_ + iDof] += 
            factor * shapeFnc[iFnc] * loadVec[iDof]; 
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
  }

  void AddStressRHSInt::CalcElemVector( Vector<Double> & elemVec,
					 EntityIterator& ent )
  {

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
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
      
        Double jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_,
                                                    ent.GetElem());
      
        if (isaxi_) {
	  Vector<Double> shpFncAtIp;
	  Vector<Double> coordAtIP;
	  ptelem->GetShFncAtIp(shpFncAtIp, actIntPt, ent.GetElem() );
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
