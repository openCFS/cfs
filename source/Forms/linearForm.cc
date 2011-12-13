// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <math.h>
#include <complex>
#include <iostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "Forms/curlCurlEdgeInt.hh"
#include "Forms/curlCurlNodeInt.hh"
#include "Forms/linElastInt.hh"
#include "Forms/linGradBDBInt.hh"
#include "Forms/nLinElastInt.hh"
#include "Forms/nLincurlCurlEdgeInt.hh"
#include "General/exception.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "Materials/baseMaterial.hh"
#include "Utils/biotSavart.hh"
#include "Utils/coordSystem.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/tools.hh"
#include "linMagStrictInt.hh"
#include "linearForm.hh"
#include "math.h"
#include "nLincurlCurlNodeInt.hh"

namespace CoupledField {

class ApproxData;

DECLARE_LOG(linForm)
DEFINE_LOG(linForm, "linForm")

  LinearForm::LinearForm( BaseMaterial* matData) : BaseForm( matData )
  {}

  LinearForm::~LinearForm()
  {}

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

  LinearEdgeSrcInt::LinearEdgeSrcInt( UInt numDof, 
                                       const std::string& phase,
                                       bool isaxi) 
     :  VolForceInt( numDof, phase, isaxi )
   {
     name_ = "LinearEdgeSrcInt";
   
   }

   LinearEdgeSrcInt::~LinearEdgeSrcInt()
   {
   }

   void LinearEdgeSrcInt::CalcElemVector( Vector<Double> & elemVec,
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

   void LinearEdgeSrcInt::CalcElemVector( Vector<Complex> & elemVec,
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
   void LinearEdgeSrcInt::CalcPartVector( Vector<TYPE>& elemVec, 
                                          Vector<TYPE>& loadVec,
                                          EntityIterator& ent )
   {
   
     // Extract pointer to reference element and get coordinates
     ExtractElemInfo( ent );

     const UInt nrIntPts = ptelem->GetNumIntPoints();
     const UInt nrEdges  = ptelem->GetNumEdges();
     const Vector<Double> & intWeights = ptelem->GetIntWeights();  
   
     Double jacDet;  
     
     // derivation of shape functions after global coordinates 
     Matrix<Double> shapeEdge;
     Vector<TYPE> partElemVec;
     partElemVec.Resize(nrEdges);
     partElemVec.Init();
   
   
     // set vector to desired size and set all elements to zero    
     elemVec.Resize(nrEdges); 
     elemVec.Init();

     for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {    
       ptelem->CalcEdgeShapeFncAtIp(shapeEdge, actIntPt, ptCoord_,
                                    ent.GetElem());

       partElemVec = shapeEdge * loadVec;
       
       jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem());

       for(UInt i=0; i<partElemVec.GetSize(); i++)
         elemVec[i] += partElemVec[i] * intWeights[actIntPt-1] * jacDet; 
     }
   
/*
     (*debug) << "CalcElemVector:  "  << std::endl
              << partElemVec << std::endl
              << "\n jacDet " << jacDet << std::endl;
 */
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
    curlOp_->SetAnsatzFct( ansatzFct1_ );

    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    Vector<Double> CoordAtIP;
    Matrix<Double> bMatCurl, bMatDiv, bMatCurlTrans; 
    Vector<Double> helpVec;
    Double jacDet;
    // Double factor; // TODO: Unused variable factor

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
        EXCEPTION("CurlCurlNode3DInt::CalcElementMatrix: Encountered "
                  << "negative Jacobian determinant!");
      }

      // factor = intWeights[actIntPt-1] * jacDet * reluctivity_;
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
    ptMaterial->GetScalar( startmatVal_, MAG_RELUCTIVITY,Global::REAL);

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
    nlinFnc_ = ptMaterial->GetNonlinFnc(MAG_PERMEABILITY);

    BaseForm * curlcurl2D;
    if (nlinFnc_== NULL) 
      //define the linear element matrix
      curlcurl2D = new CurlCurlNode2DInt( ptMaterial, isaxi_, coordUpdate_ );
    else {
      //define the nonlinear element matrix
      curlcurl2D = new nLinCurlCurlNode2DInt( ptMaterial, isaxi_, coordUpdate_ );
      //important to set method to FixPoint, since we compute the RHS!!
      curlcurl2D->SetNonLinMethod(FIXEDPOINT);
        
      //set the solution class to the operator
      curlcurl2D->SetSolution( * sol_ );
    }
    curlcurl2D->SetAnsatzFct( ansatzFct1_ );

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
    ptMaterial->GetScalar( startmatVal_, MAG_RELUCTIVITY,Global::REAL);

  }
  
  nLinMagNode3D_linFormInt::~nLinMagNode3D_linFormInt()
  {
  }

  void nLinMagNode3D_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                                 EntityIterator& ent )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    
    ptelem->SetAnsatzFct( ansatzFct1_ );

    // get pointer to nonlinear BH curve approximation
    nlinFnc_ = ptMaterial->GetNonlinFnc(MAG_PERMEABILITY);

    BaseForm * curlcurl3D;
    if (nlinFnc_== NULL) 
      //define the linear element matrix
      curlcurl3D = new CurlCurlNode3DInt( ptMaterial, coordUpdate_);
    else {
      //define the nonlinear element matrix
      curlcurl3D = new nLinCurlCurlNode3DInt( ptMaterial, coordUpdate_ );
      //important to set method to FixPoint, since we compute the RHS!!
      curlcurl3D->SetNonLinMethod(FIXEDPOINT);
        
      //set the solution class to the operator
      curlcurl3D->SetSolution( *sol_ );
    }
    curlcurl3D->SetAnsatzFct( ansatzFct1_ );

    // Get element solution
    Vector<Double> magPot;
    sol_->GetElemSolution( magPot, ent  );

    Matrix<Double> elemmat;
    curlcurl3D->CalcElementMatrix(elemmat, ent, ent);

    elemVec = -(elemmat * magPot);

    delete curlcurl3D;
  }

  

  nLinMagEdge_linFormInt::nLinMagEdge_linFormInt( BaseMaterial*matData,
                                                  bool coordUpdate)
    : LinearForm( matData )
  {
    name_ = "nLinMagEdge_linFormInt";
    isSolDependent_ = true;

    isaxi_       = false;
    coordUpdate_ = coordUpdate;

  }
  
  nLinMagEdge_linFormInt::~nLinMagEdge_linFormInt()
  {
  }

  void nLinMagEdge_linFormInt::CalcElemVector( Vector<Double> & elemVec,
                                               EntityIterator& ent )
  {
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );
    ptelem->SetAnsatzFct( ansatzFct1_ );

    // get pointer to nonlinear BH curve approximation
    ApproxData* nlinFnc_ = ptMaterial->GetNonlinFnc(MAG_PERMEABILITY);

    BaseForm * curlcurl3D;
    if ( nlinFnc_ == NULL ) 
      //define the linear element matrix
      curlcurl3D = new CurlCurlEdgeInt( ptMaterial, coordUpdate_);
    else {
      //define the nonlinear element matrix
      curlcurl3D = new nLinCurlCurlEdgeInt( ptMaterial, coordUpdate_ );
      //important to set method to FixPoint, since we compute the RHS!!
      curlcurl3D->SetNonLinMethod(FIXEDPOINT);
        
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
      EXCEPTION("Wrong space dimension of elements! ");
    }  
    stressBiformInt->SetAnsatzFct( ansatzFct1_ );
    
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
    EXCEPTION( "PreStressLinFormInt::CalcElemVector: not working");

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

    //   if (!elemDisp_.GetNumRows() || !elemDisp_.GetNumCols()) 
    //     EXCEPTION("Undefined displacements! ");

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
          EXCEPTION("Incorrect dimension of normal. ");
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
    Integer dimelem = ptCoord.GetNumRows();

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

  void LinearFlowNoiseInt::CalcElemVec_withCFSVel(const Matrix<Double>& ptCoord,
                                                  const Matrix<Double> & NodalVel,
                                                  Vector<Double> & Result,
                                                  Double density)
 {


   // This functions computes the element RHS vector by integrating 
   // the gradient of the shape function times the divergence of the tensor T.
 
   // Source term =  integral(grad[Sf].div[T])

   Integer l = ptelem->GetNumIntPoints(); 
   Integer n = ptelem->GetNumNodes();
   Integer dimelem = ptCoord.GetNumRows();

   Matrix<Double> xiDx;      
   Vector<Double>  Sf;

   Vector<Double>  VelAtIP;
   Matrix<Double> VelDerAtIP;

   Double jacDet;
   Integer actInt;

   Result.Resize(n);
   Result.Init();

   Vector<Double> partResult;
   partResult.Resize(n);
   partResult.Init();
   
   Vector<Double> helpVec;
   helpVec.Resize(dimelem);
   helpVec.Init();

   const Vector<Double> & intWeights = ptelem->GetIntWeights();
   
   // Loop over all integration points 
   for (actInt=1; actInt<=l; actInt++)
     {
	ptelem->GetShFncAtIp(Sf, actInt, NULL  );
	ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord,jacDet, NULL );

	// velocity at integration point: (vx  vy)^T  (2x1)
	VelAtIP = NodalVel * Sf;
      
	//first derivative of velocity at integration point: (2x2)
	//  vx,x   vx,y
       //  vy,x   vy,y
	//
       VelDerAtIP = NodalVel * xiDx;
       
       
       if(ptelem->GetDim() == 2)
         { 
           helpVec[0] = 2.0 * VelAtIP[0] * VelDerAtIP[0][0] 
             + VelAtIP[1] *  VelDerAtIP[0][1] 
             + VelAtIP[0] *  VelDerAtIP[1][1]; 
           
           helpVec[1] = 2.0 * VelAtIP[1] * VelDerAtIP[1][1] 
             + VelAtIP[0] *  VelDerAtIP[1][0] 
             + VelAtIP[1] *  VelDerAtIP[0][0]; 
         }else if(ptelem->GetDim() == 3)
         {
           helpVec[0] = 2.0 * VelAtIP[0] * VelDerAtIP[0][0] 
             + VelAtIP[1] *  VelDerAtIP[0][1] 
             + VelAtIP[0] *  VelDerAtIP[1][1]
             + VelAtIP[0] *  VelDerAtIP[2][2]
             + VelAtIP[2] *  VelDerAtIP[0][2]; 
           
           helpVec[1] = 2.0 * VelAtIP[1] * VelDerAtIP[1][1] 
             + VelAtIP[0] *  VelDerAtIP[1][0] 
             + VelAtIP[1] *  VelDerAtIP[0][0]
             + VelAtIP[1] *  VelDerAtIP[2][2]
             + VelAtIP[2] *  VelDerAtIP[1][2]; 
           
           helpVec[1] = 2.0 * VelAtIP[2] * VelDerAtIP[2][2] 
             + VelAtIP[0] *  VelDerAtIP[2][0] 
             + VelAtIP[2] *  VelDerAtIP[0][0]
             + VelAtIP[2] *  VelDerAtIP[1][1]
             + VelAtIP[1] *  VelDerAtIP[2][1]; 
         }
	  
	helpVec *= density;

       
       //Since we do an integration here, we have to integrate with the 
       //weighting factor here
       //helpVec *= intWeights[actInt-1];

       // Multiplication with the derivatives of the shape functions
	partResult  = xiDx * helpVec;
       partResult *=  intWeights[actInt-1] * jacDet;
       Result += partResult;
     }
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
  
    int dimelem=ptCoord.GetNumRows();

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

  void LinearFlowNoiseInt::CalcElemVec4QuadwithVelCentre(const Matrix<Double>& ptCoord,
                                                   const Matrix<Double> & NodalVel,
                                                   Vector<Double> & Result,
                                                   Vector<Double> & nodalLoadDensity,
                                                   Vector<Double>& divLHTensor,
                                                   const Elem* elem,
                                                   Double density)
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVec4QuadwithVelCentre" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
    // But uses the cell centre velocity for all nodes.
  
    // Source term =  integral(grad[Sf].div[T])

    const Integer numIntPoints = ptelem->GetNumIntPoints(); 
    const Integer numNodes = ptelem->GetNumNodes();
    const UInt dimelem = ptCoord.GetNumRows();

    Matrix<Double> xiDx;      
    Vector<Double> VelAtCentre;
    Matrix<Double> VelDerivAtCentre;
    Vector<Double> intWeights = ptelem->GetIntWeights();
    Vector<Double> divLHTensorTmp;

    Double jacDet;

    Result.Resize(numNodes);
    nodalLoadDensity.Resize(numNodes);
    Result.Init(0.0);
    nodalLoadDensity.Init(0.0);

    Vector<Double> partResult;
    partResult.Resize(numNodes);
    
    Vector<Double> helpVec;
    helpVec.Resize(dimelem);

    Double volume = 0.0;
    divLHTensor.Resize(dimelem);
    divLHTensorTmp.Resize(dimelem);
    divLHTensor.Init(0.0);
    divLHTensorTmp.Init(0.0);
    
    Vector<Double> centreNode(dimelem, 0.0);
    Matrix<Double> derivCentre;
    Vector<Double> ansatzCentre;
    const UInt dof = NodalVel.GetNumCols();

    ptelem->GetGlobDerivShFnc(derivCentre, centreNode, ptCoord, elem, dof);
    ptelem->GetShFnc(ansatzCentre, centreNode, elem, dof);
    VelDerivAtCentre = NodalVel * derivCentre;
    VelAtCentre = NodalVel * ansatzCentre;

    helpVec[0] = 2.0 * VelAtCentre[0] * VelDerivAtCentre[0][0] 
                     + VelAtCentre[1] * VelDerivAtCentre[0][1]
                     + VelAtCentre[0] * VelDerivAtCentre[1][1];
    helpVec[1] = 2.0 * VelAtCentre[1] * VelDerivAtCentre[1][1] 
                     + VelAtCentre[0] * VelDerivAtCentre[1][0] 
                     + VelAtCentre[1] * VelDerivAtCentre[0][0]; 
    helpVec *= density;
    // Loop over all integration points 
    for(Integer actInt=1; actInt <= numIntPoints; actInt++)
    {
      ptelem->GetGlobDerivShFncAtIp(xiDx, actInt, ptCoord, jacDet, elem);
      
      volume += jacDet * intWeights[actInt -1];
      divLHTensorTmp = (helpVec * jacDet * intWeights[actInt -1]);

      partResult  = xiDx * divLHTensorTmp;
      Result     += partResult;

      divLHTensor += divLHTensorTmp;
    } // end integration loop

    nodalLoadDensity = Result / volume;
    divLHTensor /= volume;

  } // end of method

  void LinearFlowNoiseInt::CalcElemVec4QuadwithVel(const Matrix<Double>& ptCoord,
                                                   const Matrix<Double> & NodalVel,
                                                   Vector<Double> & Result,
                                                   Vector<Double> & nodalLoadDensity,
                                                   Vector<Double>& divLHTensor,
                                                   const Elem* elem,
                                                   Double density)
  {
#ifdef TRACE
    (*trace) << "entering LinearFlowNoiseInt::CalcElemVector4Quad" << std::endl;
#endif

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
  
    // Source term =  integral(grad[Sf].div[T])

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dimelem = ptCoord.GetNumRows();

    Matrix<Double> xiDx;      
    Vector<Double> Sf;

    Vector<Double> VelAtIP;
    Matrix<Double> VelDerAtIP;
    Vector<Double> divLHTensorTmp;

    Double jacDet;

    Result.Resize(n);
    nodalLoadDensity.Resize(n);
    Result.Init(0.0);
    nodalLoadDensity.Init(0.0);

    Vector<Double> partResult;
    partResult.Resize(n);
    
    Vector<Double> helpVec;
    helpVec.Resize(dimelem);

    Double volume = 0.0;
    divLHTensorTmp.Resize(dimelem);
    divLHTensorTmp.Init(0.0);
    divLHTensor.Resize(dimelem);
    divLHTensor.Init(0.0);
    
    Vector<double> intWeights = ptelem->GetIntWeights();

    // Loop over all integration points 
    for(Integer actInt=1; actInt <= l; actInt++)
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
                       + VelAtIP[1] * VelDerAtIP[0][1]
                       + VelAtIP[0] * VelDerAtIP[1][1];


      helpVec[1] = 2.0 * VelAtIP[1] * VelDerAtIP[1][1] 
                       + VelAtIP[0] * VelDerAtIP[1][0] 
                       + VelAtIP[1] * VelDerAtIP[0][0]; 

      if (dimelem == 3)
      {
        helpVec[0] +=      VelAtIP[2] * VelDerAtIP[0][2]
                         + VelAtIP[0] * VelDerAtIP[2][2]; 
        helpVec[1] +=      VelAtIP[1] * VelDerAtIP[2][2]
                         + VelAtIP[2] * VelDerAtIP[1][2]; 
        helpVec[2] = 2.0 * VelAtIP[2] * VelDerAtIP[2][2] 
                         + VelAtIP[0] * VelDerAtIP[2][0] 
                         + VelAtIP[2] * VelDerAtIP[0][0]
                         + VelAtIP[1] * VelDerAtIP[2][1]
                         + VelAtIP[2] * VelDerAtIP[1][1]; 
      }

      helpVec *= density;
      divLHTensorTmp = (helpVec * jacDet * intWeights[actInt -1]);
       
      // Multiplication with the derivatives of the shape functions
      partResult  = xiDx * divLHTensorTmp;
      Result     += partResult;
      
      divLHTensor += divLHTensorTmp;
      volume += jacDet * intWeights[actInt-1];
    }

    nodalLoadDensity = Result / volume;
    divLHTensor /= volume;

  } // end of method


//===================================================================================
//=========================== Combustion Noise ======================================
//===================================================================================


void LinearFlowNoiseInt::CalcElemVec4CombustionTij(const Matrix<Double>& ptCoord,
                                                   const Matrix<Double>& NodalVel,
                                                   const Vector<Double>& NodalRho,
                                                   Vector<Double>& Result,
                                                   const Elem* elem) {
    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
    // T_ij = rho v_i v_j
    // Source term =  integral(grad[Sf].div[T])

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetNumRows();

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

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
    Vector<Double> helpVec;
    helpVec.Resize(dim);

    
    // Loop over all integration points 
    for (Integer actInt=1; actInt<=l; actInt++)
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

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times  vec
    // vec = partial / partial t ( (1-rho) * velocity )
    // Source term =  integral(grad[Sf].vec)

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetNumRows();

    Matrix<Double> xiDx;      
    Vector<Double>  Sf;
    Vector<Double> VecAtIP;

    Vector<double> intWeights = ptelem->GetIntWeights();

    Double jacDet;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
    Vector<Double> helpVec;
    helpVec.Resize(dim);

    
    // Loop over all integration points 
    for (Integer actInt=1; actInt<=l; actInt++)
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

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;

    Vector<Double> partResult;
    partResult.Resize(n);
    
    // Loop over all integration points 
    for (Integer actInt=1; actInt<=l; actInt++)
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

    // This functions computes the element RHS vector by integrating 
    // the gradient of the shape function times the divergence of the tensor T.
    // T_ij = rho v_i v_j
    // Source term =  integral(grad[Sf].div[T])

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetNumRows();

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

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;
    
    Vector<Double> helpVec;
    helpVec.Resize(dim);

    Vector<Double> nVec(dim);
    ComputeNormalVec( ptCoord, nVec);
    
    // Loop over all integration points 
    for (Integer actInt=1; actInt<=l; actInt++)
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

    // This functions computes the element RHS vector by integrating 
    // the shape function times scalar product of the vector with he normal vector
    // vector = \partial / \partial t ( (1-rho) v )

    Integer l = ptelem->GetNumIntPoints(); 
    Integer n = ptelem->GetNumNodes();
    Integer dim = ptCoord.GetNumRows();

    Matrix<Double> xiDx;      
    Vector<Double>  Sf;
    Vector<Double> VelAtIP;
    Vector<double> intWeights = ptelem->GetIntWeights();

    Double jacDet;

    Result.Resize(n);
    for (Integer i=0; i<n; i++)
      Result[i]=0.0;
    
    Vector<Double> nVec(dim);
    ComputeNormalVec( ptCoord, nVec);
    
    // Loop over all integration points 
    for (Integer actInt=1; actInt<=l; actInt++)
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

    // compute the normal vector
    UInt dim = nVec.GetSize();
    Vector<Double> vec1(dim);
    Vector<Double> vec2(dim);

    if ( dim == 3) {
      UInt idx=0;
      if ( ptCoord.GetNumCols() == 4) {
        idx = 3;
      }
      else if ( ptCoord.GetNumCols() == 3) {
        idx = 2;
      }
      else {
        EXCEPTION("Surface Element is no quadrilateral or triangle");
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

  nLinKuznetsovRHSInt::~nLinKuznetsovRHSInt() {}

  void nLinKuznetsovRHSInt::CalcElemVector( Vector<Double>& elemVec,
                                            EntityIterator& ent ) 
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt spaceDim = ptelem->GetDim();  

    elemVec.Resize(numFncs);
    elemVec.Init();

    solGradAtIp_.Resize(numFncs);
    solDeriv1GradAtIp_.Resize(numFncs);

    Double jacDet, totalfactor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = 0;
      ptelem->GetShFncAtIp(ShpFncAtIp_,actIntPt, ent.GetElem() );
      ptelem->GetGlobDerivShFncAtIp(xiDx_, actIntPt, ptCoord_, 
                                    jacDet, ent.GetElem() );
        
      if (isaxi_) {
        CoordAtIp_ = ptCoord_ * ShpFncAtIp_;
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
//       solGradAtIp_       = xiDxTransp * sol_;
//       solDeriv1GradAtIp_ = xiDxTransp * solderiv1_;


      //get 1st and 2nd derivative of solution at integration point
      solDeriv1AtIp_ = solderiv1_ * ShpFncAtIp_;
      solDeriv2AtIp_ = solderiv2_ * ShpFncAtIp_;
        
      totalfactor=0;
      const unsigned int xiDxcols(xiDx_.GetNumCols());
      for (UInt j=0; j<xiDxcols; j++)
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

  nLinWesterveltRHSInt::~nLinWesterveltRHSInt() {}

  void nLinWesterveltRHSInt::CalcElemVector( Vector<Double>& elemVec,
                                             EntityIterator& ent )
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
       
    elemVec.Resize(numFncs);
    elemVec.Init();

    Double jacDet, totalfactor;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent.GetElem() );
      ptelem->GetShFncAtIp(ShpFncAtIp_,actIntPt, ent.GetElem() );
        
      if (isaxi_) {
        CoordAtIp_ = ptCoord_ * ShpFncAtIp_;
        jacDet *= 2 * PI * CoordAtIp_[0];
      }
        
      //get solution and derivatives at integration point
      solAtIp_ = ShpFncAtIp_ * sol_;
      solDeriv1AtIp_ =  ShpFncAtIp_ * solderiv1_;
      solDeriv2AtIp_ =  ShpFncAtIp_ * solderiv2_;
        
      totalfactor = jacDet * 2.0 * factor_ * (solDeriv1AtIp_ * solDeriv1AtIp_ +
                                              solAtIp_ * solDeriv2AtIp_);
        
      for (UInt i=0; i< numFncs; i++)
        elemVec[i] += ShpFncAtIp_[i] * totalfactor;
    }
  }

  // =============================================================================
  // acoustic power source in heat conduction PDE
  // =============================================================================
  linAcouPowerSourceInt::linAcouPowerSourceInt( bool isaxi, 
                                                bool isharmonic,
                                                const std::string& readerId,
                                                const std::string& regionName,
                                                Double density ) :
    LinearForm(),
    readerId_(readerId),
    isharmonic_(isharmonic),
    density_(density),
    actStep_(0),
    lastStep_(0),
    sequenceStep_(1) // assume, that acoustic results do not come from multi sequence analysis
  {
    name_ = "linAcouPowerSourceInt";
    isaxi_ = isaxi;
    
    regionNames_.Push_back(regionName);

    if ( isharmonic_ )
    {
      // in harmonic case, results are from single frequency step
      UInt stepValue = 1;
      
      res_C_ = domain->GetResultHandler()->GetStoreSol<Complex>( readerId_,
          sequenceStep_, stepValue, ACOU_POTENTIAL, regionNames_ );
      resD1_C_ = domain->GetResultHandler()->GetStoreSol<Complex>( readerId_, 
          sequenceStep_, stepValue, ACOU_POTENTIAL_DERIV_1, regionNames_ );
      
    } else {
      // Specify that we want to use the current step number.
      mParser_->SetExpr( mHandle_, "step" );
      
      // assume, that first calculated time step is 1
      lastStep_ = 0;
    }
  }

  linAcouPowerSourceInt::~linAcouPowerSourceInt() {}
  
  void linAcouPowerSourceInt::CalcElemVectorDouble( Vector<Double>& elemPower,
                                                    EntityIterator& ent )
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const UInt spacedim = ptelem->GetDim();
    //const UInt nrNodes  = ptelem->GetNumNodes();

    Vector<Double> elemSol, elemSolD1, solGradAtIp, solD1GradAtIp, N2;
    Double solD1AtIp, N1;

    // get current time step of transient analysis
    actStep_ = (UInt) mParser_->Eval(mHandle_);
    
    if ( lastStep_ < actStep_ ) {
      // get result only once per time step
      res_R_ = domain->GetResultHandler()->GetStoreSol<Double>( readerId_,
           sequenceStep_, actStep_, ACOU_POTENTIAL, regionNames_ );
      resD1_R_ = domain->GetResultHandler()->GetStoreSol<Double>( readerId_, 
           sequenceStep_, actStep_, ACOU_POTENTIAL_DERIV_1, regionNames_ );
      
      lastStep_ = actStep_;
    }
    
    // retrieve acoustic solution, i.e. acouPotential and acouPotentialD1 of element
    res_R_->GetElemSolution(elemSol, ent);
    resD1_R_->GetElemSolution(elemSolD1, ent);
    
    // resize and initialize output vector
    elemPower.Resize(numFncs);
    elemPower.Init();
    
    Double jacDet, factor;
    solGradAtIp.Resize(spacedim);
    solD1GradAtIp.Resize(spacedim);
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  

      jacDet = 0;
      ptelem->GetShFncAtIp( ShpFncAtIp_, actIntPt ,  ent.GetElem() );
      ptelem->GetGlobDerivShFncAtIp( xiDx_ , actIntPt, ptCoord_, 
                                     jacDet, ent.GetElem() );
       
      if (isaxi_) {
        CoordAtIp_ = ptCoord_ * ShpFncAtIp_;
        jacDet *= 2 * PI * CoordAtIp_[0];
      }
      
      // gradient of acouPotential and acouPotentialD1 at integration point
      solGradAtIp.Init();
      solD1GradAtIp.Init();
      for ( UInt i=0; i<spacedim; i++ ) {
        for ( UInt j=0; j<numFncs; j++ ) {
          solGradAtIp[i] += xiDx_[j][i] * elemSol[j];
          solD1GradAtIp[i] += xiDx_[j][i] * elemSolD1[j];
        }
      }
      // acouPotentialD1 at integration point
      solD1AtIp = ShpFncAtIp_ * elemSolD1;
      
      // calculate helper variables N1 and N2
      N1 = 0;
      for (UInt k=0; k<spacedim; k++)
        N1 += solGradAtIp[k] * solD1GradAtIp[k];
      
      N2 = solGradAtIp * solD1AtIp;
      // std::cerr << "N1:\n" << N1 << "\nN2:\n" << N2 << "\n";
      
      factor = density_ * jacDet * intWeights[actIntPt-1];
      // source element vector
      for (UInt i=0; i< numFncs; i++) {
        
        elemPower[i] += factor * ShpFncAtIp_[i] * N1;
        
        for (UInt j=0; j<spacedim; j++)
          elemPower[i] -= factor * xiDx_[i][j] * N2[j];
      }
    }
    // std::cerr << "elemPower = " << elemPower << "\n\n";
  }

  
  void linAcouPowerSourceInt::CalcElemVectorComplex( Vector<Double>& elemPower,
                                                     EntityIterator& ent )
  {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();  
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const UInt spacedim = ptelem->GetDim();  
    // const UInt nrNodes  = ptelem->GetNumNodes();
    
    Vector<Complex> elemSol, elemSolD1, solGradAtIp, solD1GradAtIp, N2;
    Complex solD1AtIp, N1;
    
    // get element solution
    res_C_->GetElemSolution(elemSol, ent);
    resD1_C_->GetElemSolution(elemSolD1, ent);
    
    // ----------------------------------------------------------------------
//    std::cerr << "elemSol:\n";
//    for (UInt i=0; i<elemSol.GetSize();i++)
//      std::cerr << "  abs = " << std::abs(elemSol[i]) << "  angle = "
//                << std::atan2( elemSol[i].imag(), elemSol[i].real() )*180/PI << "\n";
//    std::cerr << "elemSolD1:\n";
//    for (UInt i=0; i<elemSolD1.GetSize();i++)
//      std::cerr << "  abs = " << std::abs(elemSolD1[i]) << "  angle = "
//                << std::atan2( elemSolD1[i].imag(), elemSolD1[i].real() )*180/PI << "\n";
    // ----------------------------------------------------------------------
    
    
    // resize and initialize output vector
    elemPower.Resize(numFncs);
    elemPower.Init();
    
    Double jacDet, factor;
    N2.Resize(spacedim);
    solGradAtIp.Resize(spacedim);
    solD1GradAtIp.Resize(spacedim);
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {

      jacDet = 0.0;
      ptelem->GetShFncAtIp( ShpFncAtIp_, actIntPt ,  ent.GetElem() );
      ptelem->GetGlobDerivShFncAtIp( xiDx_ , actIntPt, ptCoord_, 
                                     jacDet, ent.GetElem() );
      
      if (isaxi_) {
        CoordAtIp_ = ptCoord_ * ShpFncAtIp_;
        jacDet *= 2 * PI * CoordAtIp_[0];
      }
      
      // gradient of acouPotential and acouPotentialD1 at integration point
      solGradAtIp.Init();
      solD1GradAtIp.Init();
      for ( UInt i=0; i<spacedim; i++ ) {
        for ( UInt j=0; j<numFncs; j++ ) {
          solGradAtIp[i] += xiDx_[j][i] * elemSol[j];
          solD1GradAtIp[i] += xiDx_[j][i] * elemSolD1[j];
        }
      }
      // acouPotentialD1 at integration point
      solD1AtIp = ShpFncAtIp_ * elemSolD1;
      
      // calculate helper variables N1 and N2
      N1 = Complex(0.0, 0.0);
//      for (UInt k=0; k<spacedim; k++)
//        N1 += 0.5 * solD1GradAtIp[k] * std::conj( solGradAtIp[k] );
      // --------------------------------------------------------------------
//      std::cerr << "N1:\n" << "  abs = " << std::abs(N1) << "  angle = "
//                << std::atan2( N1.imag(), N1.real() )*180/PI << "\n";
      // --------------------------------------------------------------------

      for (UInt k=0; k<spacedim; k++)
        N2[k] = 0.5 * std::conj( solGradAtIp[k] ) * solD1AtIp;
      // --------------------------------------------------------------------
//      std::cerr << "N2:\n";
//      for (UInt i=0; i<N2.GetSize(); i++)
//        std::cerr << "  abs = " << std::abs(N2[i]) << "  angle = "
//                  << std::atan2( N2[i].imag(), N2[i].real() )*180/PI << "\n";
      // --------------------------------------------------------------------
        
      factor = density_ * jacDet * intWeights[actIntPt-1];
      for (UInt i=0; i< numFncs; i++) {
        
        elemPower[i] += factor * ShpFncAtIp_[i] * N1.real();
        
        for (UInt j=0; j<spacedim; j++)
          elemPower[i] -= factor * xiDx_[i][j] * N2[j].real();
      }
    }
    // std::cerr << "elemPower:\n" << elemPower << "\n\n";
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
   
    // const Vector<Double> & intWeights = ptelem->GetIntWeights(); // TODO: Unused variable intWeights 
    Vector<Double> ShpFncAtIp, CoordAtIP;
    Matrix<Double> xiDx;

    elemVec.Resize(numFncs);
    elemVec.Init(0);

    // Double factor; // TODO: Unused variable factor
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

      // factor = intWeights[actIntPt-1] * jacDet;
      for (UInt i=0; i<numFncs; i++) {
        for (UInt j=0; j<D_.GetSize(); j++) {
          elemVec[i] += xiDx[i][j] * D_[j];
        }
      }
    }

  }


  // =========================================================================
  // maxwell homogenization rhs volume integrator
  // =========================================================================


  VolChargeHomInt::VolChargeHomInt(BaseMaterial* matData, Global::ComplexPart matDataType, UInt numDof,
      const std::string& phase,
      bool isaxi) {

    matData_ = matData;
    name_ = "VolChargeHomInt";
    isaxi_ = isaxi;
    numDofs_ = numDof;
    phase_ = phase;
    locMatDataType_ = matDataType;


  }

  VolChargeHomInt::~VolChargeHomInt() {
    delete bilinearStiff_;

  }

  void VolChargeHomInt::SetVolChargeVector(StdVector<std::string> & volChargex,
      StdVector<std::string> & volChargey,
      StdVector<std::string> & volChargez,
      const CoordSystem * coordSys,
      bool isUnit, Double volume, const int dim ) {


    locChargex_ = volChargex;
    locChargey_ = volChargey;
    locChargez_ = volChargez;
    coordSys_ = coordSys;
    isUnitValue_ = isUnit;
    volume_ = volume;
    dim_ = dim;


    // flag for updatedLagrange formulation
    bool upLagrangeForm = true;

    //transform the type
    SubTensorType tensorType;

    if ( dim_ == 3 ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }

    bilinearStiff_  = new linGradBDBInt(matData_, ELEC_PERMITTIVITY, tensorType, upLagrangeForm);
    this->SetMatDataType(locMatDataType_);
    bilinearStiff_->SetMatDataType(locMatDataType_);

  }

  void VolChargeHomInt::CalcElemVector( Vector<double> & result,
      EntityIterator& ent )
  {
    Vector<Complex> help;
    VolChargeHomInt::CalcElemVector(help, ent);
    result = (help.GetPart(Global::REAL));
  }


  void VolChargeHomInt::CalcElemVector(Vector<Complex> &elemVec, EntityIterator &ent)
  {
    // get global coordinate system and math parser
    MathParser * parser = domain->GetMathParser();

    // Now evaluate each entry
    // Vector<T> locLoadVec( dim_ );
    Vector<Double> locLoadVec( dim_ );
    for ( UInt i = 0; i < locChargex_.GetSize(); i++ ) {
      parser->SetExpr( mHandle_, locChargex_[i] );
      locLoadVec[0] = parser->Eval( mHandle_ );
      parser->SetExpr( mHandle_, locChargey_[i] );
      locLoadVec[1] = parser->Eval( mHandle_ );
      if (dim_ == 3){
        parser->SetExpr( mHandle_, locChargez_[i] );
        locLoadVec[2] = parser->Eval( mHandle_ );
      }
    }

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo(ent);
    ptelem->SetAnsatzFct(ansatzFct1_);
    const Vector<double> &intWeights = ptelem->GetIntWeights();

    // extract pointer and get coords again for the integrator
    bilinearStiff_->SetAnsatzFct( ansatzFct1_ );
    bilinearStiff_->ExtractElemInfo(ent);

    // calc dMat and immediately calculate the element stress
    Matrix<double> dMat;
    bilinearStiff_->calcDMat(dMat, ent.GetElem());
    Vector<double> stress;
    stress = dMat * locLoadVec;

    elemVec.Resize(ptelem->GetNumNodes(), 0.0);
    Matrix<double> linBMat;
    Vector<double> partElemVec;
    Vector<Complex> tmp_vec(elemVec.GetSize());

    for(unsigned int aIP = 1, nr = ptelem->GetNumIntPoints(); aIP <= nr; ++aIP)
    {
      bilinearStiff_->calcBMatOnly(linBMat, aIP, ptelem, ptCoord_);
      linBMat.MultT(stress, partElemVec);

      const double jacDet(ptelem->CalcJacobianDetAtIp(aIP, ptCoord_, ent.GetElem()));

      partElemVec *= jacDet * intWeights[aIP-1];
      tmp_vec.SetPart(matDataType_, partElemVec);
      elemVec += tmp_vec;
    }

    LOG_DBG2(linForm) << "VolChargeHomInt:CEV el=" << ent.GetElem()->elemNum << " -> " << elemVec.ToString();

  }



  // =========================================================================
  // mechanic rhs volume integrator
  // =========================================================================
  
  
  VolForceInt::VolForceInt(UInt numDof,
                           const std::string& phase,
                           bool isaxi) {


    name_ = "VolForceInt";
    isaxi_ = isaxi;
    numDofs_ = numDof;
    phase_ = phase;
  }
    
  VolForceInt::~VolForceInt() {


  }
    
  void VolForceInt::SetVolForceVector(StdVector<std::string> & volForce,
                                      const CoordSystem * coordSys,
                                      bool isUnit, Double volume ) {


    locForce_ = volForce;
    coordSys_ = coordSys;
    isUnitValue_ = isUnit;
    volume_ = volume;

  }

  void VolForceInt::CalcElemVector( Vector<Double> & elemVec,
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

    //    std::cout << "ElemeVec: " << elemVec << std::endl;
  }

  void VolForceInt::CalcElemVector( Vector<Complex> & elemVec,
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
  void VolForceInt::CalcPartVector( Vector<TYPE>& elemVec, 
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
    elemVec.Init();
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
  : LinearForm(), matData_(matData), subTensorType_(type)
  {
    name_ = "AddStressRHSInt";

    switch(subTensorType_)
    {
    case AXI:
      isaxi_ = true;
      addStress_.Resize(4, 0.0);
      addStress_[0] = stressVec[0];
      addStress_[1] = stressVec[1];
      addStress_[2] = stressVec[5];
      addStress_[3] = stressVec[2];
      break;
    case PLANE_STRAIN:
      addStress_.Resize(3, 0.0);
      addStress_[0] = stressVec[0];
      addStress_[1] = stressVec[1];
      addStress_[2] = stressVec[5];
      break;
    case FULL:
      addStress_.Resize(6, 0.0);
      addStress_[0] = stressVec[0];
      addStress_[1] = stressVec[1];
      addStress_[2] = stressVec[2];
      addStress_[3] = stressVec[3];      
      addStress_[4] = stressVec[4];
      addStress_[5] = stressVec[5];
      break;
    default:
      EXCEPTION("Subtensortype not implemented");
    }
    
    bilinearStiff_ = new linElastInt(matData_, subTensorType_);
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

    // probably need to extract element info for bilinearStiff_ too! 
    // see AddStrainRHSInt::CalcElemVector

    partElemVec.Resize(nrNodes * nrDofs);
    partElemVec.Init();

    elemVec.Resize(nrNodes*nrDofs);
    elemVec.Init();


    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
    {    
      bilinearStiff_->CalcBMatOnly(linBMat, actIntPt, ptelem, ptCoord_);

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

  }

  
  
  // ==================================================================
  // linMech: add mechanical strain as force to RHS
  // ==================================================================

  AddStrainRHSInt::AddStrainRHSInt(BaseMaterial* matData,
             const Vector<double>& strainVec, SubTensorType type)
    : LinearForm(), matData_(matData), subTensorType_(type)
  {
    name_ = "AddStrainRHSInt";
    
    switch(subTensorType_)
    {
    case PLANE_STRAIN:
    case PLANE_STRESS: // TODO Verify!
      addStrain_.Resize(3, 0.0);
      addStrain_[0] = strainVec[0];
      addStrain_[1] = strainVec[1];
      addStrain_[2] = strainVec[5]; // voigt notation
      break;
    case FULL:
      addStrain_.Resize(6, 0.0);
      addStrain_[0] = strainVec[0];
      addStrain_[1] = strainVec[1];
      addStrain_[2] = strainVec[2];
      addStrain_[3] = strainVec[3];      
      addStrain_[4] = strainVec[4];
      addStrain_[5] = strainVec[5];
      break;
    default:
      EXCEPTION("Subtensortype not implemented");
    }
    
    bilinearStiff_  = new linElastInt(matData_, subTensorType_);
  }

  AddStrainRHSInt::~AddStrainRHSInt()
  {
    delete bilinearStiff_;
  }

  void AddStrainRHSInt::CalcElemVector(Vector<double> &elemVec, EntityIterator &ent, const Vector<Double>* test_strain)
  {
    assert(test_strain == NULL || ((domain->GetGrid()->GetDim() == 2 && test_strain->GetSize() == 3) || test_strain->GetSize() == 6));

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo(ent);
    ptelem->SetAnsatzFct(ansatzFct1_);
    const Vector<double> &intWeights = ptelem->GetIntWeights();
  
    // extract pointer and get coords again for the integrator
    bilinearStiff_->SetAnsatzFct( ansatzFct1_ );
    bilinearStiff_->ExtractElemInfo(ent);
    
    // calc dMat and immediately calculate the element stress
    Matrix<double> dMat;
    bilinearStiff_->calcDMat(dMat, ent.GetElem());
    Vector<double> stress;
    if(test_strain != NULL)
      stress = dMat * (*test_strain);
    else
      stress = dMat * addStrain_;

    elemVec.Resize(ptelem->GetNumNodes() * ptelem->GetDim(), 0.0);
    Matrix<double> linBMat;
    Vector<double> partElemVec;
        
    for(unsigned int aIP = 1, nr = ptelem->GetNumIntPoints(); aIP <= nr; ++aIP)
    {    
      bilinearStiff_->CalcBMatOnly(linBMat, aIP, ptelem, ptCoord_);
      linBMat.MultT(stress, partElemVec);

      const double jacDet(ptelem->CalcJacobianDetAtIp(aIP, ptCoord_, ent.GetElem()));

      partElemVec *= jacDet * intWeights[aIP-1];
      elemVec += partElemVec;
    }

    LOG_DBG2(linForm) << "AddStrainRHSInt:CEV el=" << ent.GetElem()->elemNum << " -> " << elemVec.ToString();

  }
  // =============================================================================
  // MagSourceInt (RHS integrator for Biot-Savart coil of magnetic Scalar PDE)
  // =============================================================================

  //! Calculates the RHS for Biot-Savart coils for scalar potential formulation

  //! This class implements the RHS excitation for a magnetic-scalar formulation
  //! of the magnetic field. The integrator calcualates
  //! \f [ \int \nabla N_a (\mu - \mu_0) \mathbf{H}_s d\Omega \f ]
  BiotSavartSourceInt::
  BiotSavartSourceInt( BaseMaterial* matData,
                       SubTensorType tensorType,
                       shared_ptr<BiotSavart>& bs,
                       bool isaxi, bool coordUpdate ) 
  : LinearForm( matData ) {
    name_ = "BiotSavartSourceInt";
    isaxi_ = isaxi;
    coordUpdate_ = coordUpdate;
    biotSavart_ = bs;
    mParser_->SetExpr(mHandle_, "1.0");
    subTensorType_ = tensorType; 
  }

  BiotSavartSourceInt::~BiotSavartSourceInt() {

  }
  
  void BiotSavartSourceInt::SetFactor( const std::string& factor ) {
    mParser_->SetExpr(mHandle_, factor);
  }

  void BiotSavartSourceInt::CalcElemVector( Vector<Double> & elemVec,
                                            EntityIterator& ent ) {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    Matrix<Double> permeability;
    const Double mu0=1.25664e-06;

    ptMaterial->GetTensor(permeability,MAG_PERMEABILITY,
                          Global::REAL, subTensorType_ );
    
    

    // as we add the contribution of the Biot-Savart field afterwards in
    // the magnetic scalar potential PDE, we just calculate for the "difference",
    // i.e. we have to subtract mu_0 from the permeability here.
    UInt muSize = permeability.GetNumRows();
    for( UInt i = 0; i < muSize; i++ ) {
      permeability[i][i] -= mu0;
    }
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    const Vector<Double> * intPoints = ptelem->GetIntPoints();

    Vector<Double> ShpFncAtIp, helpVec, CoordAtIP, hField, bField;
    Matrix<Double> xiDx;

    elemVec.Resize(numFncs);
    elemVec.Init(0);  
    helpVec.Resize(numFncs);
    helpVec.Init(0);

    Double jacDet = 0.0;
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      // get the derivative of the global shape function 
      ptelem->GetGlobDerivShFncAtIp( xiDx, actIntPt, ptCoord_, 
                                     jacDet, ent.GetElem() );
      ptelem->Local2GlobalCoord( CoordAtIP, intPoints[actIntPt-1],
                                 ptCoord_, ent.GetElem() );
      
      // Calculate source field H_s using Biot-Savart's law
      biotSavart_->CalcFieldAtPoint(hField, CoordAtIP );
      bField = (permeability * hField);
      helpVec = xiDx * bField;
      helpVec *=  intWeights[actIntPt-1] * jacDet;
      elemVec += helpVec   ;           
    }
    // multiply element vector by factor (e.g. -1 due to magnetostrictive 
    // coupling)
    elemVec *= mParser_->Eval( mHandle_ );//for the magnetostriction coupling case
  }

  // ============================================================================
  // BiotSavarMechCouplingRHSInt
  // =============================================================================
  BiotSavartMechCouplingInt::
  BiotSavartMechCouplingInt( BaseMaterial* matData,
                             SubTensorType type,
                             shared_ptr<BiotSavart> & bs,
                             bool isaxi,
                             bool coordUpdate )
  : LinearForm( matData ) {
    
    name_ = "BiotSavartMechCouplingInt";
    subTensorType_ = type;
    isaxi_ = false;
    biotSavart_ = bs;
    
    // create new magnetostriction bilinear form
    magStrictForm_ = new LinMagStrictInt( ptMaterial, subTensorType_ );
  }

  BiotSavartMechCouplingInt::~BiotSavartMechCouplingInt() {
    delete magStrictForm_;
  }

  void BiotSavartMechCouplingInt::CalcElemVector( Vector<Double> & elemVec,
                                                  EntityIterator& ent ) {
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent );

    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrIntPts = ptelem->GetNumIntPoints();
    const UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt dim = magStrictForm_->getNumDofsA();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    const Vector<Double> * intPoints = ptelem->GetIntPoints();
    Double jacDet;
    Matrix<Double> aMat, dMat, bMatTrans,aMatTrans;
    Matrix<Double> adMat;

    elemVec.Resize( numFncs * dim);
    elemVec.Init(0);

    // compute the correct material tensor
    Vector<Double> lCoord, hField;
    Vector<Double>  ShpFncAtIp;

    magStrictForm_->ExtractElemInfo( ent );
    for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {     

      // Setup the A matrix for current integration point
      magStrictForm_->calcAMat( aMat, actIntPt, ptCoord_ );
      magStrictForm_->calcDMat( dMat,ent.GetElem());
      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_,ent.GetElem() );

      ptelem->Local2GlobalCoord( lCoord, intPoints[actIntPt-1],ptCoord_, ent.GetElem() );
      
      // Calculate Biot-Savart field H_s at integration point
      biotSavart_->CalcFieldAtPoint( hField, lCoord );
      ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent.GetElem() );

      if( isaxi_ ) {
        EXCEPTION("axi-symmetric implementation missing");
      }
      jacDet *= intWeights[actIntPt-1];
      // compute intermediate matrix (B^T * D)
      adMat = aMat * dMat;

      // We now compute  (B^T * D) * H_{s} and scale it by the determinant
      // of the Jacobian 
      for ( UInt i = 0; i < adMat.GetNumRows(); i++ ) {
        for ( UInt j = 0; j < adMat.GetNumCols(); j++ ) {
          elemVec[i] +=  adMat[i][j]*hField[j]*jacDet;
        }
      }
    }
  }
} // end of namespace
