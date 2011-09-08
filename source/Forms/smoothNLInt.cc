#include <iostream>
#include <fstream>

#include "Utils/nodestoresol.hh"

#include "smoothNLInt.hh"
#include "mechStressStrain.hh"


namespace CoupledField
{

  // ===================================================================================
  // =================== standard con- and destructors (just for tracing) ==============
  // ===================================================================================

  SmoothNLInt::SmoothNLInt(BaseMaterial* matData, 
                           Matrix<Double> & elastWeights,
                           Matrix<Double> & couplingNodes,
                           Double AMax,
                           Double AMin,
                           Double characteristicLength,
                           Double expo,
                           SubTensorType type,
                           bool coordUpdate ) 
    : BDBInt(matData, type), elastWeights_(elastWeights), couplingNodes_(couplingNodes), AMax_(AMax), AMin_(AMin),characteristicLength_(characteristicLength),exponent_(expo), elastFactorMax_(0.0), elastFactorMin_(1.0e99), eMean_(1.0)
  {

    name_ = "SmoothNLInt";
    updateDMatInEveryIP_ = false;
    coordUpdate_ = coordUpdate;
    isSolDependent_ = true;

    //file.open("ElastWeights.dat", std::ios::out);

    if ( type == FULL ) {
      dimD_   = 6;
      nrDofs_ = 3;
    }
    else if ( type == PLANE_STRAIN ) {
      dimD_   = 3;
      nrDofs_ = 2;
    }
    else if ( type == AXI ) {
      dimD_   = 4;
      nrDofs_ = 2;
      isaxi_  = true;
    }
  }
 

  SmoothNLInt::~SmoothNLInt()
  {
  }


  void SmoothNLInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                       EntityIterator& ent1, 
                                       EntityIterator& ent2 ) {
    

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    UInt elemNumber = ent1.GetElem()->elemNum;
    Double elastFactor;
    // TODO: Check if this is still needed
    // Double elastFactor2, elastFactor3;
    //Double LMaxMax;
    //Double LMinMax;
    bool firstTime, firstIter;
    Double deltaX, deltaY;
    Double dist, distNew;
    Double maxDist, exponent;
    //ptelem->GetMaxMinEdgeLength(ptCoord_,maxEdgeLength_,minEdgeLength_);

    firstTime=false;
    firstIter=false;

    if(elastWeights_[elemNumber-1][4] < -0.001 && elastWeights_[elemNumber-1][3] < -0.001){
      firstTime=true;
      firstIter=true;
    }
    else if(elastWeights_[elemNumber-1][4] > 0.001 && elastWeights_[elemNumber-1][3] < -0.001){
      firstIter=true;
    }

//     file << elastWeights_[elemNumber-1][4]
//          << " " << elastWeights_[elemNumber-1][3];
// TODO: szoerner: check the smoothing, maybe better if no distance. Volume
// calculation also needs a check or rather should be removed
      
    if(firstTime && firstIter){
      Vector<Double> intPoint, globIntPoint;
    
      ptelem->GetCoordMidPoint(intPoint);
      ptelem->Local2GlobalCoord(globIntPoint, intPoint,ptCoord_,ent1.GetElem() );

      deltaX=couplingNodes_[0][0]-globIntPoint[0];
      deltaY=couplingNodes_[0][1]-globIntPoint[1];

      dist=std::sqrt(deltaX*deltaX + deltaY*deltaY);
        
      for(UInt i=1; i<couplingNodes_.GetNumRows(); i++){
        deltaX=couplingNodes_[i][0]-globIntPoint[0];
        deltaY=couplingNodes_[i][1]-globIntPoint[1];
        //deltaZ=couplingNodes_[i][2]-zMid;
        distNew=std::sqrt(deltaX*deltaX + deltaY*deltaY);
        if(distNew<dist)
          dist=distNew;
      }

      Double auxFac=1.0;
      if(globIntPoint[0]<0.061 && abs(globIntPoint[1])<0.002)
        auxFac=1.0;
      maxDist=characteristicLength_;
      exponent=exponent_;
      elastFactor=std::pow((dist/maxDist),exponent);

      if(elastFactor<1e-10)
        elastFactor=1e-10;

      elastWeights_[elemNumber-1][1]=(elastFactor);
      elastWeights_[elemNumber-1][2]=(elastFactor);
      elastWeights_[elemNumber-1][3]=1.0;
      elastWeights_[elemNumber-1][4]=1.0;
    } else if(!firstTime && firstIter){
      elastFactor=elastWeights_[elemNumber-1][2];
      elastWeights_[elemNumber-1][1]=elastFactor;
      elastWeights_[elemNumber-1][3]=1.0;
    } else if(!firstTime && !firstIter) {

      Matrix<Double> elemResult_;
    
      sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );
      
      std::auto_ptr<MechStressStrain<Double> > strain =
        std::auto_ptr<MechStressStrain<Double> >(new MechStressStrain<Double>(ptMaterial,subTensorType_));
      
      strain->SetActElemSol(elemResult_);
      
      Vector<Double> actStrain;

      Vector<Double> intPoint;

      ptelem->GetCoordMidPoint(intPoint);
      
      strain->SetIntPoint(intPoint);
      
      strain->CalcStrainVec(actStrain,1,ent1);

      Matrix<Double> newCoord;
      newCoord=ptCoord_+elemResult_;

      elastFactor=ComputeElastFactor(actStrain[0],actStrain[1],actStrain[2] );

      //3000 for Voice Model 
      elastFactor*=3000.0;

      elastFactor*=elastWeights_[elemNumber-1][2];

      elastFactor+=(elastWeights_[elemNumber-1][1]);
      //if( Ai/A0 < 0.05 ){
        //elastFactor=(1.0*(1.0-elastFactor))+1.0;
        //elastFactor=A0/Ai;
        //if (elastFactor>1.1)
          //elastFactor=1.1;
        //elastFactor*=1.1;
        //elastFactor3=elastFactor;
        //elastFactor*=elastWeights_[elemNumber-1][1];
        //elastWeights_[elemNumber-1][1]=elastFactor;
      //}
      //else{
        //elastFactor=elastWeights_[elemNumber-1][1];
      //}
    }

    // First of all, set ansatz function to element
    ptelem->SetAnsatzFct( ansatzFct1_ );
    
    UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
    const UInt nrDofs   = getNrDofs();  
    
    double jacDet;
    
    Matrix<Double> bMat; 
    Matrix<Double> dMat; 
    Matrix<Double> dbMat; 
    Double aux, fac, *ptr1, *ptr2;

    elemMat.Resize( nrFncs * nrDofs);
    elemMat.Init();

    dbMat.Resize( getDimD(), nrFncs * nrDofs);

    const UInt nrIntPts = ptelem->GetNumIntPoints(); 
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  

    // // Check if material has to be rotated
    calcDMat( dMat );

    // Loop over all integration points
    for ( UInt actIntPt = 1; actIntPt <= nrIntPts; actIntPt++ ) {

      // Setup the B matrix for current integration point
      CalcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        EXCEPTION("SmoothNLInt::CalcElementMatrix: Encountered "
                   << "negative Jacobian determinant!");
      }

      // Special things must be done in the axi-symmetric case
      if ( isaxi_ ) {
        Vector<Double> ShpFncAtIp;
        Vector<Double> CoordAtIP;
        ptelem->GetShFncAtIp( ShpFncAtIp, actIntPt, ent1.GetElem() );

        CoordAtIP = ptCoord_ * ShpFncAtIp;
        jacDet *= 2 * PI * CoordAtIP[0];
      }

      // Compute the matrix product D * B and store as intermediate matrix
      dMat.Mult( bMat, dbMat );

      // We now compute B^T * D * B and scale it by the determinant
      // of the Jacobian and the weight of the current integration
      // point. The result is added to the element matrix.
      fac = jacDet * intWeights[actIntPt-1] * elastFactor;
      for ( UInt k = 0; k < bMat.GetNumRows(); k++ ) {
        ptr1 =  bMat[k];
        ptr2 = dbMat[k];
        for ( UInt i = 0; i < bMat.GetNumCols(); i++ ) {
          aux = fac * ptr1[i];
          for ( UInt j = 0; j < dbMat.GetNumCols(); j++ ) {
            elemMat[i][j] += aux * ptr2[j];
          }
        }
      }
    }
  }

  // returns B - matrix for BDB
  void SmoothNLInt::CalcBMat(Matrix<Double> & bMat, UInt ip, const Matrix<Double> & ptCoord)
  {
    
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  

    UInt actDim, actNode, j, k;
    
    
    bMat.Resize(getDimD(), numFncs * nrDofs);
    bMat.Init();
    
    // local shape functions derived after global coords (format: numFncs x spaceDim)
    Matrix<Double> xiDx;

    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem());


    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < numFncs; actNode++)
        bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
        j = 1;
        k = 0;
        
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
            bMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
          }
        break;


      case 3:
        UInt actDim=spaceDim;
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
          }

        actDim++;
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
          }

        actDim++;
        for (actNode = 0; actNode < numFncs; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
          }
        break;
      }
  }


  void SmoothNLInt::calcDMat(Matrix<Double> & dMat)
  {

    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,matDataType_,subTensorType_);

  }  

  // calculates the D-matrix 
  void SmoothNLInt::calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord)
  {
    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,matDataType_,subTensorType_);
  }

  void SmoothNLInt::CalcMinMaxStrain( EntityIterator& ent1, 
                                      EntityIterator& ent2 ){
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );
    Double elastFactor;

    Matrix<Double> elemResult_;
    
    sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

    std::auto_ptr<MechStressStrain<Double> > strain =
      std::auto_ptr<MechStressStrain<Double> >(new MechStressStrain<Double>(ptMaterial,subTensorType_));
    
    strain->SetActElemSol(elemResult_);
    
    Vector<Double> actStrain;
    
    Vector<Double> intPoint;
    
    ptelem->GetCoordMidPoint(intPoint);
    
    strain->SetIntPoint(intPoint);
    
    strain->CalcStrainVec(actStrain,1,ent1);

    elastFactor = ComputeElastFactor(actStrain[0], actStrain[1], actStrain[2]);

    if (elastFactor>elastFactorMax_)
      elastFactorMax_=elastFactor;
    if (elastFactor<elastFactorMin_)
      elastFactorMin_=elastFactor;

  }

  Double SmoothNLInt::ComputeElastFactor( const Double& e11, const Double& e22, const Double& e12 )
  {
    Double fac;

    //Strategy 5
    //fac=std::sqrt( ((e11*e11)+(e22*e22))/(2.0) )/eMean_;
    //Strategy 5 mod
    fac=std::sqrt( ((e11*e11)+(e22*e22)+(e12*e12))/(3.0) )/eMean_;

    return fac;
  }

} // end namespace CoupledField
