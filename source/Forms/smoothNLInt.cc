#include <iostream>
#include <fstream>

#include "smoothNLInt.hh"
#include "mechStressStrain.hh"

#include "Utils/nodestoresol.hh"

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
    Double elastFactor, elastFactor2, elastFactor3, e11=0.0, e22=0.0, e12=0.0, A_elem;
    //Double LMaxMax;
    //Double LMinMax;
    bool firstTime, firstIter;
    Double deltaX, deltaY;
    Double dist, distNew;
    Double maxDist, exponent;
    A_elem = ptelem->CalcVolume(ptCoord_, isaxi_);
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
      
    if(firstTime && firstIter){
      Vector<Double> intPoint, globIntPoint;
    
      ptelem->GetCoordMidPoint(intPoint);
      ptelem->Local2GlobalCoord(globIntPoint, intPoint,ptCoord_,ent1.GetElem() );

      deltaX=couplingNodes_[0][0]-globIntPoint[0];
      deltaY=couplingNodes_[0][1]-globIntPoint[1];

      dist=std::sqrt(deltaX*deltaX + deltaY*deltaY);
        
      for(UInt i=1; i<couplingNodes_.GetSizeRow(); i++){
        deltaX=couplingNodes_[i][0]-globIntPoint[0];
        deltaY=couplingNodes_[i][1]-globIntPoint[1];
        //deltaZ=couplingNodes_[i][2]-zMid;
        distNew=std::sqrt(deltaX*deltaX + deltaY*deltaY);
        if(distNew<dist)
          dist=distNew;
      }
//         dist*
//     2_____________3
//     |             |
//     |_____________|
//     5             4
      //dist2
      deltaX=0.061-globIntPoint[0];
      deltaY=0.002-globIntPoint[1];
      Double dist2=std::sqrt(deltaX*deltaX + deltaY*deltaY);
      //dist3
      deltaX=0.071-globIntPoint[0];
      deltaY=0.002-globIntPoint[1];
      Double dist3=std::sqrt(deltaX*deltaX + deltaY*deltaY);
      //dist4
      deltaX=0.071-globIntPoint[0];
      deltaY=-0.002-globIntPoint[1];
      Double dist4=std::sqrt(deltaX*deltaX + deltaY*deltaY);
      //dist5
      deltaX=0.061-globIntPoint[0];
      deltaY=-0.002-globIntPoint[1];
      Double dist5=std::sqrt(deltaX*deltaX + deltaY*deltaY);

      deltaX=0.066-globIntPoint[0];
      deltaY=globIntPoint[1];
      Double distM=std::sqrt(deltaX*deltaX + deltaY*deltaY);


      Double auxFac=1.0;
      if(globIntPoint[0]<0.061 && abs(globIntPoint[1])<0.002)
        auxFac=1.0;
      maxDist=characteristicLength_;
      exponent=exponent_;
      //elastFactor=1f.0/(dist/maxDist);
      //if(globIntPoint[0]>0.071)
      //  exponent=2*exponent_;
        //maxDist=0.213;
      //else
      elastFactor=std::pow((dist/maxDist),exponent);
      //elastFactor=auxFac*(std::pow((dist/maxDist),exponent)) 
      //           +3.0*(std::pow((dist2/(0.1*maxDist)),2.0*exponent))
      //           +30.0*(std::pow((dist3/(0.1*maxDist)),2.0*exponent))
      //           +30.0*(std::pow((dist4/(0.1*maxDist)),2.0*exponent))
      //           +3.0*(std::pow((dist5/(0.1*maxDist)),2.0*exponent))
      //           +(std::pow((distM/(0.6*maxDist)),2.0*exponent));
      //if(dist<0.011){
        //elastFactor=-1.258e5*dist+1845.2;
        // potenz=(-2):  elastFactor=-1.0856e4*dist+179.13;
      //}
      //else{
        //elastFactor=std::pow((dist/maxDist),-3); 
        // potenz=(-2):  elastFactor=std::pow((dist/maxDist),-2);
      //}
      //if(globIntPoint[0]>0.011 && globIntPoint[0]<0.061)
      //Double elastFactor2 = globIntPoint[0];

      if(elastFactor<1e-10)
        elastFactor=1e-10;

      elastWeights_[elemNumber-1][1]=(elastFactor);
      elastWeights_[elemNumber-1][2]=(elastFactor);
      elastWeights_[elemNumber-1][3]=1.0;
      elastWeights_[elemNumber-1][4]=1.0;
    }
    else if(!firstTime && firstIter){
      elastFactor=elastWeights_[elemNumber-1][2];
      elastWeights_[elemNumber-1][1]=elastFactor;
      elastWeights_[elemNumber-1][3]=1.0;
    }
    
    else if(!firstTime && !firstIter) {

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
      Double A0, Ai;
      A0=ptelem->CalcVolume(ptCoord_, isaxi_);
      Ai=ptelem->CalcVolume(newCoord, isaxi_);

      //std::cout << "A rel = " << Ai/A0 << std::endl;

//       e11=abs(actStrain[0]);
//       e22=abs(actStrain[1]);
//       e12=abs(actStrain[2]);
//       e11=actStrain[0];
//       e22=actStrain[1];
//       e12=actStrain[2];

      e11=actStrain[0];
      e22=actStrain[1];
      e12=actStrain[2];


      elastFactor=ComputeElastFactor(e11,e22,e12 );
      //elastFactor2=elastFactor;
      //elastFactor/=elastFactorMax_;

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

    else{
      (*error) << " Wrong Time iteration flags of smooother!\n"
               << "firstTime=" 
               << firstTime
               << "; firstIter" 
               << firstIter;
      Error( __FILE__, __LINE__ );
    }

//     file << " " << elastFactor
//          << " " << elastFactorMax_
//          << " " << elastFactorMin_
//          << " " << e11
//          << " " << e22
//          << " " << e12
//          << std::endl;

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
      calcBMat( bMat, actIntPt, ptCoord_ );

      // Compute Jacobian for integration point
      jacDet = ptelem->CalcJacobianDetAtIp( actIntPt, ptCoord_, ent1.GetElem() );

      // Perform a safety check
      if ( jacDet < 0.0 ) {
        (*error) << "SmoothNLInt::CalcElementMatrix: Encountered "
                 << "negative Jacobian determinant!";
        Error( __FILE__, __LINE__ );
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
      for ( UInt k = 0; k < bMat.GetSizeRow(); k++ ) {
        ptr1 =  bMat[k];
        ptr2 = dbMat[k];
        for ( UInt i = 0; i < bMat.GetSizeCol(); i++ ) {
          aux = fac * ptr1[i];
          for ( UInt j = 0; j < dbMat.GetSizeCol(); j++ ) {
            elemMat[i][j] += aux * ptr2[j];
          }
        }
      }
    }
  }

  // returns B - matrix for BDB
  void SmoothNLInt::calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
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
    Double elastFactor, e11=0.0, e22=0.0, e12=0.0;

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

    e11=abs(actStrain[0]);
    e22=abs(actStrain[1]);
    e12=abs(actStrain[2]);

    elastFactor=ComputeElastFactor(e11,e22,e12 );

    if (elastFactor>elastFactorMax_)
      elastFactorMax_=elastFactor;
    if (elastFactor<elastFactorMin_)
      elastFactorMin_=elastFactor;

  }

  Double SmoothNLInt::ComputeElastFactor( Double & e11, Double & e22, Double & e12 ){
    Double fac;
    Double e11Sqrd, e22Sqrd, e12Sqrd;
    //Strategy 0
    //if(Sx<=5.5e-3)
    //  elastFactor = ((64.516*Sx) + 0.6452);
    //else
    //  elastFactor = ((-40.816*Sx)+1.2245);
    //elastFactor*=2.0;
    //elastFactor+=1.0;

    

    //Strategy 5
    //fac=std::sqrt( ((e11*e11)+(e22*e22))/(2.0) )/eMean_;
    //Strategy 5 mod
    fac=std::sqrt( ((e11*e11)+(e22*e22)+(e12*e12))/(3.0) )/eMean_;

    //Strategy 6
    //elastFactor=(1.0-nu)/(2.0*(1.0-2.0*nu)*(1.0+nu)) * ((e11*e11)+(e22*e22)+(2*nu/(1.0-nu)*e11*e22 ));

    //Strategy 7
    //elastFactor=(e11-e22)*(e11-e22)/(4.0*(1.0+nu));

    //Strategy 8
    //elastFactor=((e11*e11)+(e22*e22))/(2.0*eMean_*eMean_);


    //if(elastFactor<1.0e-5)
    //  elastFactor=1.0e-5;
    //else if(elastFactor>1.0)
    //  elastFactor=1.0;
    
    //if(e11>0.8)
    //  elastFactor*=(e11+0.25);
    //else if(e22>0.8)
    //  elastFactor*=(e22*0.25);
    //if(e12>0.5)
    //  elastFactor*=(e12+0.7);

    return fac;
  }

} // end namespace CoupledField
