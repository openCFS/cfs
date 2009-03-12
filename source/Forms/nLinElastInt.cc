// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "Forms/nLinElastInt.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField
{



  // =============================================================================
  // base class for nonlinear mechanics  
  // =============================================================================




  // base class for nonlinear calculation of elasticity
  nLinElastInt::nLinElastInt(BaseMaterial*  matData) 
    : linElastInt(matData)
  {
    name_ = "nLinElastInt";
    isSolDependent_ = true;
  }
 

  nLinElastInt::~nLinElastInt()
  {
  }


  void nLinElastInt::CalcElementMatrix( Matrix<Double>& elemMat,
                                        EntityIterator& ent1, 
                                        EntityIterator& ent2 ) {

    // get displacements of element

    sol_->GetElemSolutionAsMatrix( elemDisp_, ent1  );
    
    // call method of base class
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
  }


  void nLinElastInt::SetActEntities( EntityIterator ent1,
                                     EntityIterator ent2 ) {
    it1_ = ent1;
    it2_ = ent2;
  }

  // returns nonlinear B - matrix (first part) for BDB
  void nLinElastInt::calcBMat(Matrix<Double> & bMat, UInt ip, 
                              Matrix<Double> & ptCoord )
  {

    
    if (!elemDisp_.GetNumRows() || !elemDisp_.GetNumCols()) 
      EXCEPTION("Undefined displacements!");

    ptelem->SetAnsatzFct( ansatzFct1_ );
    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
   
    const UInt nrDofs  = getNrDofs();    
    const UInt dimD    = getDimD();

    
    
    Matrix<Double> linBMat;    
    linElastInt::calcBMat( linBMat, ip, ptCoord);


    bMat.Resize( dimD, numFncs * nrDofs );
    bMat.Init();

    
    // local shape functions derived after global coords 
    // (format: numFncs x nrDofs)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem() );


    // calculate derivates of global displacements at element nodes
    // elemDisp_ holds displacement of the actual element (dimension: nrDofs x numFncs)
    // xiDx holds derivatives of shape functions after global coords in one IP (dimension: numFncs x nrDofs)
    // displDeriv dimension: nrDofs x nrDofs
    Matrix<Double>  displDeriv;  
    displDeriv = elemDisp_ * xiDx;

    // we need transposed of displacement derivatives
    Matrix<Double> displDerivTransp;
    displDeriv.Transpose(displDerivTransp);  


    Matrix<Double> bMatOneNode;
    // size_row() = nr of rows!!
    bMatOneNode.Resize(linBMat.GetNumRows(), nrDofs );
    bMatOneNode.Init();

    
    Matrix<Double> helpMat;
    for( UInt actNode=0; actNode < numFncs; actNode++)
      {
        linBMat.GetSubMatrix(bMatOneNode, 0, actNode*nrDofs);

        helpMat = bMatOneNode * displDerivTransp;

        bMat.SetSubMatrix(helpMat, 0, actNode*nrDofs);
      }


    if (isaxi_)
      {
        const UInt spaceDim  = ptelem->GetDim();

        Vector<Double> shpFncAtIp;
        ptelem->GetShFncAtIp(shpFncAtIp, ip, it1_.GetElem() );

        Vector<Double> coordAtIP;
        coordAtIP = ptCoord * shpFncAtIp;

        Double  l33=0;  // for l33, see Bathe, page 552
        for (UInt actPos=0; actPos < shpFncAtIp.GetSize(); actPos++)
          l33 += elemDisp_[0][actPos] * shpFncAtIp[actPos] / coordAtIP[0];
        
        
        for (UInt actNode = 0; actNode < numFncs; actNode++)      
          {         
            //  (N_a/r) * (u_x/r)
            bMat[dimD - 1][actNode * spaceDim]     = l33 * shpFncAtIp[actNode] / coordAtIP[0];
            bMat[dimD - 1][actNode * spaceDim + 1] = 0;
          }
      }

  }
  

  // =============================================================================
  // 3D nonlinear mechanics (part with nonlinear B-mat)
  // =============================================================================


  nLinMech3dInt_BNonLin::nLinMech3dInt_BNonLin(BaseMaterial* matData) 
    : nLinElastInt(matData)
  {
    name_ = "nLinMech3dInt_BNonLin";
  }
 

  nLinMech3dInt_BNonLin::~nLinMech3dInt_BNonLin()
  {
  }


  // calculates the D-matrix of a 3d-problem 
  // (see mech3dInt::calcDMat)
  void nLinMech3dInt_BNonLin::calcDMat(Matrix<Double> & dMat)
  {
    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::REAL);
  }



  // =============================================================================
  // 3D nonlinear mechanics (material nonlinearity)
  // =============================================================================


  // class for mechanical material nonlinearities ...
 nLinMech3dInt_Material::nLinMech3dInt_Material(ApproxData *nlinFnc,
                                                BaseMaterial* matData) 
   : nLinElastInt(matData)
  {
    name_ = "nLinMech3dInt_Material";
    updateDMatInEveryIP_ = 1;
    isSolDependent_ = true;
    nLinFnc_ = nlinFnc;
  }
 

  nLinMech3dInt_Material::~nLinMech3dInt_Material()
  {
  }


  // calculates the D-matrix of a 3d-problem 
  // (see mech3dInt::calcDMat)
  void nLinMech3dInt_Material::calcDMat(Matrix<Double> & dMat, 
                                        UInt ip, 
                                        Matrix<Double> & ptCoord)
  {
    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::REAL);
   
    Matrix<Double>  displDeriv;  
    Vector<Double> stressVec;

    Vector<Double> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);

    Matrix<Double> linBMat;    
    calcBMat( linBMat, ip, ptCoord);

    std::string nonLinDepend;
    ptMaterial->GetScalar(nonLinDepend, NONLIN_DEPENDENCY);

    Integer nonLinCoeff;
    ptMaterial->GetScalar(nonLinCoeff, NONLIN_COEFFICIENT);

    std::string nonLinApproxType;
    ptMaterial->GetScalar(nonLinApproxType,NONLIN_APPROXIMATION_TYPE);


    if(nonLinDepend=="mechStress"){

      //Compute element stress
      
      Vector<Double> linStrain(linBMat.GetNumRows());
      linStrain.Init();
      stressVec.Init();
      linStrain = linBMat * displVec;
      stressVec = dMat * linStrain;
      
      Double stressSum=0.0;
      for(UInt i=0; i<stressVec.GetSize();i++)
        stressSum+=std::abs(stressVec[i]);

      Double c11, c13;
      
      if (nonLinApproxType=="polynomial"){
        if(nonLinCoeff==33)
          dMat[2][2]+=0.01*(stressSum + 0.02*stressSum*stressSum)*dMat[2][2];
        else if(nonLinCoeff==11){
          c11 = dMat[0][0]+0.01*(stressSum + 0.02*stressSum*stressSum)*dMat[0][0];
          dMat[0][0] = c11;
          dMat[1][1] = c11;
          dMat[5][5] = 0.5 * (c11 - dMat[0][1]);
        }
        else if(nonLinCoeff==13){
          c13 = dMat[0][2]+0.01*(stressSum + 0.02*stressSum*stressSum)*dMat[0][2];
          dMat[0][2] = c13;
          dMat[1][2] = c13;
          dMat[2][0] = c13;
          dMat[2][1] = c13;
        }
        else
          EXCEPTION("The nonlinear parameter dependency is not known");
      }
      else if (nonLinApproxType=="smoothSplines")
        {
          Double approxCoeff;

          approxCoeff = nLinFnc_->EvaluateFunc(stressSum);
          if(nonLinCoeff==33){
            dMat[2][2] = approxCoeff * 1.0e+11;
          }
        }
      else
        EXCEPTION("The nonlinear approximation type is not known");
    }
    else{
      std::cout<<"The data dependency you have chosen is " << nonLinDepend <<std::endl;
      std::cout<<"(If this is true, check if there is a blank in front of choice) " <<std::endl;
      EXCEPTION("Nonlinear dependency not implemented here");
    }
  }

  void nLinMech3dInt_Material::calcBMat( Matrix<Double> &bMat, UInt ip,
                                         Matrix<Double> &ptCoord ) {


    const UInt nrNodes  = ptelem->GetNumNodes();
    const UInt spaceDim = ptelem->GetDim();  
    const UInt nrDofs   = getNrDofs();  

    UInt actDim, actNode, j, k;
    
    
    bMat.Resize(getDimD(), nrNodes * nrDofs);
    bMat.Init();
    
    // local shape functions derived after global coords
    // (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;

    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem() );

    for(actDim=0; actDim < spaceDim; actDim++)
      for(actNode=0; actNode < nrNodes; actNode++)
        bMat[actDim][actNode * spaceDim + actDim] = xiDx[actNode][actDim];

    switch(spaceDim)
      {
      case 2:
        j = 1;
        k = 0;
        
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[spaceDim][actNode * spaceDim + 1] = xiDx[actNode][0];
            bMat[spaceDim][actNode * spaceDim]     = xiDx[actNode][1];
          }

        if (isaxi_)
          {
            UInt idxtheta = getDimD();
            Vector<Double> ShpFncAtIp;
            Vector<Double> CoordAtIP;

            if (isSetIntPoint_) 
              ptelem->GetShFnc(ShpFncAtIp,intPoint_,it1_.GetElem());
            else
              ptelem->GetShFncAtIp(ShpFncAtIp,ip,it1_.GetElem());

            CoordAtIP = ptCoord * ShpFncAtIp;

            for (actNode = 0; actNode < nrNodes; actNode++)          
              bMat[idxtheta-1][actNode * spaceDim] =
                ShpFncAtIp[actNode] / CoordAtIP[0];
          }

        break;


      case 3:
        UInt actDim=spaceDim;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][1];
          }

        actDim++;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][2];
            bMat[actDim][actNode * spaceDim + 2] = xiDx[actNode][0];
          }

        actDim++;
        for (actNode = 0; actNode < nrNodes; actNode++)
          {
            bMat[actDim][actNode * spaceDim]     = xiDx[actNode][1];
            bMat[actDim][actNode * spaceDim + 1] = xiDx[actNode][0];
          }
        break;
      }

    isSetIntPoint_ = false;
  }






  // =============================================================================
  // 3D nonlinear mechanics (part of Piola-Kirchhoff stress)
  // =============================================================================

  nLinMechInt_PiolaStress::nLinMechInt_PiolaStress(BaseMaterial* matData) 
    : nLinElastInt(matData)
  {

    name_ = "nLinMechInt_PiolaStress";
    updateDMatInEveryIP_ = 1;
    //    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMechInt_PiolaStress::~nLinMechInt_PiolaStress()
  {
  }






  /// calculates Piola-Kirchoff-stresses T (vector notation)
  // T = c . V with c the stiffness matrix and V the Green-Lagrange strain vector
  // V = B_lin * u  +  1/2 * B_nonLin * u
  // see Habil. M. Kaltenbacher Eq. (3.33)
  void nLinMechInt_PiolaStress::
  CalcStressVec(Vector<Double>& piolaStressVec, UInt ip, BaseFE *elem, 
                Matrix<Double> & ptCoord, Matrix<Double>& elemDisp )
  {

    // store information

    elemDisp_ = elemDisp;
    ptelem = elem;
    

    Matrix<Double> dMat;
    calcMaterialDMat(dMat);
  
    // convert displacement of all elem nodes into one vector: 
    // (uNode1X, uNode1Y, uNode2X, uNode2Y, ...)
    Vector<Double> displVec;
    elemDisp_.ConvertToVec_AppendCols(displVec);


    // linear differential operator B_lin
    Matrix<Double> linBMat;    
    calcLinBMat( linBMat, ip, ptelem, ptCoord);


    // nonlinear differential operator B_nonLin
    Matrix<Double> nLinBMat;    
    calcNonLinBMat( nLinBMat, ip, ptelem, ptCoord, elemDisp_);


    // special construction: direct addition of these two terms does not work, 
    // because of temporary vectors ... :O(    
    Vector<Double> part1;
    part1 = (linBMat * displVec );
    Vector<Double> part2;
    part2 = (nLinBMat * displVec);
    part2 *= 0.5;
    Vector<Double> nonLinStrain;
  
    nonLinStrain = part1 + part2;
    piolaStressVec = dMat * nonLinStrain;
  }




  void nLinMechInt_PiolaStress::
  calcMaterialDMat(Matrix<Double> & dMat)
  {

    //full matrix, axi as well as plane case overwrites this method

    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::REAL);
  }




  // returns linear B - matrix
  void nLinMechInt_PiolaStress::
  calcLinBMat(Matrix<Double> & bMat, UInt ip, BaseFE *elem, Matrix<Double> & ptCoord)
  {

    ptelem = elem;

    // dirty trick: dimension of d-matrix is set to 6 (as it should be in 3d mechanics)
    // this is done, because the special Piola-Kirchhoff-BDB uses a d-matrix of size 9
    // this means, the calculation of the "standard" (linear mech.) and the nonlinear b-matrix
    // would be of wrong dimension
    setPiolaDimD(  getMaterialDMatSize() );
  
    // linear differential operator B_lin
    linElastInt::calcBMat(bMat, ip, ptCoord);

    setPiolaDimD( getFullPiolaDMatSize() );
  }





  /// returns nonlinear B - matrix
  void nLinMechInt_PiolaStress::
  calcNonLinBMat(Matrix<Double> & bMat, UInt ip, BaseFE * elem, 
                 Matrix<Double> & ptCoord, Matrix<Double>& elemDisp)
  {

    ptelem = elem;
    elemDisp_ = elemDisp;

    // dirty trick: dimension of d-matrix is set to 6 (as it should be in 3d mechanics)
    // this is done, because the special Piola-Kirchhoff-BDB uses a d-matrix of size 9
    // this means, the calculation of the "standard" (linear mech.) and the nonlinear b-matrix
    // would be of wrong dimension
    setPiolaDimD(  getMaterialDMatSize() );

  
    // linear differential operator B_lin
    nLinElastInt::calcBMat(bMat, ip, ptCoord);
    setPiolaDimD( getFullPiolaDMatSize() );
  }

  

  /// conversion of stress vector to stress tensor
  void nLinMechInt_PiolaStress::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {

    UInt indexRow[] = {1, 2, 3, 2, 1, 1}; // first index of tensor notation
    UInt indexCol[] = {1, 2, 3, 3, 3, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs() );
    stressTensor.Init();
    
    // build symmetrical tensor
    for (UInt i=0; i<piolaStress.GetSize(); i++)
      {
        stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];      
        stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];      
      }

  }


  // calculates the D-matrix needed for the Piola-Kirchhoff-Stress part
  // This matrix is equal to a block diagonal matrix with Piola-Stresses in tensor 
  // notation used as diagonal blocks (nrDim times)
  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void nLinMechInt_PiolaStress::calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord)
  {

    const UInt dimD = getDimD();
    const UInt nrDofs = getNrDofs();

    Vector<Double> piolaStressVec;
    Matrix<Double> stressTensor;
    
    CalcStressVec(piolaStressVec, ip, ptelem, ptCoord, elemDisp_ );
    convertStressVecToTensor(stressTensor, piolaStressVec);

    dMat.Resize(dimD);
    dMat.Init();
    

    if (isaxi_)
      // the stressTensor has already the correct shape due to the 
      // method "convertStressVecToTensor"
      dMat = stressTensor;
    else
      for (UInt i=0; i<nrDofs; i++)
        dMat.SetSubMatrix(stressTensor, i*nrDofs, i*nrDofs);    
  }



  // returns B - matrix for piola stresses
  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void nLinMechInt_PiolaStress::
  calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
  {

    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrDofs   = getNrDofs();    


    bMat.Resize(getDimD(), numFncs * nrDofs);
    bMat.Init();

    
    // local shape functions derived after global coords 
    // (format: numFncs x nrDofs)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem() );


    for( UInt actNode=0; actNode < numFncs; actNode++ )
      for( UInt globPos=0; globPos < nrDofs; globPos++ )
        for( UInt actDof=0; actDof < nrDofs; actDof++ )
          bMat[globPos*nrDofs + actDof][actNode * nrDofs + globPos] = xiDx[actNode][actDof];      


    if (isaxi_)
      {
        const UInt spaceDim = ptelem->GetDim();
        Vector<Double> shpFncAtIp;
        Vector<Double> coordAtIP;
        
        ptelem->GetShFncAtIp(shpFncAtIp,ip, it1_.GetElem() );
        coordAtIP = ptCoord * shpFncAtIp;
        
        for ( UInt  actNode = 0; actNode < numFncs; actNode++)          
          bMat[getDimD() -1][actNode * spaceDim] = shpFncAtIp[actNode] / coordAtIP[0];
      }
  }


  // =============================================================================
  // class for nonlinear 3d piola stresses
  // =============================================================================


  nLinMech3dInt_PiolaStress::nLinMech3dInt_PiolaStress(BaseMaterial* matData) 
    :nLinMechInt_PiolaStress(matData)
  {

    name_ = "nLinMech3dInt_PiolaStress";
    updateDMatInEveryIP_ = 1;
    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMech3dInt_PiolaStress::~nLinMech3dInt_PiolaStress()
  {

  }


  // =============================================================================
  // class regarding mechanical prestresses
  // =============================================================================
 
  // base class for regarding prestress
  PreStressInt::PreStressInt(BaseMaterial* matData, Vector<Double> aPreStressVal) 
    : nLinMechInt_PiolaStress(matData), preStressVal_(aPreStressVal)
  {

    name_ = "PreStressInt";

    updateDMatInEveryIP_ = 1;
    setPiolaDimD( getFullPiolaDMatSize() );
  }


  PreStressInt::~PreStressInt()
  {

  }


  // calculates the D-matrix needed for regarding pre-stresses
  void PreStressInt::CalcStressVec(Vector<Double>& preStressVec, UInt ip, Matrix<Double> & ptCoord)
  {

    // ???????????????????????????????????????????????
    //    const UInt dimD = nLinMechInt_BNonLin::getDimD();
    const UInt dimD = getDimD();
    
    Vector<Double> preStrainVec(dimD);
    preStressVec.Resize(dimD);
    preStressVec.Init(0);

    Matrix<Double> matData;
    if (isaxi_) {
      //axi
      ptMaterial->GetTensor(matData,MECH_STIFFNESS_TENSOR,Global::REAL,AXI);
      preStrainVec[0]= preStressVal_[0] / matData[0][0];
      preStrainVec[1]= preStressVal_[1] / matData[1][1];
    }
    else if (getNrDofs() == 3) {
      //3D
      ptMaterial->GetTensor(matData,MECH_STIFFNESS_TENSOR,Global::REAL);
      preStrainVec[0]= preStressVal_[0] / matData[0][0];
      preStrainVec[1]= preStressVal_[1] / matData[1][1];
      preStrainVec[2]= preStressVal_[2] / matData[2][2];
    }
    else {
      //2D-plane strain
      ptMaterial->GetTensor(matData,MECH_STIFFNESS_TENSOR,Global::REAL, PLANE_STRAIN);
      preStrainVec[0]= preStressVal_[0] / matData[0][0];
      preStrainVec[1]= preStressVal_[1] / matData[1][1];
    }
   
    preStressVec = matData * preStrainVec;
  }




  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void PreStressInt::calcDMat(Matrix<Double> & dMat, UInt ip, Matrix<Double> & ptCoord)
  {

    const UInt dimD = getFullPiolaDMatSize();
    const UInt nrDofs = getNrDofs();

    Vector<Double> piolaStressVec;
    Matrix<Double> stressTensor;
    
    CalcStressVec(piolaStressVec, ip, ptCoord );
    convertStressVecToTensor(stressTensor, piolaStressVec);

    dMat.Resize(dimD);
    dMat.Init();
    

    if (isaxi_)
      // the stressTensor has already the correct shape due to the 
      // method "convertStressVecToTensor"
      dMat = stressTensor;
    else
      for (UInt i=0; i<nrDofs; i++)
        dMat.SetSubMatrix(stressTensor, i*nrDofs, i*nrDofs);    
  }


  // returns B - matrix for piola stresses
  // (see e.g. Bathe: "Finite Element Procedures" p. 556)
  void PreStressInt::calcBMat(Matrix<Double> & bMat, UInt ip, Matrix<Double> & ptCoord)
  {

    UInt numFncs = ptelem->GetNumFncs( ansatzFct1_ );
    ptelem->SetAnsatzFct( ansatzFct1_ );
    const UInt nrDofs   = getNrDofs();    

    UInt dimD = getFullPiolaDMatSize();
    bMat.Resize(dimD, numFncs * nrDofs);
    bMat.Init();

    
    // local shape functions derived after global coords 
    // (format: numFncs x nrDofs)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem() );


    for( UInt actNode=0; actNode < numFncs; actNode++ )
      for( UInt globPos=0; globPos < nrDofs; globPos++ )
        for( UInt actDof=0; actDof < nrDofs; actDof++ )
          bMat[globPos*nrDofs + actDof][actNode * nrDofs + globPos] = xiDx[actNode][actDof];      


    if (isaxi_)
      {
        const UInt spaceDim = ptelem->GetDim();
        Vector<Double> shpFncAtIp;
        Vector<Double> coordAtIP;
        
        ptelem->GetShFncAtIp(shpFncAtIp,ip, it1_.GetElem() );
        coordAtIP = ptCoord * shpFncAtIp;
        
        for ( UInt actNode = 0; actNode < numFncs; actNode++)          
          bMat[getDimD() -1][actNode * spaceDim] = shpFncAtIp[actNode] / coordAtIP[0];
      }
  }


 
  // class for regarding 3d prestress
  PreStressInt3D::PreStressInt3D(BaseMaterial* matData, Vector<Double> aPreStressVal)
    : PreStressInt(matData, aPreStressVal)
  {

    name_ = "PreStressInt3D";
  }


  PreStressInt3D::~PreStressInt3D()
  {

  }

  /// conversion of stress vector to stress tensor in 3D
  void PreStressInt3D::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {

    UInt indexRow[] = {1, 2, 3, 2, 1, 1}; // first index of tensor notation
    UInt indexCol[] = {1, 2, 3, 3, 3, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs());
    stressTensor.Init();
    // build symmetrical tensor
    for (UInt i=0; i<piolaStress.GetSize(); i++)
      {
        stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];      
        stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];      
      }
  }


 
  // class for regarding 2d prestress in plane strain case
  PreStressIntPlaneStrain::PreStressIntPlaneStrain(BaseMaterial* matData, 
                                                   Vector<Double> aPreStressVal)
    : PreStressInt(matData, aPreStressVal)
  {

    name_ = "PreStressIntPlaneStrain";
  }


  PreStressIntPlaneStrain::~PreStressIntPlaneStrain()
  {

  }

  /// conversion of stress vector to stress tensor
  void PreStressIntPlaneStrain::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {

    UInt indexRow[] = {1, 2, 1}; // first index of tensor notation
    UInt indexCol[] = {1, 2, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs());
    stressTensor.Init();
    
    // build symmetrical tensor
    for (UInt i=0; i<piolaStress.GetSize(); i++)
      {
        stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];      
        stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];      
      }
  }


  
  // class for regarding 2d axi prestress 
  PreStressIntAxi::PreStressIntAxi(BaseMaterial* matData, 
                                   Vector<Double> aPreStressVal)
    : PreStressInt(matData, aPreStressVal)
  {

    name_ = "PreStressIntAxi";
    isaxi_ = true;
  }


  PreStressIntAxi::~PreStressIntAxi()
  {

  }

  /// conversion of stress vector to stress tensor
  void PreStressIntAxi::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {

    // indizes see Bathe: "Finite Element Procedures", Sec. 6.3, 
    // page 553: "2. Piola-Kirchhoff stress matrix & vector"
    const UInt indexSize = 9;
    UInt indexRow[]   = {1, 2, 1, 2, 3, 4, 3, 4, 5}; // first index of tensor notation
    UInt indexCol[]   = {1, 2, 2, 1, 3, 4, 4, 3, 5}; // second index of tensor notation
    UInt indexPiola[] = {1, 2, 3, 3, 1, 2, 3, 3, 4}; // index in piola stress vec

    const UInt axiTensSize = 5;
    
    stressTensor.Resize( axiTensSize);
    stressTensor.Init();
    
    // build symmetrical tensor
    for (UInt i=0; i<indexSize; i++)
      {
        stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[ indexPiola[i] -1 ];     
        stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[ indexPiola[i] -1 ];     
      }

  }

  // =============================================================================
  // 2D nonlinear plane strain mechanics
  // =============================================================================

  
  // calculated the D-matrix for the plain strain state
  void  nLinMechPlaneStrainInt_BNonLin::calcDMat(Matrix<Double> & dMat)
  {

    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::REAL,PLANE_STRAIN);
  }

  




  void nLinMechPlaneStrainInt_PiolaStress::calcMaterialDMat(Matrix<Double> & dMat)
  {

    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::REAL,PLANE_STRAIN);
  }




  /// conversion of stress vector to stress tensor
  void nLinMechPlaneStrainInt_PiolaStress::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {

    UInt indexRow[] = {1, 2, 1}; // first index of tensor notation
    UInt indexCol[] = {1, 2, 2}; // second index of tensor notation

    stressTensor.Resize( getNrDofs());
    stressTensor.Init();
    
    // build symmetrical tensor
    for (UInt i=0; i<piolaStress.GetSize(); i++)
      {
        stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[i];      
        stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[i];      
      }

  }


  // =============================================================================
  // 2D nonlinear axisymmetric mechanics
  // =============================================================================

  
  // calculated the D-matrix for the axi state
  void  nLinMechAxiInt_BNonLin::calcDMat(Matrix<Double> & dMat)
  {

     ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::REAL,AXI);
  }



  // nonlinear calculation of elasticity for axi 
  nLinMechAxiInt_PiolaStress::nLinMechAxiInt_PiolaStress(BaseMaterial* matData) 
    : nLinMechInt_PiolaStress(matData)

  {

    name_ = "nLinMechAxiInt_PiolaStress";
    isaxi_ = true;
    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMechAxiInt_PiolaStress::~nLinMechAxiInt_PiolaStress()
  {

  }




  void nLinMechAxiInt_PiolaStress::calcMaterialDMat(Matrix<Double> & dMat)
  {

    ptMaterial->GetTensor(dMat,MECH_STIFFNESS_TENSOR,Global::REAL,AXI);
 
  }


  /// conversion of stress vector to stress tensor
  void nLinMechAxiInt_PiolaStress::
  convertStressVecToTensor(Matrix<Double>& stressTensor, Vector<Double>& piolaStress)
  {

    // indizes see Bathe: "Finite Element Procedures", Sec. 6.3, 
    // page 553: "2. Piola-Kirchhoff stress matrix & vector"
    const UInt indexSize = 9;
    UInt indexRow[]   = {1, 2, 1, 2, 3, 4, 3, 4, 5}; // first index of tensor notation
    UInt indexCol[]   = {1, 2, 2, 1, 3, 4, 4, 3, 5}; // second index of tensor notation
    UInt indexPiola[] = {1, 2, 3, 3, 1, 2, 3, 3, 4}; // index in piola stress vec

    const UInt axiTensSize = 5;
    
    stressTensor.Resize( axiTensSize);
    stressTensor.Init();
    
    // build symmetrical tensor
    for (UInt i=0; i<indexSize; i++)
      {
        stressTensor[ indexRow[i] -1 ][ indexCol[i] -1 ] = piolaStress[ indexPiola[i] -1 ];     
        stressTensor[ indexCol[i] -1 ][ indexRow[i] -1 ] = piolaStress[ indexPiola[i] -1 ];     
      }

  }


  // ===================================================================================
  // nonlinear calculation of elasticity in plane strain state 
  // ===================================================================================

  nLinMechPlaneStrainInt_BNonLin::nLinMechPlaneStrainInt_BNonLin(BaseMaterial* matData) 
    : nLinElastInt(matData)
  {
    name_ = "nLinMechPlaneStrainInt_BNonLin";

  }
 

  nLinMechPlaneStrainInt_BNonLin::~nLinMechPlaneStrainInt_BNonLin()
  {

  }

  
  // nonlinear calculation of elasticity in 3d
  nLinMechPlaneStrainInt_PiolaStress::
  nLinMechPlaneStrainInt_PiolaStress(BaseMaterial* matData) 
    : nLinMechInt_PiolaStress(matData)

  {

    name_ = "nLinMechPlaneStrainInt_PiolaStress";
    setPiolaDimD( getFullPiolaDMatSize() );
  }
 

  nLinMechPlaneStrainInt_PiolaStress::~nLinMechPlaneStrainInt_PiolaStress()
  {
  }





  // ===================================================================================
  // axisymmetric nonlinear calculation of elasticity
  // ===================================================================================

  nLinMechAxiInt_BNonLin::nLinMechAxiInt_BNonLin(BaseMaterial* matData) 
    : nLinElastInt(matData)
  {

    name_ = "nLinMechAxiInt_BNonLin";
    isaxi_ = true;
  }
 

  nLinMechAxiInt_BNonLin::~nLinMechAxiInt_BNonLin()
  {

  }


} // end namespace CoupledField
