// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "General/environment.hh"
#include "nLinElecMaterial.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField {
  
  
  // class for electrical material nonlinearities ...
  nLinElec3dInt_Material::nLinElec3dInt_Material(ApproxData *nlinFnc, 
                                                 BaseMaterial* matData,
                                                 SubTensorType type) 
    : BDBInt(matData, type)
  {
    
    name_ = "nLinElec3dInt_Material";

    updateDMatInEveryIP_ = 1;
    isSolDependent_ = true;

    nLinFnc_ = nlinFnc;

    if ( type == FULL ) {
      dim_ = 3;
    }
    else {
      dim_ = 2;
    }
    
    if ( type == AXI ) {
      isaxi_     = true;
    }
    
  }
  
  nLinElec3dInt_Material::~nLinElec3dInt_Material()
  {
  }


  // ====================
  //   set 
  // ====================
  void nLinElec3dInt_Material::Set4NonLinMaterial(Grid* ptGrid, 
                                                  StdPDE* ptPDE,
                                                  shared_ptr<EqnMap> eqnMap,
                                                  shared_ptr<ResultInfo> result) 
  {
    
    EfieldOp_ =  new GradientFieldOp<Double>(ptGrid, ptPDE, eqnMap, *sol_,
                                             result->fctType, isaxi_, coordUpdate_);
  }


  // ============
  //   calcBMat
  // ============
  void nLinElec3dInt_Material::CalcBMat( Matrix<Double> &bMat, UInt ip,
			     const Matrix<Double> &ptCoord ) {


    // Obtain info on number of element's nodes
    const UInt numNodes = ptelem->GetNumNodes();


    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( dim_, numNodes );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord, it1_.GetElem() );

    if ( subTensorType_ == FULL ) {
      // The matrix bMat can be seen as a 3 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the BDB product evaluated at the k-th node of the finite
      // element. 
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
        bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    else  {
      // The matrix bMat can be seen as a 1 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the ADB product evaluated at the k-th node of the finite
      // element. We assume that the first coordinate equals y and the
      // second z.
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
      }
    }
    
  }

  // ====================
  //   calcElementMatrix
  // ====================
  void nLinElec3dInt_Material::CalcElementMatrix( Matrix<Double>& elemMat,
                                                  EntityIterator& ent1, 
                                                  EntityIterator& ent2 ) {

    // get displacements of element
    sol_->GetElemSolution( elemPot_, ent1 );

    Vector<Double> LCoord, Efield;

    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 ); 

    ent1_ = ent1;
    
     // call method of base class
    BDBInt::CalcElementMatrix( elemMat, ent1, ent2 );
  }

 

  // calc nonlinear DMAt
  void nLinElec3dInt_Material::calcDMat(Matrix<Double> & dMat, 
                                        UInt ip, 
                                        Matrix<Double> & ptCoord)
  {
    ptMaterial->GetTensor(dMat,ELEC_PERMITTIVITY,matDataType_,subTensorType_);
   
    Matrix<Double> xiDx;
    ptelem->GetGlobDerivShFncAtIp(xiDx, ip, ptCoord, it1_.GetElem());

    Vector<Double> LCoord, Efield;

    //compute electric field in the midpoint of element
    ptelem->GetCoordMidPoint(LCoord);
    
    EfieldOp_->CalcElemGradField( Efield, ent1_, LCoord, 1 );

    std::string nonLinDepend;
    ptMaterial->GetScalar(nonLinDepend, NONLIN_DEPENDENCY);
    
    Integer nonLinCoeff;
    ptMaterial->GetScalar(nonLinCoeff, NONLIN_COEFFICIENT);
    
    std::string nonLinApproxType;
    ptMaterial->GetScalar(nonLinApproxType,NONLIN_APPROXIMATION_TYPE);
    
    
    if(nonLinDepend=="elecField"){
      
      //       //Compute element stress
      
      
       Double eFieldSum=0.0;
       for(UInt i=0; i<Efield.GetSize();i++)
         eFieldSum+=std::abs(Efield[i]);

       Double eps11;
      
        if (nonLinApproxType=="polynomial"){
          if(nonLinCoeff==33)
            dMat[2][2]+=0.01*(eFieldSum + 0.02*eFieldSum*eFieldSum)*dMat[2][2];
          else if(nonLinCoeff==11){
            eps11 = dMat[0][0]+0.01*(eFieldSum + 0.02*eFieldSum*eFieldSum)*dMat[0][0];
            dMat[0][0] = eps11;
            dMat[1][1] = eps11;
          }
          
          else
            EXCEPTION("The nonlinear parameter dependency is not known");
        }
        else if (nonLinApproxType=="smoothSplines")
         {
           Double approxCoeff;

           approxCoeff = nLinFnc_->EvaluateFunc(eFieldSum);
        
           if(nonLinCoeff==33){
             dMat[2][2] = approxCoeff;
           }
         }
       else
         EXCEPTION("The nonlinear approximation type is not known");
     }
     else{
       std::cout<<"The data dependency you have chosen is " << nonLinDepend <<std::endl;
       std::cout<<"(If this is true, check if there is a blank in font of choice) " <<std::endl;
       EXCEPTION("Nonlinear dependency not implemented here");
     }

  }

}







 
