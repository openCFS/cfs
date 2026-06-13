// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AeroacousticBase.cc
 *       \brief    <Description>
 *
 *       \date     Jan 7, 2017
 *       \author   kroppert
 */
//================================================================================================

#include "AeroacousticBase.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "FeBasis/H1/H1Elems.hh"


namespace CFSDat{

AeroacousticBase::AeroacousticBase(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){


}

AeroacousticBase::~AeroacousticBase(){

}



void AeroacousticBase::OmegaVectorProductU(const Vector<Double>& inVecVel,
                                          const Vector<Double>& inVecOmega,
                                          Vector<Double>& retVec,
                                          const UInt& dim){

  retVec.Resize(inVecVel.GetSize());
  CF::Matrix<Double> vel , omega;
  vel.Resize(inVecVel.GetSize() / dim, dim);
  vel.Init();

if(dim == 3){
  omega.Resize(inVecVel.GetSize() / dim, dim);
  omega.Init();
//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(CF::UInt j = 0; j < inVecVel.GetSize() / dim; ++j) {
    for(CF::UInt k = 0; k < dim; ++k){
      vel[j][k] = inVecVel[dim * j + k];
      omega[j][k] = inVecOmega[dim * j + k];
    }
   }
  for(CF::UInt j = 0; j < inVecVel.GetSize() / dim; ++j) {
    retVec[j * dim] =     vel[j][2] * omega[j][1] - vel[j][1] * omega[j][2];
    retVec[j * dim + 1] = vel[j][0] * omega[j][2] - vel[j][2] * omega[j][0];
    retVec[j * dim + 2] = vel[j][1] * omega[j][0] - vel[j][0] * omega[j][1];
  }
}else{
  // 2D case
  omega.Resize(inVecVel.GetSize() / dim, 1);
  omega.Init();
//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(CF::UInt j = 0; j < inVecVel.GetSize() / dim; ++j) {
    for(CF::UInt k = 0; k < dim; ++k){
      vel[j][k] = inVecVel[dim * j + k];
    }
    omega[j][0] = inVecOmega[j];
   }
//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(CF::UInt j = 0; j < inVecVel.GetSize() / dim; ++j) {
    retVec[j * dim] =     -vel[j][1] * omega[j][0];
    retVec[j * dim + 1] =  vel[j][0] * omega[j][0];
   }
}

}




void AeroacousticBase::ScalarProduct(Vector<Double>& retVec,
                                    const Vector<Double>& inVec1,
                                    const Vector<Double>& inVec2,
                                    const UInt& numEquPerEnt,
                                    const Double& scalarFactor){

  if(inVec1.GetSize() != inVec2.GetSize()) EXCEPTION("ScalarProduct: vectors don't have the same size!!!");

  retVec.Resize(inVec1.GetSize() / numEquPerEnt);
  retVec.Init();

  UInt maxNumTrgEntities = retVec.GetSize();

  for(UInt i = 0; i < maxNumTrgEntities; ++i){
    Double t = 0;
    for(UInt k = 0; k < numEquPerEnt; ++k){
      t += inVec1[numEquPerEnt * i + k] * inVec2[numEquPerEnt * i + k];
    }
    retVec[i] = t * scalarFactor;
  }
}

void AeroacousticBase::TensorProduct(Vector<Double>& retVec,
                                    const Vector<Double>& tensor1,
                                    const Vector<Double>& tensor2,
                                    const UInt& numEquPerEnt,
                                    const Double& scalarFactor,
                                    const Vector<Double>& scalar){

  if(tensor1.GetSize() != tensor2.GetSize()) EXCEPTION("TensorProduct: input tensors don't have the same size!!!");
//  if(tensor1.GetSize()/numEquPerEnt != scalar.GetSize()) EXCEPTION("TensorProduct: input tensors and scalar dimension don't match!!!");

  int a = 1;
  if (numEquPerEnt==3) a = 3;

  retVec.Resize(tensor1.GetSize() * (numEquPerEnt + a) / numEquPerEnt);
  retVec.Init();

  //We use Voigt notation to number the tensor
  for(UInt i = 0; i < scalar.GetSize(); ++i){
      if (a == 1) {
        retVec[(numEquPerEnt + a)*i + 0] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 0];
        retVec[(numEquPerEnt + a)*i + 1] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 1] * tensor2[numEquPerEnt * i + 1];
        retVec[(numEquPerEnt + a)*i + 2] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 1];
      }
      else if  (a == 3){
        retVec[(numEquPerEnt + a)*i + 0] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 0];
        retVec[(numEquPerEnt + a)*i + 1] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 1] * tensor2[numEquPerEnt * i + 1];
        retVec[(numEquPerEnt + a)*i + 2] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 2] * tensor2[numEquPerEnt * i + 2];

        retVec[(numEquPerEnt + a)*i + 3] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 1] * tensor2[numEquPerEnt * i + 2];

        retVec[(numEquPerEnt + a)*i + 4] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 2];
        retVec[(numEquPerEnt + a)*i + 5] = scalarFactor * scalar[i] * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 1];


      }

  }
}

void AeroacousticBase::TensorProduct(Vector<Double>& retVec,
                                    const Vector<Double>& tensor1,
                                    const Vector<Double>& tensor2,
                                    const UInt& numEquPerEnt,
                                    const Double& scalarFactor){

  if(tensor1.GetSize() != tensor2.GetSize()) EXCEPTION("TensorProduct: input tensors don't have the same size!!!");
//  if(tensor1.GetSize()/numEquPerEnt != scalar.GetSize()) EXCEPTION("TensorProduct: input tensors and scalar dimension don't match!!!");

  int a = 1;
  if (numEquPerEnt==3) a = 3;

  retVec.Resize(tensor1.GetSize() * (numEquPerEnt + a) / numEquPerEnt);
  retVec.Init();

  UInt iLoop = tensor1.GetSize()/numEquPerEnt;

  //We use Voigt notation to number the tensor
  for(UInt i = 0; i < iLoop; ++i){
      if (a == 1) {
        retVec[(numEquPerEnt + a)*i + 0] = scalarFactor * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 0];
        retVec[(numEquPerEnt + a)*i + 1] = scalarFactor * tensor1[numEquPerEnt * i + 1] * tensor2[numEquPerEnt * i + 1];
        retVec[(numEquPerEnt + a)*i + 2] = scalarFactor * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 1];
      }
      else if  (a == 3){
        retVec[(numEquPerEnt + a)*i + 0] = scalarFactor * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 0];
        retVec[(numEquPerEnt + a)*i + 1] = scalarFactor * tensor1[numEquPerEnt * i + 1] * tensor2[numEquPerEnt * i + 1];
        retVec[(numEquPerEnt + a)*i + 2] = scalarFactor * tensor1[numEquPerEnt * i + 2] * tensor2[numEquPerEnt * i + 2];

        retVec[(numEquPerEnt + a)*i + 3] = scalarFactor * tensor1[numEquPerEnt * i + 1] * tensor2[numEquPerEnt * i + 2];

        retVec[(numEquPerEnt + a)*i + 4] = scalarFactor * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 2];
        retVec[(numEquPerEnt + a)*i + 5] = scalarFactor * tensor1[numEquPerEnt * i + 0] * tensor2[numEquPerEnt * i + 1];


      }

  }
}



void AeroacousticBase::Node2Cell(Vector<Double>& returnVec,
    const Vector<Double>& inVec,
    const UInt& numEquPerEnt,
    const StdVector< StdVector<CF::UInt> >& targetSourceIndex,
    const UInt& maxNumTrgEntities,
    const UInt& gridDim){

  returnVec.Resize(maxNumTrgEntities * numEquPerEnt);
  returnVec.Init();


//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(UInt i = 0; i < maxNumTrgEntities; ++i){
    UInt patchSize = targetSourceIndex[i].GetSize();
    const StdVector<UInt>& sM = targetSourceIndex[i];
    CF::UInt targetIndex = i * numEquPerEnt;
    for (CF::UInt j = 0; j < patchSize; j++) {
      CF::UInt sourceIndex = sM[j];
      for(UInt k = 0; k < numEquPerEnt; ++k){
        returnVec[targetIndex + k] += inVec[numEquPerEnt * (sourceIndex-1) + k] / patchSize;
      }
    }
  }
}


Vector<Double> AeroacousticBase::ScalarToTwoD(const Vector<Double>& inVec){
  Vector<Double> r;
  r.Resize(2 * inVec.GetSize());
  for(UInt i = 0; i < inVec.GetSize(); ++i){
    r[2 * i] = inVec[i];
    r[2 * i + 1] = 0.0;
  }
  return r;
}

Vector<Double> AeroacousticBase::TwoDToScalar(const Vector<Double>& inVec){
  Vector<Double> r;
  r.Resize(inVec.GetSize() / 2);
  for(UInt i = 0; i < r.GetSize(); ++i){
    r[i] = inVec[2 * i];
  }
  return r;
}


}

