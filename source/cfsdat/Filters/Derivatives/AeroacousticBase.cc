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

AeroacousticBase::AeroacousticBase(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
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
//#pragma omp parallel for num_threads(NUM_CFS_THREADS)
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
//#pragma omp parallel for num_threads(NUM_CFS_THREADS)
  for(CF::UInt j = 0; j < inVecVel.GetSize() / dim; ++j) {
    for(CF::UInt k = 0; k < dim; ++k){
      vel[j][k] = inVecVel[dim * j + k];
    }
    omega[j][0] = inVecOmega[j];
   }
//#pragma omp parallel for num_threads(NUM_CFS_THREADS)
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




void AeroacousticBase::Node2Cell(Vector<Double>& returnVec,
      const Vector<Double>& inVec,
      const UInt& numEquPerEnt,
      const StdVector<CF::UInt>& targetSource,
      const StdVector<CF::UInt>& targetSourceIndex,
      const UInt& numNN,
      const UInt& maxNumTrgEntities,
      const UInt& gridDim){

    returnVec.Resize(maxNumTrgEntities * numEquPerEnt);
    returnVec.Init();

//#pragma omp parallel for num_threads(NUM_CFS_THREADS)
    for(UInt i = 0; i < maxNumTrgEntities; ++i){
      const CF::UInt jEnd = targetSourceIndex[i + 1];
      CF::UInt targetIndex = i * numEquPerEnt;
      for (CF::UInt j = targetSourceIndex[i]; j < jEnd; j++) {
        CF::UInt sourceIndex = targetSource[j];
        for(UInt k = 0; k < numEquPerEnt; ++k){
          returnVec[targetIndex + k] += inVec[numEquPerEnt * sourceIndex + k] / numNN;
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

