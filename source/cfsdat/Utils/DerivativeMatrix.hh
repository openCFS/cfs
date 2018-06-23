// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     DerivativeMatrix.hh
 *       \brief    Struct for a sparse derivatve matrix in R3, hence each entry is a vector of length 3
 * Furthermore some optimized functions to compute sveral types of derivatives
 *
 *       \date     Apr 4, 2018
 *       \author   matutz
 */
//================================================================================================



#ifndef DERIVATIVEMATRIX_H
#define	DERIVATIVEMATRIX_H

#include "Utils/StdVector.hh"
#include "General/Environment.hh"

namespace CFSDat {
  
struct DerivativeMatrix {
  
  //! number of output values
  UInt outCount;
  
  //! number of input values
  UInt inCount;
  
  //! indices to find the vectors/coefficients for the derivative computation
  StdVector<UInt> outInIndex;

  //! indices of the input values
  StdVector<UInt> outIn;

  //! vectorial factors for the input values
  StdVector<Double> outInVector;
  
  //! Clear the output values
  template<typename T, int S>
  void ClearOut(T* out) {
    const UInt outSize = outCount * S;
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt i = 0; i < outSize; i++) {
      out[i] = 0;
    }
  }
  
  //! Calculate divergence of a vector
  template<typename T>
  void VectorDivergence(T* out, T* in, bool add) {
    if (!add) {
      ClearOut<T,1>(out);
    }
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt iOut = 0; iOut < outCount; iOut++) {
      T& outSum = out[iOut];
      UInt iOutInIndex = outInIndex[iOut];
      UInt iOutInCount = outInIndex[iOut + 1] - iOutInIndex;
      UInt* iOutIn = &outIn[iOutInIndex];
      Double* iOutInVector = &outInVector[iOutInIndex * 3];
      for (UInt iIn = 0; iIn < iOutInCount; iIn++) {
        T* inValue = &in[iOutIn[iIn] * 3];
        outSum += inValue[0] * iOutInVector[iIn * 3 + 0]
                + inValue[1] * iOutInVector[iIn * 3 + 1]
                + inValue[2] * iOutInVector[iIn * 3 + 2];
      }
    }
  }

  //! Calculate divergence of a tensor
  template<typename T>
  void TensorDivergence(T* out, T* in, bool add) {
    if (!add) {
      ClearOut<T,3>(out);
    }
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt iOut = 0; iOut < outCount; iOut++) {
      T* outSum = &out[iOut * 3];
      UInt iOutInIndex = outInIndex[iOut];
      UInt iOutInCount = outInIndex[iOut + 1] - iOutInIndex;
      UInt* iOutIn = &outIn[iOutInIndex];
      Double* iOutInVector = &outInVector[iOutInIndex * 3];
      for (UInt iIn = 0; iIn < iOutInCount; iIn++) {
        T* inValue = &in[iOutIn[iIn] * 9];
        outSum[0] += inValue[0] * iOutInVector[iIn * 3 + 0]
                   + inValue[1] * iOutInVector[iIn * 3 + 1]
                   + inValue[2] * iOutInVector[iIn * 3 + 2];
        outSum[1] += inValue[3] * iOutInVector[iIn * 3 + 0]
                   + inValue[4] * iOutInVector[iIn * 3 + 1]
                   + inValue[5] * iOutInVector[iIn * 3 + 2];
        outSum[2] += inValue[6] * iOutInVector[iIn * 3 + 0]
                   + inValue[7] * iOutInVector[iIn * 3 + 1]
                   + inValue[8] * iOutInVector[iIn * 3 + 2];
      }
    }
  }

  //! Calculate gradient of a scalar
  template<typename T>
  void ScalarGradient(T* out, T* in, bool add) {
    if (!add) {
      ClearOut<T,3>(out);
    }
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt iOut = 0; iOut < outCount; iOut++) {
      T* outSum = &out[iOut * 3];
      UInt iOutInIndex = outInIndex[iOut];
      UInt iOutInCount = outInIndex[iOut + 1] - iOutInIndex;
      UInt* iOutIn = &outIn[iOutInIndex];
      Double* iOutInVector = &outInVector[iOutInIndex * 3];
      for (UInt iIn = 0; iIn < iOutInCount; iIn++) {
        T inValue = in[iOutIn[iIn]];
        outSum[0] += iOutInVector[iIn * 3 + 0] * inValue;
        outSum[1] += iOutInVector[iIn * 3 + 1] * inValue;
        outSum[2] += iOutInVector[iIn * 3 + 2] * inValue;
      }
    }
  }
      
  //! Calculate gradient of a vector
  template<typename T>
  void VectorGradient(T* out, T* in, bool add) {
    if (!add) {
      ClearOut<T,9>(out);
    }
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt iOut = 0; iOut < outCount; iOut++) {
      T* outSum = &out[iOut * 9];
      UInt iOutInIndex = outInIndex[iOut];
      UInt iOutInCount = outInIndex[iOut + 1] - iOutInIndex;
      UInt* iOutIn = &outIn[iOutInIndex];
      Double* iOutInVector = &outInVector[iOutInIndex * 3];
      for (UInt iIn = 0; iIn < iOutInCount; iIn++) {
        T* inValue = &in[iOutIn[iIn] * 3];
        outSum[0] += inValue[0] * iOutInVector[iIn * 3 + 0];
        outSum[1] += inValue[0] * iOutInVector[iIn * 3 + 1];
        outSum[2] += inValue[0] * iOutInVector[iIn * 3 + 2];
        outSum[3] += inValue[1] * iOutInVector[iIn * 3 + 0];
        outSum[4] += inValue[1] * iOutInVector[iIn * 3 + 1];
        outSum[5] += inValue[1] * iOutInVector[iIn * 3 + 2];
        outSum[6] += inValue[2] * iOutInVector[iIn * 3 + 0];
        outSum[7] += inValue[2] * iOutInVector[iIn * 3 + 1];
        outSum[8] += inValue[2] * iOutInVector[iIn * 3 + 2];
      }
    }
  }
      
  //! Calculate curl of a vector
  template<typename T>
  void VectorCurl(T* out, T* in, bool add) {
    if (!add) {
      ClearOut<T,3>(out);
    }
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt iOut = 0; iOut < outCount; iOut++) {
      T* outSum = &out[iOut * 3];
      UInt iOutInIndex = outInIndex[iOut];
      UInt iOutInCount = outInIndex[iOut + 1] - iOutInIndex;
      UInt* iOutIn = &outIn[iOutInIndex];
      Double* iOutInVector = &outInVector[iOutInIndex * 3];
      for (UInt iIn = 0; iIn < iOutInCount; iIn++) {
        T* inValue = &in[iOutIn[iIn] * 3];
        outSum[0] += inValue[2] * iOutInVector[iIn * 3 + 1]
                   - inValue[1] * iOutInVector[iIn * 3 + 2];
        outSum[1] += inValue[0] * iOutInVector[iIn * 3 + 2]
                   - inValue[2] * iOutInVector[iIn * 3 + 0];
        outSum[2] += inValue[1] * iOutInVector[iIn * 3 + 0]
                   - inValue[0] * iOutInVector[iIn * 3 + 1];
      }
    }
  }

  
};

}

#endif	/* DERIVATIVEMATRIX_H */

