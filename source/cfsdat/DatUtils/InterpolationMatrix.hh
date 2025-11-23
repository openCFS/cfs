// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     InterpolationMatrix.hh
 *       \brief    Struct for a sparse interpolation matrix
 *
 *       \date     Apr 4, 2018
 *       \author   matutz
 */
//================================================================================================



#ifndef INTERPOLATIONMATRIX_H
#define	INTERPOLATIONMATRIX_H

#include "Utils/StdVector.hh"
#include "General/Environment.hh"

namespace CFSDat {
  
struct InterpolationMatrix {
  
  //! number of output values
  UInt outCount;

  //! number of input values
  UInt inCount;

  //! indices to find the coefficients for the interpolation
  StdVector<UInt> outInIndex;

  //! indices of the input values
  StdVector<UInt> outIn;

  //! vectorial factors for the input values
  StdVector<Double> outInFactor;

  //! Interpolate some values
  template<typename T>
  void Interpolate(T* out, T* in, UInt numDofs) {
    if (numDofs == 1) {
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (UInt iOut = 0; iOut < outCount; iOut++) {
        T sum(0.0);
        UInt iOutInIndex = outInIndex[iOut];
        UInt iOutInCount = outInIndex[iOut + 1] - iOutInIndex;
        if (iOutInCount > 0) {
          UInt* iOutIn = &outIn[iOutInIndex];
          Double* iOutInFactor = &outInFactor[iOutInIndex];
          for (UInt iIn = 0; iIn < iOutInCount; iIn++) {
            sum += in[iOutIn[iIn]] * iOutInFactor[iIn];
          }
        }
        out[iOut] = sum;
      }
    } else {
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (UInt iOut = 0; iOut < outCount; iOut++) {
        T* sum = &out[iOut * numDofs];
        for (UInt d = 0; d < numDofs; d++) {
          sum[d] = 0.0;
        }
        UInt iOutInIndex = outInIndex[iOut];
        UInt iOutInCount = outInIndex[iOut + 1] - iOutInIndex;
        if (iOutInCount > 0) {
          UInt* iOutIn = &outIn[iOutInIndex];
          Double* iOutInFactor = &outInFactor[iOutInIndex];
          for (UInt iIn = 0; iIn < iOutInCount; iIn++) {
            Double iFactor = iOutInFactor[iIn];
            T* iInVector = &in[iOutIn[iIn] * numDofs];
            for (UInt d = 0; d < numDofs; d++) {
              sum[d] += iInVector[d] * iFactor;
            }
          }
        }
      }
    }
  }
  
};

}

#endif	/* INTERPOLATIONMATRIX_H */

