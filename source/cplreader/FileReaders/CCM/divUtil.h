/* 
 * File:   divUtil.h


 * Author: tz
 *
 * Created on 17. August 2012, 13:29
 */

#include <algorithm>

#include "General/environment.hh"
#include "Domain/resultInfo.hh"


#ifndef _DIVUTIL_H
#define	_DIVUTIL_H

namespace CCM {

  //typedef unsigned int uint;

  void WaitReturn();

  void Exit();

  template<class T>
  inline T GetMin(T argA, T argB) {
    return argA < argB ? argA : argB;
  }

  template<class T>
  inline T GetMax(T argA, T argB) {
    return argA < argB ? argB : argA;
  }

  template<class T>
  inline void ApplyMin(T& min, T possibleMin) {
    min = min < possibleMin ? min : possibleMin;
  }

  template<class T>
  inline void ApplyMax(T& max, T possibleMax) {
    max = possibleMax < max ? max : possibleMax;
  }

  template<class T>
  inline bool IsEqualApprox(T valueA, T valueB, T tolerance) {
    T diff = valueA - valueB;
    return diff <= tolerance && diff >= -tolerance;
  }
  
  template<class T>
  inline T* CreateHeapArray(uint size, T& value) {
    T* array = new T[size];
    std::fill(array, array + size, value);
    return array;
  }

  template<class T>
  inline T* CreateZeroArray(uint size) {
    T bla = 0;
    return CreateHeapArray(size, bla);
  }
  
}

#endif	/* _DIVUTIL_H */

