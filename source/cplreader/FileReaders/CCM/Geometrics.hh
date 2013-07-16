/* 
 * File:   geoUtil.h
 * Author: tz
 *
 * Created on 16. August 2012, 14:48
 */

#include <math.h>
#include "divUtil.h"

#ifndef _GEOUTIL_H
#define	_GEOUTIL_H

namespace CCM {

  template<class T>
  T HexaOpposite(T index, T* faces) {
    if (faces[0] == index) {
      return faces[1];
    } else if (faces[1] == index) {
      return faces[0];
    } else if (faces[2] == index) {
      return faces[3];
    } else if (faces[3] == index) {
      return faces[2];
    } else if (faces[4] == index) {
      return faces[5];
    } else if (faces[5] == index) {
      return faces[4];
    }
    return 0;
  }

  template<typename T>
  inline void CreateVector(T* source, T* target, T* vector) {
    vector[0] = target[0] - source[0];
    vector[1] = target[1] - source[1];
    vector[2] = target[2] - source[2];
  }

  template<typename T>
  inline T Length(T* vector) {
    return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
  }

  template<typename T>
  inline T ScalarProduct(T* vectorA, T* vectorB) {
    return vectorA[0] * vectorB[0] + vectorA[1] * vectorB[1] + vectorA[2] * vectorB[2];
  }

  template<typename T>
  inline void Product(T* vector, T value) {
    vector[0] *= value;
    vector[1] *= value;
    vector[2] *= value;
  }

  template<typename T>
  inline void Copy(T* target, T* source) {
    target[0] = source[0];
    target[1] = source[1];
    target[2] = source[2];
  }

  template<typename T>
  inline void Add(T* vectorA, T* vectorB) {
    vectorA[0] += vectorB[0];
    vectorA[1] += vectorB[1];
    vectorA[2] += vectorB[2];
  }

  template<typename T>
  inline void Substract(T* vectorA, T* vectorB) {
    vectorA[0] -= vectorB[0];
    vectorA[1] -= vectorB[1];
    vectorA[2] -= vectorB[2];
  }

  template<typename T>
  inline void CrossProduct(T* vectorA, T* vectorB, T* product) {
    product[0] = vectorA[1]*vectorB[2] - vectorA[2]*vectorB[1];
    product[1] = vectorA[2]*vectorB[0] - vectorA[0]*vectorB[2];
    product[2] = vectorA[0]*vectorB[1] - vectorA[1]*vectorB[0];
  }

  template<typename T>
  inline void Normalize(T* vector) {
    T length = Length(vector);
    vector[0] = vector[0] / length;
    vector[1] = vector[1] / length;
    vector[2] = vector[2] / length;
  }

  template<typename T>
  inline T Angle(T* vectorA, T* vectorB) {
    return acos(ScalarProduct(vectorA, vectorB) / Length(vectorA) / Length(vectorB));
  }

  template<typename T>
  inline void Inverse(T* vector) {
    vector[0] = -vector[0];
    vector[1] = -vector[1];
    vector[2] = -vector[2];
  }

  template<typename T>
  inline void Center(T* coordinateA, T* coordinateB, T* center) {
    center[0] = (coordinateA[0] + coordinateB[0]) / 2;
    center[1] = (coordinateA[1] + coordinateB[1]) / 2;
    center[2] = (coordinateA[2] + coordinateB[2]) / 2;
  }

  template<typename T>
  inline void Center(T* coordinateA, T* coordinateB, T* coordinateC, T* center) {
    center[0] = (coordinateA[0] + coordinateB[0] + coordinateC[0]) / 3.0;
    center[1] = (coordinateA[1] + coordinateB[1] + coordinateC[1]) / 3.0;
    center[2] = (coordinateA[2] + coordinateB[2] + coordinateC[2]) / 3.0;
  }

  template<typename T>
  inline void Center(T* coords, uint coordCount, T* center) {
    center[0] = 0;
    center[1] = 0;
    center[2] = 0;
    for (uint i=0; i < coordCount; i++) {
      center[0] += coords[i*3];
      center[1] += coords[i*3+1];
      center[2] += coords[i*3+2];
    }
    center[0] /= coordCount;
    center[1] /= coordCount;
    center[2] /= coordCount;
  }

  template<typename T>
  inline T Distance(T* coordsA, T* coordsB) {
    T vect[3];
    CreateVector(coordsA, coordsB, vect);
    return Length(vect);
  }

  template<typename T>
  inline bool IsNull(T* vect, T& distanceTolerance) {
    return IsEqualApprox(vect[0], 0, distanceTolerance) &&
            IsEqualApprox(vect[1], 0, distanceTolerance) &&
            IsEqualApprox(vect[2], 0, distanceTolerance);
  }

  template<typename T>
  inline bool IsNull(T* vect) {
    return vect[0] == 0 && vect[1] == 0 && vect[2] == 0;
  }

}

#endif	/* _GEOUTIL_H */

