// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionMapping.hh
 *       \brief    <Description>
 *
 *       \date     Mar 12, 2016
 *       \author   sschoder
 */
//================================================================================================

#ifndef COEFFUNCTIONMAPPING_HH_
#define COEFFUNCTIONMAPPING_HH_

#define _USE_MATH_DEFINES
#include <math.h>

#include "CoefFunctionPML.hh"

namespace CoupledField{


template<typename T>
class CoefFunctionMapping : public CoefFunctionPML<T> {

public:


  CoefFunctionMapping(PtrParamNode pmlDef, PtrCoefFct speedOfSound,
                  shared_ptr<EntityList> EntList,
                  StdVector<RegionIdType> pdeDomains,
                  bool isVector );

  virtual ~CoefFunctionMapping();

  virtual string GetName() const { return "CoefFunctionMapping"; }

  //! Return real-valued tensor at integration point
  void GetTensor(Matrix<Double>& tensor,
                 const LocPointMapped& lpm );
  void GetTensor(Matrix<Complex>& tensor,
                 const LocPointMapped& lpm );

  //! Return real-valued vector at integration point
  void GetVector(Vector<Double>& vec,
                 const LocPointMapped& lpm ) ;
  void GetVector(Vector<Complex>& vec,
                 const LocPointMapped& lpm ) ;

 //! Return real-valued scalar at integration point
  // this is little bit of a hack,
  // seeing that the jacobian is transformed according to the changed
  // derivatives, we pass this function as a scalar function to the bilinearform
  // an transform the jacobian with it....
 void GetScalar(Double& val,
                const LocPointMapped& lpm ) ;
 void GetScalar(Complex& val,
                const LocPointMapped& lpm ) ;

};

}
#endif /* COEFFUNCTIONMAPPING_HH_ */
