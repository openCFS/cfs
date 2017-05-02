/*
 * CoefFunctionVisco.hh
 *
 *  Created on: 22 Feb 2016
 *      Author: INTERN\ftoth
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONVISCO_HH_
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONVISCO_HH_

#include "CoefFunction.hh"

namespace CoupledField {
  class CoefFunctionVisco : public CoefFunction {

  public:

    //! Constructor
    CoefFunctionVisco(StdVector<PtrCoefFct> relaxationTensors, StdVector<PtrCoefFct> pronyFunctions);

    //! Destructor
    virtual ~CoefFunctionVisco();

    //! Return real-valued vector at integration point
    void GetVector(Vector<Double>& vec, const LocPointMapped& lpm );

    //! return vector size
    UInt GetVecSize() const;

  protected:
    UInt vecSize_;
    StdVector<PtrCoefFct> Cs_; // relaxation tensors
    StdVector<PtrCoefFct> ps_; // Prony functions
    UInt N_; // length of the prony series
  };

}


#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONVISCO_HH_ */
