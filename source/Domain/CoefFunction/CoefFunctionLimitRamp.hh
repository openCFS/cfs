/*
 * CoefFunctionLimitRamp.hh
 *
 *  Created on: 30 April 2024
 *      Author: Dominik Mayrhofer
 */

#ifndef SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONLIMITRAMP_HH_
#define SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONLIMITRAMP_HH_

#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"


namespace CoupledField {

// ============================================================================
//  Coef Function Limit Ramp
// ============================================================================
//! Coefficient function for calculating a ramp based on a (primary) feFunction.
//! The ramp is defined by the values primaryOffset (offset of the primary 
//! function from zero), primaryPeakVal (peak value for the primary which is
//! used to flatten the ramp again) and coefFuncPeakVal (peak value of the 
//! coefFunction for the value given by primaryPeakVal), where the function 
//! interpolates linearly between zero and coefFuncPeakVal value. 
//! If primaryPeakVal is larger than primaryOffset, all values of the primary 
//! below the offset lead to a zero interpolation value, whereas values above 
//! the defined peak value lead to the value defined by coefFuncPeakVal. If the 
//! peak value primaryPeakVal is lower than primaryOffset, everything is
//! inversed.
//! One possibility to extend the ramp function would be to consider not only 
//! the current value of a surface element but an averaged value. This can be
//! done by passing the additional bool useMeanPres and setting it to true.
//! 
//! \note This class only works for real-valued scalar data on surfaces, but could be extended.

  class CoefFunctionLimitRamp : public CoefFunction {

  public:

    //! Constructor
    CoefFunctionLimitRamp(shared_ptr<BaseFeFunction> feFct, 
                          StdVector<std::string> surfList,
                          StdVector<Double> primaryOffset, 
                          StdVector<Double> primaryPeakVal, 
                          StdVector<Double> coefFuncPeakVal, 
                          StdVector<bool> useMeanPres);

    //! Destructor
    virtual ~CoefFunctionLimitRamp();

    virtual string GetName() const { return "CoefFunctionLimitRamp"; }

    //! Return real-valued vector at integration point
    void GetScalar(Double& scal, const LocPointMapped& surflpm );

    //! Calculate the mean pressure
    void CalcMeanPressure();


  protected:
    //! feFunction
    shared_ptr<BaseFeFunction> feFct_;

    //! Ramp parameters
    StdVector<bool> useMeanPres_;
    StdVector<std::string> surfList_;
    std::set<RegionIdType> surfRegions_;
    StdVector<Double> primaryOffset_;
    StdVector<Double> primaryPeakVal_;
    StdVector<Double> coefFuncPeakVal_;
  };

}


#endif /* SOURCE_DOMAIN_COEFFUNCTION_COEFFUNCTIONLIMITRAMP_HH_ */
