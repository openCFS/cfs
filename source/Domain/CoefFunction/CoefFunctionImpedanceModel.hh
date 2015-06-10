#ifndef COEF_FUNCTION_IMPEDANCEMODEL_HH
#define COEF_FUNCTION_IMPEDANCEMODEL_HH

#include "CoefFunctionTimeFreq.hh"


namespace CoupledField{

template<typename T>
class CoefFunctionImpedanceModel : public CoefFunctionTimeFreq<Complex> {

};

//! Coefficient function for calculating the impedance
template<>
class CoefFunctionImpedanceModel<Complex> : public CoefFunctionTimeFreq<Complex>,
                           public boost::enable_shared_from_this<CoefFunctionImpedanceModel<Complex> > {
  public:
    //! Constructor with given field
    CoefFunctionImpedanceModel(MathParser* mp, Double c0, Double density, Double innerR, Double outerR);

    //! Destructor
    virtual ~CoefFunctionImpedanceModel();

    //! Return real-valued scalar at integration point
    virtual void GetScalar(Double& scal,
        const LocPointMapped& lpm ) {
      EXCEPTION( "CoefFunction::GetScalar<Double> called: This may not happen. "
          << "Most likely this method is called with a complex-valued "
          << "CoefFunction object." );
    }
    //! \copydoc CoefFunction::GetScalar
    void GetScalar(Complex& coefScalar, const LocPointMapped& lpm ) {
      assert(dimType_ == SCALAR);
      this->Recalculate();
      coefScalar =  constCoefScalar_;
    }

    // =========================================================================
    // STRING REPRESENTATION
    // =========================================================================
    //! \copydoc CoefFunctionAnalytic::GetStrScalar
    virtual void GetStrScalar(std::string& real, std::string& imag);
  protected:

    //! Recalculate impedance
    void Recalculate();

  private:

    Double c0_; // speed of sound
    Double density_; // density of medium
    Double innerR_; // inner radius of duct
    Double outerR_; // outer radius, duct and expansion chamber

  
};


} // end of namespace
#endif

