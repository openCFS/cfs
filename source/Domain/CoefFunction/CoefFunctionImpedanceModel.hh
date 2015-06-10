#ifndef COEF_FUNCTION_IMPEDANCEMODEL_HH
#define COEF_FUNCTION_IMPEDANCEMODEL_HH

#include "CoefFunctionTimeFreq.hh"

typedef enum {IMP_NONE, IMP_ZYL_MPP, IMP_INTERPOL, IMP_MUFFLER} IMPEDANCE_TYPE;

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
    CoefFunctionImpedanceModel(MathParser* mp);

    //! Destructor
    virtual ~CoefFunctionImpedanceModel();

    /**
     *  Generates necessary data for MPP impedance model (c_0, density etc.)
     *  @material The material data concerning the acoustic properties and MPP properties
     *  @innerR Inner radius of pipe
     *  @outerR Outer radius of pipe = innerR + volume behind MPP
     */
    void GenerateZylindricMpp(BaseMaterial* const material, Double innerR, Double outerR);
    /**
     * Generates necessary data for impedance model with interpolated data of Z, the impedance
     *  @material The material data concerning the acoustic properties and MPP properties
     */
    void GenerateInterpolImpedance(BaseMaterial* const material);

    /**
     *  Generates necessary data for muffler impedance model (c_0, density etc.)
     *  @material The material data concerning the acoustic properties and muffler properties
     */
    void GenerateMuffler(BaseMaterial* const material);


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
      if (impedanceType_ == IMP_ZYL_MPP) {
        this->Recalculate_zylMpp();
      } else if (impedanceType_ == IMP_MUFFLER) {
        this->Recalculate_muffler();
      } else if (impedanceType_ == IMP_INTERPOL) {
        this->Recalculate_interpol();
      }
      coefScalar =  constCoefScalar_;
    }

    // =========================================================================
    // STRING REPRESENTATION
    // =========================================================================
    //! \copydoc CoefFunctionAnalytic::GetStrScalar
    virtual void GetStrScalar(std::string& real, std::string& imag);
  protected:

    //! Recalculate impedance
    void Recalculate_zylMpp();
    void Recalculate_interpol();
    void Recalculate_muffler();

  private:

    IMPEDANCE_TYPE impedanceType_; // type of impedance

    Double c0_; // speed of sound
    Double density_; // density of medium
    Double nu_; // kinematic viscosity
    Double holeDiam_; // hole diameter or slit width of MPP
    Double plateThick_; // plate thickness
    Double sigma_; // plate porosity
    Double flowMachNr_; // flow mach number
    Double endCorrection_; // end correction term (additional mass around hole)
    Double beta_; // parameter for the non-linear flow term
    Double innerR_; // inner radius of duct
    Double outerR_; // outer radius, duct and expansion chamber
    Double currFrequ_; // current, already calculated frequency

    PtrCoefFct impedanceCoef_real_; // Impedance CoefFunction
    PtrCoefFct impedanceCoef_imag_; // Impedance CoefFunction
};


} // end of namespace
#endif

