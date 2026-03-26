#ifndef COEF_FUNCTION_IMPEDANCEMODEL_HH
#define COEF_FUNCTION_IMPEDANCEMODEL_HH

#include "CoefFunctionTimeFreq.hh"
#include "Materials/BaseMaterial.hh"

typedef enum {IMP_NONE, IMP_SLIT_MPP, IMP_CIRC_MPP, IMP_INTERPOL, IMP_FCT} IMPEDANCE_TYPE;

namespace CoupledField{

template<typename T>
class CoefFunctionImpedanceModel : public CoefFunctionTimeFreq<Complex> {
  virtual string GetName() const { return "CoefFunctionImpedanceModel"; }
};

//! Coefficient function for calculating the impedance
template<>
class CoefFunctionImpedanceModel<Complex> : public CoefFunctionTimeFreq<Complex>,
                           public enable_shared_from_this<CoefFunctionImpedanceModel<Complex> > {
  public:
    //! Constructor with given field
    CoefFunctionImpedanceModel(MathParser* mp, \
        BaseMaterial* const material, bool isNormalised=false);

    //! Destructor
    virtual ~CoefFunctionImpedanceModel();

    virtual string GetName() const { return "CoefFunctionImpedanceModel<complex>"; }

    /**
     *  Generates necessary data for MPP impedance model (c_0, density etc.)
     *  @impNode Impedance node with all the necessary information
     */
    void GenerateSlitMpp(const PtrParamNode impNode);
    /**
     *  Generates necessary data for MPP impedance model with circular holes (c_0, density etc.)
     *  @impNode Impedance node with all the necessary information
     */
    void GenerateCircMpp(const PtrParamNode impNode);
    /**
     * Generates necessary data for impedance model with interpolated data of Z, the impedance
     *  @impNode Impedance node with all the necessary information
     */
    void GenerateInterpolImpedance(const PtrParamNode  impNode);

    /**
     *  Generates necessary data for function impedance model (c_0, density etc.)
     *  @impNode Impedance node with all the necessary information
     */
    void GenerateImpedanceFct(const PtrParamNode impNode);



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
      switch (impedanceType_) {
      case IMP_SLIT_MPP:
        this->Recalculate_slitMpp();
        break;
      case IMP_CIRC_MPP:
        this->Recalculate_circMpp();
        break;
      case IMP_INTERPOL:
        this->Recalculate_impFct();
        break;
      case IMP_FCT:
        this->Recalculate_impFct();
        break;
      default:
        EXCEPTION("Impedance type not implemented");
      }
      coefScalar =  constCoefScalar_;
    }

    // =========================================================================
    // STRING REPRESENTATION
    // =========================================================================
    //! \copydoc CoefFunctionAnalytic::GetStrScalar
    virtual void GetStrScalar(std::string& real, std::string& imag);
  protected:

    // Read and initiate necessary values
    void Init(const PtrParamNode impNode);

    //! Recalculate impedance
    void Recalculate_slitMpp();
    void Recalculate_circMpp();
    void Recalculate_impFct();

    //! Calculate cavity impedance
    inline void Calculate_cavityImpedance(Complex& Z_cav, const Double omega);

  private:

    IMPEDANCE_TYPE impedanceType_; // type of impedance

    BaseMaterial * material_; // material data
    bool isNormalised_; // is impedance normalised or not
    Double normalisedFactor_;
    Double c0_; // speed of sound
    Double density_; // density of medium
    Double nu_; // kinematic viscosity
    Double holeDiam_; // hole diameter or slit width of MPP
    Double plateThick_; // plate thickness
    Double sigma_; // plate porosity
    Double mppVolDepth_; // volume depth behind mpp
    Double flowMachNr_; // flow mach number
    // Double endCorrection_; // end correction term (additional mass around hole)
    Double beta_; // parameter for the non-linear flow term
    Double currFrequ_; // current, already calculated frequency

    PtrCoefFct impedanceCoef_real_; // Impedance CoefFunction
    PtrCoefFct impedanceCoef_imag_; // Impedance CoefFunction
};


} // end of namespace
#endif

