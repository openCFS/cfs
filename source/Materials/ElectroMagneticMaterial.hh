#ifndef ELECTROMAGNETICMATERIAL_DATA
#define ELECTROMAGNETICMATERIAL_DATA

#include "BaseMaterial.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Utils/LinInterpolate.hh"

namespace CoupledField {

  //!Class for handling electromagnetic material data
  class ElectroMagneticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectroMagneticMaterial(MathParser* mp, CoordSystem * defaultCoosy, bool isDarwin = false);


    //! Destructor
    ~ElectroMagneticMaterial();
    
    //! Trigger finalization of material
    void Finalize();

    // ======================================================================
    //  Coefficient Function Related Methods
    // ======================================================================
    //@{ \name Coefficient Function Related Method
    //! Return tensor-valued coefficient function for nonlinear function
    virtual PtrCoefFct GetTensorCoefFncNonLin( MaterialType matType,
                                               SubTensorType type,
                                               Global::ComplexPart matDataType,
                                               PtrCoefFct dependency );

    //! Return scalar-valued coefficient function for nonlinear function 
    virtual PtrCoefFct GetScalCoefFncNonLin(MaterialType matType,
                                            Global::ComplexPart matDataType,
                                            PtrCoefFct fluxCoef,
                                            PtrCoefFct tempCoef = NULL );
 
    //! only valid for magnetostrictive coupling; nu = nu(S)
    virtual PtrCoefFct GetScalCoefFncNonLin_MagStrict(MaterialType matType,
                                            Global::ComplexPart matDataType,
                                            PtrCoefFct mechStrain );
    //! Return scalar-valued coefficient function for derivative w.r.t. parameter 
    virtual PtrCoefFct GetScalCoefFncNonLinDerivParam(MaterialType matType,
                                                      Global::ComplexPart matDataType,
                                                      PtrCoefFct fluxCoef);                                                                                
    //@}

  private:

    
    //! Calculate full permeability and reluctivity tensors from scalar values
    void ComputeFullMuTensor();
    
    //! Initialize smooth spline approximation of nonlinear BH curve
    void InitSmoothSplineApprox(MatDescriptorNl& matNl, MaterialType matType);
    
    // Helper to create analytic compound coefficient functions
    PtrCoefFct CreateAnalyticCoefFct(std::string exprStr, std::map<std::string, PtrCoefFct>& symbols);

    //! CoefFunction for anisotropic material which is passed to its derivative
    //! used to calculate an approximation of the derivative with respect to the angle
    shared_ptr<CoefFunctionApproxAniso> baseCoefAniso_;

    //! CoefFunction for temperature-dependent BH-curves
    shared_ptr<CoefFunctionApproxIsotropicTemperatureDependent> baseCoefIsoTempDependBH_;

    bool isDarwin_;
  };

} // end of namespace

#endif
