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
                                            PtrCoefFct fluxCoef );
 
    //! only valid for magnetostrictive coupling; nu = nu(S)
    virtual PtrCoefFct GetScalCoefFncNonLin_MagStrict(MaterialType matType,
                                            Global::ComplexPart matDataType,
                                            PtrCoefFct mechStrain );
    //@}

//    //============================ HYSTERESIS ===================================
//
//    //Initialize hysteresis
//    //Already defined in base class; never called?
//    //virtual void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
//    //                       bool isInverse = false, bool computeInverse = false );
//
//    //set values for differential material approach
//    virtual void SetPreviousHystVal( UInt nrElem, Vector<Double>& Xval );
//
//    //! compute scalar differential parameter
//    virtual Double ComputeScalarDiffVal( UInt nrElem, Vector<Double>& Xval );
//
//    //! compute scalar differential parameters
//    virtual void ComputeScalarDiffValues( UInt nrElem, Vector<Double>& in,
//                                            Vector<Double>& scalarValues );
//
//
//    //! computes the scalar hystereis value
//    virtual Double ComputeScalarHystVal( UInt nrElem, Vector<Double>& Xval ) {
//      EXCEPTION( "ComputeScalarHystVal not implemented" );
//      return 1.0;
//   };
//
//    //! compute the vector hysteresis values
//    virtual void ComputeVectorHystVal( UInt nrElem, Vector<Double>& Xin,
//                                       Vector<Double>& Yout );
//
//    //!
//    virtual void ComputeInverseScalar( UInt idxEl, UInt comp, Double Yin,
//                               Double& Xout );
//
//    //! computes the scalar hystereis value
//    virtual Double GetScalarHystVal( UInt nrElem );
//
//    virtual void GetVectorHystVal( UInt nrElem, Vector<Double>& Val );
//
//    Double ComputeMatDiff( Vector<Double>& dX, Vector<Double>& dY, UInt idx );
  private:

    //! Calculate full permeability and reluctivity tensors from scalar values
    void ComputeFullMuTensor();
    
    //UInt dim_;

    //Matrix<Double> vecXprevious_; //! previous Xval in hysteresis
    //Matrix<Double> vecYprevious_; //! previous Yval in hysteresis
    //Matrix<Double> vecXact_; //! actual Xval in hysteresis (for inverse hysteresis)
    //Matrix<Double> vecYact_; //! actual Yval in hysteresis (for inverse hysteresis)

    //Vector<Double> matDiffprevious_;

    //! CoefFunction for anisotropic material which is passed to its derivative
    //! used to calculate an approximation of the derivative with respect to the angle
    shared_ptr<CoefFunctionApproxAniso> baseCoefAniso_;

    bool isDarwin_;
  };

} // end of namespace

#endif
