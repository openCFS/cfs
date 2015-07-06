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
    ElectroMagneticMaterial(MathParser* mp,
                            CoordSystem * defaultCoosy);

    //! Destructor
    ~ElectroMagneticMaterial();
    
    //! Trigger finalization of material
    void Finalize();

    //! set a scalar string material parameter
    void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! set a complex material tensor
    void SetTensor(const Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	
    
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

    //@}
    //============================ HYSTERESIS ===================================

    //Initialize hysteresis
    virtual void InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                           bool isInverse = false, bool computeInverse = false );

    //set values for differential material approach
    virtual void SetPreviousHystVal( UInt nrElem, Vector<Double>& Xval );

    //! compute scalar differential parameter
    virtual Double ComputeScalarDiffVal( UInt nrElem, Vector<Double>& Xval );

    //! compute scalar differential parameters
    virtual void ComputeScalarDiffValues( UInt nrElem, Vector<Double>& in,
                                            Vector<Double>& scalarValues );


    //! computes the scalar hystereis value
    virtual Double ComputeScalarHystVal( UInt nrElem, Vector<Double>& Xval ) {
      EXCEPTION( "ComputeScalarHystVal not implemented" );
      return 1.0; 
   };

    //! compute the vector hysteresis values
    virtual void ComputeVectorHystVal( UInt nrElem, Vector<Double>& Xin, 
                                       Vector<Double>& Yout );

    //!
    virtual void ComputeInverseScalar( UInt idxEl, UInt comp, Double Yin, 
                               Double& Xout );

    //! computes the scalar hystereis value
    virtual Double GetScalarHystVal( UInt nrElem );

    virtual void GetVectorHystVal( UInt nrElem, Vector<Double>& Val );

    Double ComputeMatDiff( Vector<Double>& dX, Vector<Double>& dY, UInt idx );


  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  MaterialType matType, 
			  SubTensorType subTensor) const;

    
    //! Calculate full tensor from scalar values
    void ComputeFullMuTensor();
    
    UInt dim_;

    Matrix<Double> vecXprevious_; //! previous Xval in hysteresis
    Matrix<Double> vecYprevious_; //! previous Yval in hysteresis
    Matrix<Double> vecXact_; //! actual Xval in hysteresis (for inverse hysteresis)
    Matrix<Double> vecYact_; //! actual Yval in hysteresis (for inverse hysteresis)

    Vector<Double> matDiffprevious_;

    Double Xsat_, Ysat_;

    //! CoefFunction for anisotropic material which is passed to its derivative
    //! used to calculate an approximation of the derivative with respect to the angle
    shared_ptr<CoefFunctionApproxAniso> baseCoefAniso_;

  };

} // end of namespace

#endif
