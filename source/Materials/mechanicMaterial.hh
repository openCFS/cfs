// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef MECHANICMATERIAL_DATA
#define MECHANICMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling mechanic material data
  */

  class MechanicMaterial : public BaseMaterial {

  public:

    //! Default constructor
    MechanicMaterial();

    //! Destructor
    ~MechanicMaterial();

    //! Trigger finalization of mataterial (calculation of rotated matrices)
    void Finalize();

   //! set a scalar string material parameter
    void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a real vector
    virtual void SetVector(const Vector<Double>& param, MaterialType matType,
			    Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! set a complex material tensor
    void SetTensor(const Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType );

    //! get a scalar string material parameter
    void GetScalar( std::string& param, 
                    MaterialType matType ) const;

     //! get a scalar integer material parameter
    void GetScalar( Integer& param, 
                    MaterialType matType) const;

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;

    //! get a real vector
    void GetVector( Vector<Double>& param, MaterialType matType,
			    Global::ComplexPart dataType ) const;

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

  private:


    //! compute the2 correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  MaterialType matType, 
			  SubTensorType subTensor) const;

    //! Compute elasticity tensor from given parameters
    void ComputeFullStiffTensor();

    Double density_;
    Double PoissonRatio_;
    Double RayleighAlpha_;
    Double RayleighBeta_;
    Double RayleighFrequency_;
    Double lossTangens_;

    Complex scalarEmodulus_;
    Complex scalarLameLambda_;
    Complex scalarLameMu_;

    Matrix<Complex> stiffnessTensor_;

  };

} // end of namespace

#endif
