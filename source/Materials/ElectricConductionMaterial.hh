// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ELECTRICCONDUCTIONMATERIAL_DATA
#define ELECTRICCONDUCTIONMATERIAL_DATA

#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling electric current conduction material data
  */

  class ElectricConductionMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectricConductionMaterial(MathParser* mp,
                 CoordSystem * defaultCoosy);

    //! Destructor
    ~ElectricConductionMaterial();

    //! Trigger finalization of material
    void Finalize();

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    Global::ComplexPart dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, 
                   MaterialType matType, 
                   Global::ComplexPart  dataType );

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

  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;
    

    //! Calculate full tensor from scalar values
    void ComputeFullMuTensor();

  };

} // end of namespace

#endif
