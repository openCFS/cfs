// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ELECTROMAGNETICMATERIAL_DATA
#define ELECTROMAGNETICMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling electromagnetic material data
  */

  class ElectroMagneticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectroMagneticMaterial();

    //! Destructor
    ~ElectroMagneticMaterial();

    //! set a scalar string material parameter
    void SetScalar( std::string& param, const MaterialType& matType);

    //! set a scalar real material parameter
    void SetScalar( Double& param, const MaterialType& matType, 
		    const DataType& dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex& param, const MaterialType& matType, 
		    const DataType& dataType );

    //! set a real material tensor
    void SetTensor( Matrix<Double>& param, const MaterialType& matType,
		    const DataType& dataType );

    //! set a complex material tensor
    void SetTensor( Matrix<Complex>& param, const MaterialType& matType,
		    const DataType& dataType );

    //! get a scalar real material parameter
    void GetScalar( Double& param, const MaterialType& matType, 
		    const DataType& dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, const MaterialType& matType, 
		    const DataType& dataType ) const;

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, const MaterialType& matType,
		    const DataType& dataType,
		    const SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, const MaterialType& matType,
		    const DataType& dataType,
		    const SubTensorType = FULL ) const;	

    //Initialize approximations of nonlinear curves
    void InitApproxCurves();

    ApproxData* GetNonlinFncBH() {
      return nlinFncBH_;
    }

  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;

    ApproxData* nlinFncBH_;

  };

} // end of namespace

#endif
