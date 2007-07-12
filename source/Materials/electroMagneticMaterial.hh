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
    void SetScalar(const std::string& param, MaterialType matType);

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    DataType dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    DataType dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType,
		    DataType dataType );

    //! set a complex material tensor
    void SetTensor(const Matrix<Complex>& param, MaterialType matType,
		    DataType dataType );

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    DataType dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    DataType dataType ) const;

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
		    DataType dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    DataType dataType,
		    SubTensorType = FULL ) const;	

    //Initialize approximations of nonlinear curves
    void InitApproxCurves();

    ApproxData* GetNonlinFncBH() {
      return nlinFncBH_;
    }

  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  MaterialType matType, 
			  SubTensorType subTensor) const;

    ApproxData* nlinFncBH_;

  };

} // end of namespace

#endif
