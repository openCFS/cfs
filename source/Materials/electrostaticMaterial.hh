// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ELECTROSTATICMATERIAL_DATA
#define ELECTROSTATICMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {


  //! Class for Material Data
  /*! 
     Class for handling electrostatic material data
  */

  class Hysteresis;

  class ElectroStaticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectroStaticMaterial();

    //! Destructor
    ~ElectroStaticMaterial();

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

    //! get a string material parameter
    void GetScalar( std::string& param, 
		    MaterialType matType) const;

    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    DataType dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    DataType dataType ) const;


    //! get a scalar integer material parameter
    void GetScalar( Integer& param, 
                    MaterialType matType ) const;
    
    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, MaterialType matType,
		    DataType dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    DataType dataType,
		    SubTensorType = FULL ) const;	


    //======================== for hysteresis =======================================
    //! compute scalar differential parameter
    Double ComputeScalarDiffVal( UInt nrElem, Double Xval );


    //set values for differential material approach
    void SetPreviousHystVal( UInt nrElem, Double Xval );

    Double ComputeScalarHystVal( UInt nrElem, Double Xval );

  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  MaterialType matType, 
			  SubTensorType subTensor) const;

    Complex scalarPermittivity_;

  };

} // end of namespace

#endif
