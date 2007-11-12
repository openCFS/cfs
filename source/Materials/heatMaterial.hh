// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef HEATMATERIAL_DATA
#define HEATMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling heat material data
  */

  class HeatMaterial : public BaseMaterial {

  public:

    //! Default constructor
    HeatMaterial();

    //! Destructor
    ~HeatMaterial();

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    DataType dataType );

    //! set a scalar complex material parameter
    void SetScalar( Complex param, MaterialType matType, 
		    DataType dataType );

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, 
                   MaterialType matType, 
                   DataType  dataType );


    //! set a complex material tensor
    // void SetTensor(const Matrix<Complex>& param, MaterialType matType,  DataType dataType )
    // {
    // matTypeNotAllowed( matType, "tensor");
    //}

    //    //! set a complex material tensor
    //    void SetTensor( Matrix<Complex>& param, const MaterialType& matType,
    //		    const DataType& dataType );

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
    
    

  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;
    

  };

} // end of namespace

#endif
