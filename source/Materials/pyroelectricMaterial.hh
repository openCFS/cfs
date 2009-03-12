#ifndef PYROELECTRICMATERIAL_DATA
#define PYROELECTRICMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling thermoelastic material data
  */

  class PyroelectricMaterial : public BaseMaterial {

  public:

    //! Default constructor
    PyroelectricMaterial();

    //! Destructor
    ~PyroelectricMaterial();

    //! set a real material tensor
    void SetTensor(const Matrix<Double>& param, MaterialType matType, Global::ComplexPart dataType );

    //! set a complex material tensor
    void SetTensor( const Matrix<Complex>& param, MaterialType matType,
                    Global::ComplexPart dataType );

    //! get a real material tensor
    void GetTensor( Matrix<Double>& param, const MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;	

    //! get a complex material tensor
    void GetTensor( Matrix<Complex>& param, MaterialType matType,
		    Global::ComplexPart dataType,
		    SubTensorType = FULL ) const;

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );


    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;


    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;
    


  private:

    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;

    //Complex scalarPyrocoefficient_;
  };

} // end of namespace

#endif
