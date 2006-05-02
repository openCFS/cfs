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

    //! set a scalar integer material parameter
    void SetScalar( Integer& param, const MaterialType& matType);

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

    //! get a string material parameter
    void GetScalar( Integer& param, const MaterialType& matType) const;

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

    //! Print material definition to given output stream
    void Print(std::ostream & out) const;

  private:

    Complex scalarPermeability_;
    Complex scalarReluctivity_;
    Complex scalarConductivity_;
    Complex scalarPermittivity_;

    Matrix<Complex> permeabilityTensor_;
    Matrix<Complex> permittivityTensor_;

    // hysteresis
    Double Esat_;
    Double Psat_;
    Double aJiles_;
    Double alphaJiles_;
    Double kJiles_;
    Double cJiles_;
    Integer directionP_;
    std::string hystType_;
  };

} // end of namespace

#endif
