#ifndef FLOWMATERIAL_DATA
#define FLOWMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling flow material data
  */

  class FlowMaterial : public BaseMaterial {

  public:

    //! Default constructor
    FlowMaterial();

    //! Destructor
    ~FlowMaterial();

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


    //! get a scalar real material parameter
    void GetScalar( Double& param, const MaterialType& matType, 
		    const DataType& dataType ) const;

    //! get a scalar complex real material parameter
    void GetScalar( Complex& param, const MaterialType& matType, 
		    const DataType& dataType ) const;


  private:

    Double density_;
    Double dynamicViscosity_;

  };

} // end of namespace

#endif
