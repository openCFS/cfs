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

  private:

  };

} // end of namespace

#endif
