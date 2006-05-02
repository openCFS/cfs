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


    //! compute the correct subTensor (3D, AXI, ..)
    void ComputeSubTensor(Matrix<Complex>& matMatrix,
			  const MaterialType& matType, 
			  const SubTensorType& subTensor) const;

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
