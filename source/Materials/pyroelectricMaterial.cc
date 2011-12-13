#include <complex>
#include <map>
#include <set>
#include <string>
#include <utility>

#include "General/exception.hh"
#include "MatVec/matrix.hh"
#include "Materials/baseMaterial.hh"
#include "limits.h"
#include "math.h"
#include "pyroelectricMaterial.hh"
#include "stdlib.h"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  PyroelectricMaterial::PyroelectricMaterial() : BaseMaterial() {

    materialDatabaseName_ = "Pyroelectric";

    //set the allowed material parameters
	isAllowed_.insert( PYROCOEFFICIENT_TENSOR );
  }

  PyroelectricMaterial::~PyroelectricMaterial() {

    materialDatabaseName_ = "Pyroelectric";

  }

  void PyroelectricMaterial::SetTensor(const Matrix<Double>& param, MaterialType matType, Global::ComplexPart dataType ) {
    

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "tensor";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );
      if ( dataType == Global::REAL || dataType == Global::IMAG ) {
	if ( tensorParams_[matType].GetNumRows() == 0 ) {
	  tensorParams_[matType].Resize( param.GetNumRows(), param.GetNumCols() );
          tensorParams_[matType].Init();
	}
	if ( tensorParamsOrig_[matType].GetNumRows() == 0 ) {
	  tensorParamsOrig_[matType].Resize( param.GetNumRows(), param.GetNumCols() );
          tensorParamsOrig_[matType].Init();
	}

	tensorParams_[matType].SetPart( dataType, param );
	tensorParamsOrig_[matType].SetPart( dataType, param );

	// to be consistent to old structure
	if ( dataType == Global::REAL ) {
	  scalarParams_[matType] = Complex( param[2][2], 0.0);
	}
	else {
	  scalarParams_[matType] = Complex( 0.0, param[2][2]);
	  isComplex_.insert( matType );
	}
      }
      else {
	std::string msg = "SetTensor-Double";
	dataTypeNotAllowed4SetGet ( dataType, msg );
      }
    }
  }

  void PyroelectricMaterial::SetTensor( const Matrix<Complex>& param, 
                                        MaterialType matType, 
                                        Global::ComplexPart dataType ) {
    

    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      std::string dim = "scalar";
      matTypeNotAllowed( matType, dim );
    }
    else {
      isSet_.insert( matType );

      if ( dataType != Global::COMPLEX ) {
	std::string msg = "SetTensor with Matrix<Complex>";
	setMakesNoSense( dataType, msg );
      }
      else {
	tensorParams_[matType]     = param;
	tensorParamsOrig_[matType] = param;
	// to be consistent to old structure
	scalarParams_[matType] = param[2][2];
	isComplex_.insert( matType );
      }
    }    
  }

  void PyroelectricMaterial::GetTensor( Matrix<Double>& param, 
                                        MaterialType matType, 
                                        Global::ComplexPart dataType,
                                        SubTensorType subTensor) const {


    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Matrix<Complex> matTensor;
      if ( subTensor == FULL ) {
	matTensor = pos->second;
      }
      else {
	ComputeSubTensor(matTensor, matType, subTensor);
      }

      if ( dataType == Global::REAL || dataType == Global::IMAG) {
	param = matTensor.GetPart( dataType );
      }
      else {
	std::string msg = "GetTensor-Double";
	dataTypeNotAllowed4SetGet( dataType, msg );
      }
    }
  }

  void PyroelectricMaterial::GetTensor( Matrix<Complex>& param, 
                                        MaterialType matType, 
                                        Global::ComplexPart dataType,
                                        SubTensorType subTensor) const {
    

    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    if ( pos == tensorParams_.end() ) {
      std::string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
    else {
      Matrix<Complex> matTensor;
      if ( subTensor == FULL ) {
	matTensor = pos->second;
      }
      else {
	ComputeSubTensor(matTensor, matType, subTensor);
      }

      if ( dataType == Global::REAL || dataType == Global::IMAG) {
	Matrix<Double> help; 
	help = matTensor.GetPart( dataType );
	param.Resize( matTensor.GetNumRows(), matTensor.GetNumCols() );
	param.SetPart( dataType, help );
      }
      else if ( dataType == Global::COMPLEX ) {
	param = matTensor;
      }
    }
  }


  void PyroelectricMaterial::ComputeSubTensor(Matrix<Complex>& matMatrix,
					       const MaterialType& matType, 
					       const SubTensorType& subTensor) const {


    tensorMap::const_iterator pos;
    pos = tensorParams_.find( matType );

    //2D tensor axi or plane is the same
    matMatrix.Resize(2,2);
    matMatrix.Init();
    pos->second.GetSubMatrix(matMatrix, 0, 0);
  }

  void PyroelectricMaterial::SetScalar(Double param, MaterialType matType, 
                                       Global::ComplexPart dataType ) {
    EXCEPTION("SetScalar for 'Double' not implemented");
  }

  void PyroelectricMaterial::GetScalar(Double& param, MaterialType matType, 
                                       Global::ComplexPart dataType ) const{
    EXCEPTION("GetScalar for 'Double' not implemented");
  }

  void PyroelectricMaterial::GetScalar(Complex& param, MaterialType matType, 
                                       Global::ComplexPart dataType ) const{
    EXCEPTION("GetScalar for 'Complex' not implemented");
  }


}
