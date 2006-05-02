#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "baseMaterial.hh"



namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  BaseMaterial::BaseMaterial() {

    ENTER_FCN("BaseMaterial::BaseMaterial");

    materialDatabaseName_ = "";
    matFileName_ = "";
    nonlinFileName_ = "";
    isScalar    = FALSE;
    isIsotrop   = FALSE;
    isOrthotrop = FALSE;
    isTensor    = FALSE;

  }


   BaseMaterial::~BaseMaterial() {

    ENTER_FCN("BaseMaterial::~BaseMaterial");

  }

  bool BaseMaterial::IsSet( MaterialType matType ) const {
    ENTER_FCN(" BaseMaterial::IsSet" );

    bool found = false;
    scalarMap::const_iterator scalarIt = scalarParams_.find( matType );
    tensorMap::const_iterator tensorIt = tensorParams_.find( matType );

    if( scalarIt != scalarParams_.end() || 
        tensorIt != tensorParams_.end() ) {
      found = true;
    }

    return found;
  }

  void  BaseMaterial::matTypeNotAllowed( const MaterialType& matType, 
					 std::string dim ) const {
    ENTER_FCN("BaseMaterial::matTypeNotAllowed");

    std::string msg, help;
    Enum2String(matType, help);
    msg = "Material type (" + dim + ") " + help + 
          " is not available for " + materialDatabaseName_ + " Database";
    Error(msg.c_str(), __FILE__, __LINE__);
  }
  

  void  BaseMaterial::dataTypeNotAllowed4SetGet( const DataType& dataType, 
						 const std::string msg ) const {
    ENTER_FCN("BaseMaterial::dataTypeNotAllowed4Set");

    std::string msgAll, help;
    Enum2String( dataType, help );
    msgAll = "Datatype " + help + " is not allowed in function " 
           + msg;
    Error(msgAll.c_str(), __FILE__, __LINE__);
  }


  void  BaseMaterial::dataTypeNotAllowed( const DataType& dataType, 
					  const MaterialType& matType ) const {
    ENTER_FCN("BaseMaterial::dataTypeNotAllowed");

    std::string msg, help1, help2;
    Enum2String( dataType, help1 );
    Enum2String( matType, help2 );
    msg = "Datatype " + help1 + " is not allowed for material type " 
           + help2 + " in material data base " + materialDatabaseName_;
    Error(msg.c_str(), __FILE__, __LINE__);
  }

  void BaseMaterial::matTypeNotInDataBase( const MaterialType& matType,
					   std::string dim ) const {
    ENTER_FCN("BaseMaterial::matTypeNotInFile");

    std::string msg, help;
    Enum2String(matType, help);
    msg = "Material type (" + dim + ") " + help + " was not read form/defined in material file";
    Error(msg.c_str(), __FILE__, __LINE__);
  }


  void  BaseMaterial::setMakesNoSense( const DataType& dataType, 
				       std::string msg ) const {
    ENTER_FCN("BaseMaterial::setMakesNoSense");

    std::string msgAll, help;
    Enum2String( dataType, help );
    msgAll = "Set of " + msg + " makes no sense with datatype " + help;
    Error(msgAll.c_str(), __FILE__, __LINE__);
  }


  void  BaseMaterial::subTensorNotAvailable( const MaterialType& matType, 
					     const SubTensorType subTensor ) const {
    ENTER_FCN("BaseMaterial::subTensorNotAvailable") 

    std::string msg, help1, help2;
    Enum2String(matType, help1);
    Enum2String(subTensor, help2);
    msg = "Subtensor " + help2 + " not available for material type " + help1;
    Error(msg.c_str(), __FILE__, __LINE__);
  }

}
