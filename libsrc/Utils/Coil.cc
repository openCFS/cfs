#include "Coil.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/grid.hh"
#include <fstream>



namespace CoupledField {

  // ------------------------------------
  //   Default Constructor (disallowed)
  // ------------------------------------
  Coil::Coil() {};


  // ---------------------------------
  //   Copy Constructor (disallowed)
  // ---------------------------------
  Coil::Coil( const Coil &c ) {};


  // ----------------------
  //   Default Destructor
  // ----------------------
  Coil::~Coil() {
    ENTER_FCN( "Coil::~Coil" );
    if ( fileL_ != NULL ) {
      fileL_->close();
      delete fileL_;
    }
    if ( fileU_ != NULL ) {
      fileU_->close();
      delete fileU_;
    }
  };


  // -----------------------
  //   Allowed Constructor
  // -----------------------
  Coil::Coil( RegionIdType coilRegionId, 
	      std::string pdeName, 
	      Grid * ptGrid ) {
    ENTER_FCN( "Coil::Coil" );

    // We already know our name
    coilRegionId_ = coilRegionId;

    // Set all attributes to default values
    windingCrossSection_ = 0;
    value_               = 0;
    phase_               = 0;
    resistance_          = 0;
    id_                  = 0;
    dynamicsFile_        = "none";
    isRotational_        = false;
    flowDir_             = NODIR;
    saveFileL_           = "none";
    saveFileU_           = "none";
    fileL_               = NULL;
    fileU_               = NULL;

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // **************************
    //   Determine type of coil
    // **************************
    coilRegionName_ = ptGrid->RegionIdToName(coilRegionId_);
    params->GetCoilType( coilTypeS_, coilRegionName_, pdeName );

    if ( coilTypeS_ == "measurementCoil2d" ) {
      coilType_ = MEASUREMENT2D;
    }
    else if ( coilTypeS_ == "measurementCoil3d" ) {
      coilType_ = MEASUREMENT3D;
    }
    else if ( coilTypeS_ == "voltageCoil2d" ) {
      coilType_ = VOLTAGE2D;
    }
    else if ( coilTypeS_ == "voltageCoil3d" ) {
      coilType_ = VOLTAGE3D;
    }
    else if ( coilTypeS_ == "currentCoil2d" ) {
      coilType_ = CURRENT2D;
    }
    else if ( coilTypeS_ == "currentCoil3d" ) {
      coilType_ = CURRENT3D;
    }
    else {
      Info->Error( "Encountered unsupported type of coil: " + coilTypeS_,
		   __FILE__, __LINE__ );
    }

    // *******************************************
    //   Read Parameters for 2D Measurement Coil
    // *******************************************
    if ( coilType_ == MEASUREMENT2D ) {

      attrVec = "", "", "name";
      valVec  = "", "", coilRegionName_;

      keyVec  = pdeName, "coils", "measurementCoil2d", "windingCrossSection";
      params->Get( keyVec, attrVec, valVec, windingCrossSection_ );

      keyVec  = pdeName, "coils", "measurementCoil2d", "saveFileL";
      params->Get( keyVec, attrVec, valVec, saveFileL_ );

      keyVec  = pdeName, "coils", "measurementCoil2d", "saveFileU";
      params->Get( keyVec, attrVec, valVec, saveFileU_ );

      keyVec  = pdeName, "coils", "measurementCoil2d", "id";
      params->Get( keyVec, attrVec, valVec, id_ );
    }

    // *******************************************
    //   Read Parameters for 3D Measurement Coil
    // *******************************************
    if ( coilType_ == MEASUREMENT3D ) {

      attrVec = "", "", "name";
      valVec  = "", "", coilRegionName_;

      keyVec  = pdeName, "coils", "measurementCoil3d", "windingCrossSection";
      params->Get( keyVec, attrVec, valVec, windingCrossSection_ );

      keyVec  = pdeName, "coils", "measurementCoil3d", "saveFileL";
      params->Get( keyVec, attrVec, valVec, saveFileL_ );

      keyVec  = pdeName, "coils", "measurementCoil3d", "saveFileU";
      params->Get( keyVec, attrVec, valVec, saveFileU_ );
    }

    // ***************************************
    //   Read Parameters for 2D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE2D ) {

      attrVec = "", "", "name";
      valVec  = "", "", coilRegionName_;

      keyVec  = pdeName, "coils", "voltageCoil2d", "windingCrossSection";
      params->Get( keyVec, attrVec, valVec, windingCrossSection_ );

      keyVec  = pdeName, "coils", "voltageCoil2d", "value";
      params->Get( keyVec, attrVec, valVec, value_ );

      keyVec  = pdeName, "coils", "voltageCoil2d", "dynamics";
      params->Get( keyVec, attrVec, valVec, dynamicsFile_ );

      keyVec  = pdeName, "coils", "voltageCoil2d", "phase";
      params->Get( keyVec, attrVec, valVec, phase_ );

      keyVec  = pdeName, "coils", "voltageCoil2d", "resistance";
      params->Get( keyVec, attrVec, valVec, resistance_ );

      keyVec  = pdeName, "coils", "voltageCoil2d", "id";
      params->Get( keyVec, attrVec, valVec, id_ );
    }

    // ***************************************
    //   Read Parameters for 3D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE3D ) {

      attrVec = "", "", "name";
      valVec  = "", "", coilRegionName_;

      keyVec  = pdeName, "coils", "voltageCoil3d", "windingCrossSection";
      params->Get( keyVec, attrVec, valVec, windingCrossSection_ );

      keyVec  = pdeName, "coils", "voltageCoil3d", "value";
      params->Get( keyVec, attrVec, valVec, value_ );

      keyVec  = pdeName, "coils", "voltageCoil3d", "dynamics";
      params->Get( keyVec, attrVec, valVec, dynamicsFile_ );

      keyVec  = pdeName, "coils", "voltageCoil3d", "phase";
      params->Get( keyVec, attrVec, valVec, phase_ );

      keyVec  = pdeName, "coils", "voltageCoil3d", "resistance";
      params->Get( keyVec, attrVec, valVec, resistance_ );
    }

    // ***************************************
    //   Read Parameters for 2D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT2D ) {

      attrVec = "", "", "name";
      valVec  = "", "", coilRegionName_;

      keyVec  = pdeName, "coils", "currentCoil2d", "windingCrossSection";
      params->Get( keyVec, attrVec, valVec, windingCrossSection_ );

      keyVec  = pdeName, "coils", "currentCoil2d", "value";
      params->Get( keyVec, attrVec, valVec, value_ );

      keyVec  = pdeName, "coils", "currentCoil2d", "dynamics";
      params->Get( keyVec, attrVec, valVec, dynamicsFile_ );

      keyVec  = pdeName, "coils", "currentCoil2d", "phase";
      params->Get( keyVec, attrVec, valVec, phase_ );

      keyVec  = pdeName, "coils", "currentCoil2d", "id";
      params->Get( keyVec, attrVec, valVec, id_ );
    }

    // ***************************************
    //   Read Parameters for 3D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT3D ) {

      attrVec = "", "", "name";
      valVec  = "", "", coilRegionName_;

      keyVec  = pdeName, "coils", "currentCoil3d", "windingCrossSection";
      params->Get( keyVec, attrVec, valVec, windingCrossSection_ );

      keyVec  = pdeName, "coils", "currentCoil3d", "value";
      params->Get( keyVec, attrVec, valVec, value_ );

      keyVec  = pdeName, "coils", "currentCoil3d", "dynamics";
      params->Get( keyVec, attrVec, valVec, dynamicsFile_ );

      keyVec  = pdeName, "coils", "currentCoil3d", "phase";
      params->Get( keyVec, attrVec, valVec, phase_ );
    }

    // *******************************************
    //   Open results files for measurement coil
    // *******************************************
    if ( coilType_ == MEASUREMENT2D || coilType_ == MEASUREMENT3D ) {

      // Open file stream for storing inductivity
      if ( saveFileL_ != "none" ) {

	std::string msg = " Inductivity is stored in: " + saveFileL_ + '\n';
	Info->PrintF( pdeName, "%s", msg.c_str() );

	fileL_ = new std::ofstream( saveFileL_.c_str() );

	if ( fileL_ == NULL ) {
	  Info->Error( "Could not open " + saveFileL_, __FILE__, __LINE__ );
	}
      }

      // Open file stream for storing current/voltages
      if ( saveFileL_ != "none" ) {

	std::string msg = " Currents/voltages are stored in: "
	  + saveFileU_ + '\n';
	Info->PrintF( pdeName, "%s", msg.c_str() );

	fileL_ = new std::ofstream( saveFileU_.c_str() );

	if ( fileL_ == NULL ) {
	  Info->Error( "Could not open " + saveFileU_, __FILE__, __LINE__ );
	}
      }
    }

    // ***********************
    //   Read flow direction
    // ***********************
    StdVector<std::string> aux;
    if ( coilType_ == CURRENT3D ) {
      
      // Check for currentFlow specification
      keyVec  = pdeName, "coils", "currentCoil3d", "currentFlow";
      attrVec = "", "", "name";
      valVec  = "", "", coilRegionName_;
      params->GetList( keyVec, attrVec, valVec, aux );

      if ( aux.GetSize() == 1 ) {
	if ( aux[0] == "xDir" ) {
	  flowDir_ = XDIR;
	}
	else if ( aux[0] == "yDir" ) {
	  flowDir_ = YDIR;
	}
	else if ( aux[0] == "zDir" ) {
	  flowDir_ = ZDIR;
	}
	else {
	  Info->Error( "Unknown currentFlow " + aux[0], __FILE__, __LINE__ );
	}
      }
      else if ( aux.GetSize() > 1 ) {
	Info->Error( "More than 1 currentFlow specification for coil " +
		     coilRegionName_, __FILE__, __LINE__ );
      }

      // Check for rotational specification
      else {
	isRotational_ = true;

	keyVec  = pdeName, "coils", "currentCoil3d", "midPointX";
	params->Get( keyVec, attrVec, valVec, midX_ );

	keyVec  = pdeName, "coils", "currentCoil3d", "midPointY";
	params->Get( keyVec, attrVec, valVec, midY_ );

	keyVec  = pdeName, "coils", "currentCoil3d", "midPointZ";
	params->Get( keyVec, attrVec, valVec, midZ_ );

	keyVec  = pdeName, "coils", "currentCoil3d", "rotAxisX";
	params->Get( keyVec, attrVec, valVec, rotAxisX_ );

	keyVec  = pdeName, "coils", "currentCoil3d", "rotAxisY";
	params->Get( keyVec, attrVec, valVec, rotAxisY_ );

	keyVec  = pdeName, "coils", "currentCoil3d", "rotAxisZ";
	params->Get( keyVec, attrVec, valVec, rotAxisZ_ );
      }

    }

  }

}
