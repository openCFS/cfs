
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Coil.hh"

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
    if ( fileUI_ != NULL ) {
      fileUI_->close();
      delete fileUI_;
    }
  };


  // -----------------------
  //   Allowed Constructor
  // -----------------------
  Coil::Coil( std::string coilName, std::string pdeName ) {
    ENTER_FCN( "Coil::Coil" );

    // We already know our name
    coilName_ = coilName;

    // Set all attributes to default values
    area_         = 0;
    value_        = 0;
    phase_        = 0;
    resistance_   = 0;
    id_           = 0;
    dynamicsFile_ = "undefined";
    isRotational_ = false;
    saveFileL_    = "none";
    saveFileUI_   = "none";
    fileL_        = NULL;
    fileUI_       = NULL;

    // **************************
    //   Determine type of coil
    // **************************
    std::string myCoilType;
    params->Get( "coils", myCoilType, pdeName );

    if ( myCoilType == "measurementCoil2d" ) {
      myType_ = MEASUREMENT2D;
    }
    else if ( myCoilType == "measurementCoil3d" ) {
      myType_ = MEASUREMENT3D;
    }
    else if ( myCoilType == "voltageCoil2d" ) {
      myType_ = VOLTAGE2D;
    }
    else if ( myCoilType == "voltageCoil3d" ) {
      myType_ = VOLTAGE3D;
    }
    else if ( myCoilType == "currentCoil2d" ) {
      myType_ = VOLTAGE2D;
    }
    else if ( myCoilType == "currentCoil3d" ) {
      myType_ = VOLTAGE3D;
    }
    else {
      Info->Error( "Encountered unsupported type of coil: " + myCoilType,
		   __FILE__, __LINE__ );
    }

    // *******************************************
    //   Read Parameters for 2D Measurement Coil
    // *******************************************
    if ( myType_ == MEASUREMENT2D ) {
      params->CGet( "area", area_, "name", coilName_, pdeName, "coils" );
      params->CGet( "saveFileL", saveFileL_, "name", coilName_, pdeName,
		    "coils" );
      params->CGet( "saveFileUI", saveFileUI_, "name", coilName_, pdeName,
		    "coils" );
      params->CGet( "id"         , id_,
		    "name", coilName_, pdeName, "coils" );
    }

    // *******************************************
    //   Read Parameters for 3D Measurement Coil
    // *******************************************
    if ( myType_ == MEASUREMENT3D ) {
      params->CGet( "area", area_, "name", coilName_, pdeName, "coils" );
      params->CGet( "saveFileL", saveFileL_, "name", coilName_, pdeName,
		    "coils" );
      params->CGet( "saveFileUI", saveFileUI_, "name", coilName_, pdeName,
		    "coils" );
    }

    // ***************************************
    //   Read Parameters for 2D Voltage Coil
    // ***************************************
    else if ( myType_ == VOLTAGE2D ) {
      params->CGet( "area"       , area_ ,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "value"      , value_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "dynamics"   , dynamicsFile_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "phase"      , phase_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( " resistance", phase_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "id"         , id_,
		    "name", coilName_, pdeName, "coils" );
    }

    // ***************************************
    //   Read Parameters for 3D Voltage Coil
    // ***************************************
    else if ( myType_ == VOLTAGE3D ) {
      params->CGet( "area"       , area_ ,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "value"      , value_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "dynamics"   , dynamicsFile_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "phase"      , phase_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( " resistance", phase_,
		    "name", coilName_, pdeName, "coils" );
    }

    // ***************************************
    //   Read Parameters for 2D Current Coil
    // ***************************************
    else if ( myType_ == CURRENT2D ) {
      params->CGet( "area"    , area_ ,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "value"   , value_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "dynamics", dynamicsFile_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "phase"   , phase_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "id"      , id_,
		    "name", coilName_, pdeName, "coils" );
    }

    // ***************************************
    //   Read Parameters for 3D Current Coil
    // ***************************************
    else if ( myType_ == CURRENT3D ) {
      params->CGet( "area"    , area_ ,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "value"   , value_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "dynamics", dynamicsFile_,
		    "name", coilName_, pdeName, "coils" );
      params->CGet( "phase"   , phase_,
		    "name", coilName_, pdeName, "coils" );
    }

    // *******************************************
    //   Open results files for measurement coil
    // *******************************************
    if ( myType_ == Coil::MEASUREMENT ) {

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
	  + saveFileUI_ + '\n';
	Info->PrintF( pdeName, "%s", msg.c_str() );

	fileL_ = new std::ofstream( saveFileUI_.c_str() );

	if ( fileL_ == NULL ) {
	  Info->Error( "Could not open " + saveFileUI_, __FILE__, __LINE__ );
	}
      }

    }

  }

}
