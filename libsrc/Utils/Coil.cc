#include "Coil.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
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
  Coil::Coil( std::string coilName, std::string pdeName ) {
    ENTER_FCN( "Coil::Coil" );

    // We already know our name
    coilName_ = coilName;

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

    // **************************
    //   Determine type of coil
    // **************************
    params->GetCoilType( coilTypeS_, coilName_, pdeName );

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
      params->CGet( "windingCrossSection_", windingCrossSection_, "name",
		    coilName_, 1, pdeName, "coils");
      params->CGet( "saveFileL", saveFileL_, "name", coilName_, 1, pdeName,
		    "coils" );
      params->CGet( "saveFileU", saveFileU_, "name", coilName_, 1,
		    pdeName, "coils" );
      params->CGet( "id", id_, "name", coilName_, 1, pdeName, "coils" );
    }

    // *******************************************
    //   Read Parameters for 3D Measurement Coil
    // *******************************************
    if ( coilType_ == MEASUREMENT3D ) {
      params->CGet( "windingCrossSection", windingCrossSection_, "name",
		    coilName_, 1, pdeName, "coils");
      params->CGet( "saveFileL", saveFileL_, "name", coilName_, 1, pdeName,
		    "coils" );
      params->CGet( "saveFileU", saveFileU_, "name", coilName_, 1,
		    pdeName, "coils" );
    }

    // ***************************************
    //   Read Parameters for 2D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE2D ) {
      params->CGet( "windingCrossSection", windingCrossSection_ ,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "value"     , value_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "dynamics"  , dynamicsFile_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "phase"     , phase_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "resistance", resistance_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "id"        , id_,
		    "name", coilName_, 1,  pdeName, "coils" );
    }

    // ***************************************
    //   Read Parameters for 3D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE3D ) {
      params->CGet( "windingCrossSection", windingCrossSection_ ,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "value"     , value_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "dynamics"  , dynamicsFile_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "phase"     , phase_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "resistance", resistance_,
		    "name", coilName_, 1, pdeName, "coils" );
    }

    // ***************************************
    //   Read Parameters for 2D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT2D ) {
      params->CGet( "windingCrossSection", windingCrossSection_ ,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "value"   , value_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "dynamics", dynamicsFile_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "phase"   , phase_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "id"      , id_,
		    "name", coilName_, 1, pdeName, "coils" );
    }

    // ***************************************
    //   Read Parameters for 3D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT3D ) {
      params->CGet( "windingCrossSection", windingCrossSection_ ,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "value"   , value_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "dynamics", dynamicsFile_,
		    "name", coilName_, 1, pdeName, "coils" );
      params->CGet( "phase"   , phase_,
		    "name", coilName_, 1, pdeName, "coils" );
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
      params->CGetList( "currentFlow", aux, "name", coilName_, 1, pdeName,
			"coils" );
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
	Info->Error( "More than 1 currentFlow specifications for coil " +
		     coilName_, __FILE__, __LINE__ );
      }

      // Check for rotational specification
      else {
	isRotational_ = true;
	params->CGet( "midPointX", midX_, "name", coilName_, 2, pdeName,
		      "coils" );
	params->CGet( "midPointY", midY_, "name", coilName_, 2, pdeName,
		      "coils" );
	params->CGet( "midPointZ", midZ_, "name", coilName_, 2, pdeName,
		      "coils" );
	params->CGet( "rotAxisX", rotAxisX_, "name", coilName_, 2, pdeName,
		      "coils" );
	params->CGet( "rotAxisY", rotAxisY_, "name", coilName_, 2, pdeName,
		      "coils" );
	params->CGet( "rotAxisZ", rotAxisZ_, "name", coilName_, 2, pdeName,
		      "coils" );
      }

    }

  }

}
