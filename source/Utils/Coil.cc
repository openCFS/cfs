// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Coil.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/grid.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include <fstream>



namespace CoupledField {

  // ------------------------------------
  //   Default Constructor (disallowed)
  // ------------------------------------
  Coil::Coil() {}


  // ---------------------------------
  //   Copy Constructor (disallowed)
  // ---------------------------------
  Coil::Coil( const Coil &c ) {}


  // ----------------------
  //   Default Destructor
  // ----------------------
  Coil::~Coil() {
    if ( fileL_ != NULL ) {
      fileL_->close();
      delete fileL_;
    }
    if ( fileU_ != NULL ) {
      fileU_->close();
      delete fileU_;
    }
  }


  // -----------------------
  //   Allowed Constructor
  // -----------------------
  Coil::Coil( RegionIdType coilRegionId, 
              ParamNode * coilNode,
              Grid * ptGrid ) {

    // We already know our name
    coilRegionId_ = coilRegionId;

    // Set all attributes to default values
    windingCrossSection_ = 0;
    value_               = "0";
    phase_               = "0";
    resistance_          = 0;
    id_                  = 0;
    isRotational_        = false;
    saveFileL_           = "none";
    saveFileU_           = "none";
    fileL_               = NULL;
    fileU_               = NULL;
    flowCoordSys_        = NULL;

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // **************************
    //   Determine type of coil
    // **************************
    // get region name of coil
    coilNode->Get( "name", coilRegionName_ );


    // get coil-type
    coilTypeS_ = coilNode->GetName();

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

      coilNode->Get( "windingCrossSection", windingCrossSection_ );
      coilNode->Get( "saveFileL",  saveFileL_, false );
      coilNode->Get( "saveFileU", saveFileU_, false );
      coilNode->Get( "id", id_ );
    }

    // *******************************************
    //   Read Parameters for 3D Measurement Coil
    // *******************************************
    if ( coilType_ == MEASUREMENT3D ) {
      
      coilNode->Get( "windingCrossSection", windingCrossSection_ );
      coilNode->Get( "saveFileL",  saveFileL_, false );
      coilNode->Get( "saveFileU", saveFileU_, false );
    }

    // ***************************************
    //   Read Parameters for 2D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE2D ) {

      coilNode->Get( "windingCrossSection", windingCrossSection_ );
      coilNode->Get( "value", value_ );
      coilNode->Get( "phase", phase_, false );
      coilNode->Get( "resistance", resistance_ );
      coilNode->Get( "id", id_ );
    }

    // ***************************************
    //   Read Parameters for 3D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE3D ) {

      coilNode->Get( "windingCrossSection", windingCrossSection_ );
      coilNode->Get( "value", value_ );
      coilNode->Get( "phase", phase_, false );
      coilNode->Get( "resistance", resistance_ );

    }

    // ***************************************
    //   Read Parameters for 2D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT2D ) {

      coilNode->Get( "windingCrossSection", windingCrossSection_ );
      coilNode->Get( "value", value_ );
      coilNode->Get( "phase", phase_, false );
      coilNode->Get( "saveFileL",  saveFileL_, false );
      coilNode->Get( "id", id_ );
      
    }

    // ***************************************
    //   Read Parameters for 3D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT3D ) {

      coilNode->Get( "windingCrossSection", windingCrossSection_ );
      coilNode->Get( "value", value_ );
      coilNode->Get( "phase", phase_, false );
      coilNode->Get( "saveFileL",  saveFileL_, false );

      // Try to lead flow direction
      StdVector<ParamNode *> dirNodes= 
        coilNode->Get("flowDirection")->GetList("direction");
      
      std::string refCoordSysName = 
        coilNode->Get("flowDirection")->Get("coordSysId")->AsString();
      flowCoordSys_ = domain->GetCoordSystem( refCoordSysName );

      locFlowDir_.Resize( flowCoordSys_->GetDim());
      locFlowDir_.Init();
      for( UInt i = 0; i < dirNodes.GetSize(); i++ ) {
        std::string dir = dirNodes[i]->Get("dof")->AsString();
        Double val = dirNodes[i]->Get("value")->AsDouble();
        // Get local vector component index
        UInt index = flowCoordSys_->GetVecComponent(dir);
        locFlowDir_[index-1] = val;
      }
      
      
    }

    // *******************************************
    //   Open results files for measurement coil
    // *******************************************
    if ( coilType_ == CURRENT2D || coilType_ == CURRENT3D ) {

      // Open file stream for storing inductivity
      if ( saveFileL_ != "none" ) {

        std::string msg = " Inductivity is stored in: " + saveFileL_ + '\n';
        Info->PrintF( "magneticPDE", "%s", msg.c_str() );

        fileL_ = new std::ofstream( saveFileL_.c_str() );

        if ( fileL_ == NULL ) {
          Info->Error( "Could not open " + saveFileL_, __FILE__, __LINE__ );
        }
      }
    }

    if ( coilType_ == MEASUREMENT2D || coilType_ == MEASUREMENT3D ) {
      // Open file stream for storing current/voltages
      if ( saveFileU_ != "none" ) {

        std::string msg = " Currents/voltages are stored in: "
          + saveFileU_ + '\n';
        Info->PrintF( "magneticPDE", "%s", msg.c_str() );

        fileU_ = new std::ofstream( saveFileU_.c_str() );

        if ( fileU_ == NULL ) {
          Info->Error( "Could not open " + saveFileU_, __FILE__, __LINE__ );
        }
      }
    }

    // ***********************
    //   Read flow direction
    // ***********************
    StdVector<std::string> aux;

    // Note: Not working, as a re-design of this class is issued anyway

    //    if ( coilType_ == CURRENT3D )
//  {
      
//       // Check for currentFlow specification
//       keyVec  = pdeName, "coils", "currentCoil3d", "currentFlow";
//       attrVec = "", "", "name";
//       valVec  = "", "", coilRegionName_;
//       params->GetList( keyVec, attrVec, valVec, aux );

//       if ( aux.GetSize() == 1 ) {
//         if ( aux[0] == "xDir" ) {
//           flowDir_ = XDIR;
//         }
//         else if ( aux[0] == "yDir" ) {
//           flowDir_ = YDIR;
//         }
//         else if ( aux[0] == "zDir" ) {
//           flowDir_ = ZDIR;
//         }
//         else {
//           Info->Error( "Unknown currentFlow " + aux[0], __FILE__, __LINE__ );
//         }
//       }
//       else if ( aux.GetSize() > 1 ) {
//         Info->Error( "More than 1 currentFlow specification for coil " +
//                      coilRegionName_, __FILE__, __LINE__ );
//       }

//       // Check for rotational specification
//       else {
//         isRotational_ = true;

//         keyVec  = pdeName, "coils", "currentCoil3d", "midPointX";
//         params->Get( keyVec, attrVec, valVec, midX_ );

//         keyVec  = pdeName, "coils", "currentCoil3d", "midPointY";
//         params->Get( keyVec, attrVec, valVec, midY_ );

//         keyVec  = pdeName, "coils", "currentCoil3d", "midPointZ";
//         params->Get( keyVec, attrVec, valVec, midZ_ );

//         keyVec  = pdeName, "coils", "currentCoil3d", "rotAxisX";
//         params->Get( keyVec, attrVec, valVec, rotAxisX_ );

//         keyVec  = pdeName, "coils", "currentCoil3d", "rotAxisY";
//         params->Get( keyVec, attrVec, valVec, rotAxisY_ );

//         keyVec  = pdeName, "coils", "currentCoil3d", "rotAxisZ";
//         params->Get( keyVec, attrVec, valVec, rotAxisZ_ );
//       }

//     }

  }

}
