#include "Coil.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include <fstream>



namespace CoupledField {

  // ------------------------------------
  //   Default Constructor (disallowed)
  // ------------------------------------
  Coil::Coil() {
  lastSaveStep_ = -1;
  }


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
              PtrParamNode coilNode,
              Grid * ptGrid ) {

    lastSaveStep_ = -1;
    
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

    // **************************
    //   Determine type of coil
    // **************************
    // get region name of coil
    coilNode->GetValue( "name", coilRegionName_ );


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
      EXCEPTION( "Encountered unsupported type of coil: " << coilTypeS_ );
    }

    // *******************************************
    //   Read Parameters for 2D Measurement Coil
    // *******************************************
    if ( coilType_ == MEASUREMENT2D ) {

      windingCrossSection_ = 
                coilNode->Get( "windingCrossSection")->MathParse<Double>();
      coilNode->GetValue( "saveFileL",  saveFileL_, ParamNode::PASS );
      coilNode->GetValue( "saveFileU", saveFileU_, ParamNode::PASS );
      coilNode->GetValue( "id", id_ );
    }

    // *******************************************
    //   Read Parameters for 3D Measurement Coil
    // *******************************************
    if ( coilType_ == MEASUREMENT3D ) {
      
      windingCrossSection_ = 
          coilNode->Get( "windingCrossSection")->MathParse<Double>();
      coilNode->GetValue( "saveFileL",  saveFileL_, ParamNode::PASS );
      coilNode->GetValue( "saveFileU", saveFileU_, ParamNode::PASS );
    }

    // ***************************************
    //   Read Parameters for 2D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE2D ) {

      windingCrossSection_ = 
                coilNode->Get( "windingCrossSection")->MathParse<Double>();
      coilNode->GetValue( "value", value_ );
      coilNode->GetValue( "phase", phase_, ParamNode::PASS );
      coilNode->GetValue( "resistance", resistance_ );
      coilNode->GetValue( "id", id_ );
    }

    // ***************************************
    //   Read Parameters for 3D Voltage Coil
    // ***************************************
    else if ( coilType_ == VOLTAGE3D ) {

      windingCrossSection_ = 
                coilNode->Get( "windingCrossSection")->MathParse<Double>();
      coilNode->GetValue( "value", value_ );
      coilNode->GetValue( "phase", phase_, ParamNode::PASS );
      coilNode->GetValue( "resistance", resistance_ );

    }

    // ***************************************
    //   Read Parameters for 2D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT2D ) {

      windingCrossSection_ = 
                coilNode->Get( "windingCrossSection")->MathParse<Double>();
      coilNode->GetValue( "value", value_ );
      coilNode->GetValue( "phase", phase_, ParamNode::PASS );
      coilNode->GetValue( "saveFileL",  saveFileL_, ParamNode::PASS );
      coilNode->GetValue( "id", id_ );
      
    }

    // ***************************************
    //   Read Parameters for 3D Current Coil
    // ***************************************
    else if ( coilType_ == CURRENT3D ) {

      windingCrossSection_ = 
                coilNode->Get( "windingCrossSection")->MathParse<Double>();
      coilNode->GetValue( "value", value_ );
      coilNode->GetValue( "phase", phase_, ParamNode::PASS );
      coilNode->GetValue( "saveFileL",  saveFileL_, ParamNode::PASS );

      // Try to lead flow direction
      ParamNodeList dirNodes= 
        coilNode->Get("flowDirection")->GetList("direction");
      
      std::string refCoordSysName = 
        coilNode->Get("flowDirection")->Get("coordSysId")->As<std::string>();
      flowCoordSys_ = domain->GetCoordSystem( refCoordSysName );

      locFlowDir_.Resize( flowCoordSys_->GetDim());
      locFlowDir_.Init();
      for( UInt i = 0; i < dirNodes.GetSize(); i++ ) {
        std::string dir = dirNodes[i]->Get("dof")->As<std::string>();
        Double val = dirNodes[i]->Get("value")->As<Double>();
        // Get local vector component index
        UInt index = flowCoordSys_->GetVecComponent(dir);
        locFlowDir_[index-1] = val;
      }
      
      
    }

    // *******************************************
    //   Open results files for measurement coil
    // *******************************************
    if ( coilType_ == CURRENT2D 
        || coilType_ == MEASUREMENT2D
        || coilType_ == MEASUREMENT3D
        || coilType_ == CURRENT3D ) {

      // Open file stream for storing inductivity
      if ( saveFileL_ != "none" ) {

        // TODO better use the magPDE infoNode_
        info->Get("PDE/coil/inductivity/file")->SetValue(saveFileL_);

        fileL_ = new std::ofstream( saveFileL_.c_str() );

        if ( fileL_ == NULL ) {
          EXCEPTION( "Could not open " << saveFileL_ );
        }
      }
    }

    if ( coilType_ == MEASUREMENT2D || coilType_ == MEASUREMENT3D ) {
      // Open file stream for storing current/voltages
      if ( saveFileU_ != "none" ) {

        // TODO better use the magPDE infoNode_
        info->Get("PDE/coil/currents_voltages/file")->SetValue(saveFileU_);

        fileU_ = new std::ofstream( saveFileU_.c_str() );

        if ( fileU_ == NULL ) {
          EXCEPTION( "Could not open " << saveFileU_ );
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
//           EXCEPTION( "Unknown currentFlow " + aux[0] );
//         }
//       }
//       else if ( aux.GetSize() > 1 ) {
//         EXCEPTION( "More than 1 currentFlow specification for coil " +
//                      coilRegionName_);
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
