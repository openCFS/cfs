// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//=================================
/*
 * \file   Coil.cc
 * \brief  see Coil.hh
 *
 * \date   unknown
 * \author ahauck, dperchto
 */
//=================================

#include "Coil.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/DefaultCoordSystem.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

#include "Driver/BaseDriver.hh"

#include <fstream>

// header for logging
#include "DataInOut/Logging/LogConfigurator.hh"

// declare logging stream
DECLARE_LOG(coil)
DEFINE_LOG(coil, "coil")

namespace CoupledField {

  Coil::Coil() {EXCEPTION("Must not be called");}


  Coil::Coil( const Coil &c ) {}

  Coil::~Coil() {
  }


  // -----------------------
  //   Allowed Constructor
  // -----------------------
  Coil::Coil( PtrParamNode coilNode, 
              PtrParamNode infoNode,
              Grid * ptGrid,
              MathParser * mp,
              Global::ComplexPart type ) {
    myParam_ = coilNode;
    myInfo_ = infoNode;
    ptGrid_ = ptGrid;
    mParser_ = mp;
    
    // initialize data members
    sourceType_ = NO_SOURCE_TYPE;
    coilId_ = "";

    // obtain coilId
    coilId_ = myParam_->Get("id")->As<std::string>();
    
    isMultHarm_ = false;

    // Read source type (only if present: for measurement coils,
    // no excitation is given)
    if( myParam_->Has("source") ) {
      std::string exType = myParam_->Get("source")->Get("type")->As<std::string>();
      if ( exType == "current" ) {
        sourceType_ = CURRENT;
      } else if ( exType == "voltage" ) {
        sourceType_ = VOLTAGE;
      } else if ( exType == "external" ) {
        sourceType_ = EXTERNAL;
      }

      // Read value and phase and generate scalar coefficient function
      std::string value, phase;
      value = myParam_->Get("source")->Get("value")->As<std::string>();
      phase = myParam_->Get("source")->Get("phase")->As<std::string>();
      srcVal_ = CoefFunction::Generate(mParser_, type, AmplPhaseToReal(value, phase), AmplPhaseToImag(value, phase));
    }else if( myParam_->Has("sourceMultiharmonic")){
      isMultHarm_ = true;
      sourceType_ = CURRENT_MULTHARM;
      std::string value, phase;
      ParamNodeList harmonicList = myParam_->Get("sourceMultiharmonic")->GetList("harmonic");
      //srcValMH_.Resize(harmonicList.GetSize());
      for( UInt h = 0; h < harmonicList.GetSize(); ++h ) {
        PtrParamNode harmNode = harmonicList[h];
        UInt harmVal = harmNode->Get("harmonic")->As<Integer>();
        if( (harmVal==0) && !domain->GetDriver()->IsFullSystem() ){
          EXCEPTION("Zero harmonic excitation is only allowed if the full system is considered!"
                    "\n Therefore switch to <fullSystem>true</fullSystem> in the analysis tag");
        }
        if( (harmVal%2==0) && !domain->GetDriver()->IsFullSystem() ){
          EXCEPTION("Only odd harmonics are allowed for the excitation current in the optimized version!"
                    "\n Therefore switch to <fullSystem>true</fullSystem> in the analysis tag");
        }
        value = harmNode->Get("value")->As<std::string>();
        phase = harmNode->Get("phase")->As<std::string>();
        srcValMH_[harmVal] = CoefFunction::Generate(mParser_, type, AmplPhaseToReal(value, phase), AmplPhaseToImag(value, phase));
      }
    }

    /* code from NACS
    // Print information
    std::string sourceString = "- (measurement)";
    if( sourceType_ == CURRENT)
      sourceString = "current";
    if( sourceType_ == VOLTAGE)
      sourceString = "voltage";
    Info->PrintF( "",
                  "\n\n------------------\n"
                  " COIL DESCRIPTION\n"
                  "------------------\n");
    Info->PrintF("", "Coil-ID:\t\t%s\n" , coilId_.c_str());
    Info->PrintF("", "Source-Type:\t\t%s\n", sourceString.c_str());
    Info->PrintF("", "Source-Value:\t\t%s\n", value_.c_str());
    Info->PrintF("", "Source-Phase:\t\t%s\n\n", phase_.c_str());*/

    // ============================
    //  Loop over parts
    // ============================
    ParamNodeList parts = myParam_->GetList("part");
    for( UInt iPart = 0; iPart < parts.GetSize(); ++iPart ) {
     
      //Info->PrintF( "", "=== PART %d ===\n",iPart+1);

      PtrParamNode actPartNode = parts[iPart];

      shared_ptr<Part> partPtr(new Part());
      Part& actPart = *partPtr;

      // read region lists
      ParamNodeList regions = actPartNode->Get("regionList")->GetChildren();
      //Info->PrintF("", "Regions:\n");
      for( UInt iReg = 0; iReg  < regions.GetSize(); ++iReg ) {
       std::string regionName = regions[iReg]->Get("name")->As<std::string>();
       RegionIdType regionId = ptGrid->GetRegion().Parse(regionName);
       actPart.regions.Push_back(regionId);
       //Info->PrintF("", "\t%s\n",regionName.c_str());
      }

      // read direction
      PtrParamNode dirNode = actPartNode->Get("direction");

      actPart.orientFlag = dirNode->Get("orientation")->As<Integer>();
      //Info->PrintF("", "Orientation:\t\t%d\n",actPart.orientFlag);

      // read uniformity of current density
      //actPart.uniformCurrentDens_ = dirNode->Get("uniformCurrentDensity")->As<bool>();
      //std::string conducModel = (actPart.uniformCurrentDens_) ? "stranded" : "solid";
      //Info->PrintF("", "Conductor model:\t%s\n",conducModel.c_str());

      // Check if current direction is given analytical or in the form of
      // automatic current direction calculation
      if( dirNode->Has("analytic")) {

        PtrParamNode alNode = dirNode->Get("analytic");

        // Note: in case of a 2D coil, we only have direction in phi (axi) or z (plane)
        if( ptGrid_->GetDim() == 2 ) {
          StdVector<std::string> dirReal(1), dirImag(1);
          // it is a vector so that the code for generation is the same in the PDE for 2D and 3D
          // but it has only one component
          dirReal[0] = "1";
          dirImag[0] = "0";
          if ( actPart.orientFlag < 0 ) {
            dirReal[0] = "-1";
          }
          actPart.jUnitVec = CoefFunction::Generate(mParser_, type, dirReal, dirImag );
          
        } else {

          // Read in coordinate system
          std::string coordSysId = "default";
          alNode->GetValue("coordSysId", coordSysId, ParamNode::PASS);
          CoordSystem *  coordSys = domain->GetCoordSystem(coordSysId);

          // Read in components
          ParamNodeList compList = alNode->GetList("comp");
          UInt dim = ptGrid_->GetDim();
          StdVector<std::string> real(dim), imag(dim);
          imag.Init("0.0");
          real.Init("0.0");
          std::string dof, val;
          for( UInt j = 0; j < compList.GetSize(); ++j  ) {
            compList[j]->GetValue("dof", dof);
            compList[j]->GetValue("value", val, ParamNode::PASS );
            Integer index = coordSys->GetVecComponent(dof)-1;
            real[index] = val;
          }
          // Generate coil coefficient function for current direction
          PtrCoefFct dirCoef;
          if (alNode->Get("normalise",ParamNode::PASS)->As<std::string>() == "yes") { // normalise direction by it's length
            PtrCoefFct inputDirCoef = CoefFunction::Generate(mParser_, type, real, imag );
            dirCoef = CoefFunction::Generate(mParser_, type, CoefXprVecScalOp( mParser_, inputDirCoef, CoefXprUnaryOp( mParser_, inputDirCoef, CoefXpr::OP_NORM ) , CoefXpr::OP_DIV ) );
          }
          else { // keep it as it was defined
            dirCoef = CoefFunction::Generate(mParser_, type, real, imag );
          }

          // in the end multiply by the orientation factor
          CoefXprVecScalOp orientOp = CoefXprVecScalOp( mParser_, dirCoef,
                                                        boost::lexical_cast<std::string>(actPart.orientFlag),
                                                        CoefXpr::OP_MULT );

          actPart.jUnitVec = CoefFunction::Generate(mParser_, type, orientOp );


          if( coordSysId != "default" ) {
            actPart.jUnitVec->SetCoordinateSystem( coordSys );
          }
        }
        
      } else if(dirNode->Has("automatic")) {

        EXCEPTION("Automatic determination of coil current density not implemented.");

      } else if(dirNode->Has("external")) {

        // The idea is to use the ReadUserFieldValues of the PDE because it provides everything we need.
        // However, the PDE is impossible to #include, so the PDE has to take care of the rest.
        PtrParamNode extNode = dirNode->Get("external");
        partsExtJDir_[partPtr] = extNode;

      } else {
        EXCEPTION("Only analytic or automatic current directions are allowed");
      }

      /* code from NACS
      // winding surfaces are only sensible in 3D
      if( ptGrid->GetDim() == 3 ) {

        // ---------------
        //  3D Case
        // ---------------

        // sum up areas for input and output surfaces
        Double inArea = 0.0, outArea = 0.0;

        // === INFLOW SURFACES ===
        ParamNodeList inSurfNodes = dirNode->Get("inputSurfaceList")->GetList("surface");
        for( UInt iSurf = 0; iSurf < inSurfNodes.GetSize(); ++iSurf ) {
          std::string surfName = inSurfNodes[iSurf]->Get("name")->AsString();

          // check if surface was already provided
          if( actPart.inputSurfRegions.Find(surfName) != -1  ) {
            EXCEPTION("Surface '" << surfName
                      << "' was already given as inflow surface for one part of coil '"
                      << coilId_ << "'." );
          }

          // check if all elements of the surface are neighboring a volume element
          // of this part
          StdVector<UInt> missingSurfelemNums;
          StdVector<RegionIdType> empty;
          shared_ptr<EntityList> actList =
                        ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, surfName );
          ptGrid_->CheckSurfNeighbors(actList, actPart.regions, empty,
                                      missingSurfelemNums, false );
          if( missingSurfelemNums.GetSize() > 0 ) {
            EXCEPTION("Not all surface elements of '" << surfName
                      << "' have a coil volume element neighbor for coil '"
                      << coilId_ << "'. Please check the coil definition!");

          }

          actPart.inputSurfRegions.Push_back(surfName);
          // Note: We always have to consider the "surface area" of the input
          // output surfaces. Thus we are never allowed to calculate the depth
          // of the setup, i.e. isAxi is always false
          bool isAxi = false;
          inArea += ptGrid_->CalcVolumeOfElems(actList, isAxi, false );
        }

        // === OUTFLOW SURFACES ===
        ParamNodeList outSurfNodes = dirNode->Get("outputSurfaceList")->GetList("surface");
        for( UInt iSurf = 0; iSurf < inSurfNodes.GetSize(); ++iSurf ) {
          std::string surfName = outSurfNodes[iSurf]->Get("name")->AsString();
          actPart.outputSurfRegions.Push_back(surfName);

          if( actPart.inputSurfRegions.Find(surfName) != -1 ) {
            EXCEPTION("Surface '" << surfName
                      << "' was already defined as outflow surface one part of coil '"
                      << coilId_ << "'." );
          }

          // security check: assure, that none of the output surfaces is contained in the
          // list of input surfaces
          if( actPart.inputSurfRegions.Find(surfName) != -1 ) {
            EXCEPTION("The surface region '" << surfName << "' is already used for "
                      "specifying the inflow current direction of one part of coil '"
                      << coilId_ << "'. Inflow and outflow surfaces have to be distinct!");
          }

          // check if all elements of the surface are neighboring a volume element
          // of this part
          StdVector<UInt> missingSurfelemNums;
          StdVector<RegionIdType> empty;
          shared_ptr<EntityList> actList =
              ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, surfName );
          ptGrid_->CheckSurfNeighbors(actList, actPart.regions, empty,
                                      missingSurfelemNums, false );
          if( missingSurfelemNums.GetSize() > 0 ) {
            EXCEPTION("Not all surface elements of '" << surfName
                      << "' have a coil volume element neighbor for coil '"
                      << coilId_ << "'. Please check the coil definition!");

          }
          // Note: We always have to consider the "surface area" of the input
          // output surfaces. Thus we are never allowed to calculate the depth
          // of the setup, i.e. isAxi is always false
          bool isAxi = false;
          outArea += ptGrid_->CalcVolumeOfElems(actList, isAxi, false );
        }

        // ensure that input and output areas have (almost) the same area
        if( std::abs(inArea - outArea) > EPS ) {
          std::stringstream w;
          w << "Inflow and outflow surface area of part " << iPart+1 << " of coil '"
              << coilId_ << "' have different area:\n";
          w << "\tinput: " << inArea << std::endl
          << "\toutput: " << outArea << std::endl;
          Warning(w.str().c_str());
        }

        actPart.coilCrossSect  = inArea;

        // insert correct values (if not a measurement coil )
        if( value_ == "0.0") {
          actPart.sourceVal = "0.0";
        } else {
          // Note: in 3D we have always a free coordiante system associated with
          // the coil, so we have just one DOF to set, which always the first one.
          actPart.sourceVal = "("+ value_ + "*" + GenStr(actPart.orientFlag ) + ")";
        }

        // calculate current direction using free coordinate system
        actPart.coordSys = SetupCoosy( iPart+1,
                                       actPart,
                                       actPart.regions,
                                       actPart.inputSurfRegions,
                                       actPart.outputSurfRegions );

      } else {
        // ---------------
        //  2D Case
        // ---------------
        // We just take all the regions as "surfaces" and calculate their
        // "volume" i.e. area * 1m

        for( UInt iReg = 0; iReg < actPart.regions.GetSize(); ++iReg) {
          std::string regionName =
              ptGrid_->RegionIdToName(actPart.regions[iReg]);
          shared_ptr<EntityList> actList =
              ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
          // Note: We always have to consider the "surface area" of the input
          // output surfaces. Thus we are never allowed to calculate the depth
          // of the setup, i.e. isAxi is always false
          bool isAxi = false;
          actPart.coilCrossSect +=
              ptGrid_->CalcVolumeOfElems(actList, isAxi, false );
        }

        // in case of measurement coil set value to zero
        if( value_ == "0.0") {
          actPart.sourceVal = "0.0";
        } else {
          // insert value_ * directionFlag into value vector at z-position
          //
          actPart.sourceVal ="("+ value_ + "*" + GenStr(actPart.orientFlag ) + ")";
        }

        // in 2D, we always have the global Cartesian coordinate system
        actPart.coordSys.reset(new DefaultCoordSystem(ptGrid_, false));
      }*/

      // ========================================================
      //  Intiialize setup of coil (geometry, current direction)
      // ========================================================
      // initialize geometry setup (turns, cross section, orientation) 
      // a) wire cross section / kappa -> evaluate number of turns
      if( actPartNode->Has("wireCrossSection") ) {
        std::string areaStr = actPartNode->Get("wireCrossSection")->Get("area")->As<std::string>();
        MathParser::HandleType handle = mParser_->GetNewHandle();
        mParser_->SetExpr(handle,areaStr);
        actPart.wireCrossSect = mParser_->Eval(handle);
        //actPart.fillFactor =
        //    actPartNode->Get("wireCrossSection")->Get("fillFactor")->As<Double>();
        
        //actPart.numTurns = UInt((actPart.coilCrossSect * actPart.fillFactor) 
        //                        / actPart.wireCrossSect );
      } else 

      // b) turns / kappa given -> determine wire crossSection
      if( actPartNode->Has("windingTurns") ) {
        EXCEPTION("Currently the windingTurns can not be specified")
        /*actPart.numTurns =
            actPartNode->Get("windingTurns")->Get("number")->AsDouble();
        actPart.fillFactor =
            actPartNode->Get("windingTurns")->Get("fillFactor")->AsDouble();
        // switch type of
        actPart.wireCrossSect = (actPart.coilCrossSect * actPart.fillFactor)
                                / Double(actPart.numTurns);*/
      } else {
        EXCEPTION( "Either the wireCrossSection or the number of turns has to be "
                   "specified." );
      }

      /* code from NACS
      // ==========================================
      // CHECK FOR SOLID / STRANDED COIL
      // ==========================================
      // in case of a solid conductor model, we effectively have only one
      // turn and thus set hard-coded the number of turns to 1 and the cross-
      // section of the wire to the cross-section of the whole coil
      if( !actPart.uniformCurrentDens_ ) {
        actPart.numTurns = 1;
        actPart.wireCrossSect = actPart.coilCrossSect;
      }*/

      //Info->PrintF("","Number of Turns:\t%d\n",actPart.numTurns);
      //Info->PrintF("","Fill Factor:\t\t%g\n",actPart.fillFactor);
      //Info->PrintF("","Wire Cross Section:\t%g\n",actPart.wireCrossSect);
      //Info->PrintF("","Input Cross Section:\t%g\n",actPart.coilCrossSect);

      // read coil resistance (if given)
      if( actPartNode->Has("resistance") ) {
        actPartNode->Get("resistance")->GetValue("value",actPart.resistance);
      //Info->PrintF("","Coil Resistance:\t%s\n",actPart.resistance.c_str());
      }

      // In the end, loop over all regions and add part to every region
      for( UInt i = 0; i < actPart.regions.GetSize(); ++i ) {
        parts_[actPart.regions[i]] = partPtr;
      }
      //Info->PrintF("","\n");
    } // loop over parts
  }
  
  /*Coil::Coil( const std::string& id, Grid * ptGrid ) {
    EXCEPTION("Implement me");
    ptGrid_ = ptGrid;
    coilId_ = id;

    // myParam_ is not set

    // initialize data members
    isAxi_ = false;
    sourceType_ = NO_SOURCE_TYPE;
    coilId_ = "";
    value_ = "0.0";
    phase_ = "0.0";

    // obtain math parser handle
    mParser_ = domain->GetMathParser();

    valueHandle_ = MathParser::GLOB_HANDLER;
    phaseHandle_ = MathParser::GLOB_HANDLER;
  }*/
  
  /*shared_ptr<CoordSystem>
  Coil:: SetupCoosy( UInt iPart,
                     const Part& actPart,
                     StdVector<RegionIdType>& regions,
                     const StdVector<std::string>& inSurfaces,
                     const StdVector<std::string>& outSurfaces ) {
//
//       <free id="coil-f">
//              <regionList>
//                <region name="coil"/>
//              </regionList>
//              <dirName> r </dirName>
//              <surf1> coil-back </surf1>
//              <surf2> coil-front </surf2>
//            </free>
//
//       
//       // calculate coordinate system
//    std::string idString = "_coil_"+coilId_+"_part_"+GenStr(iPart);
//    PtrParamNode root ( new ParamNode("free", ParamNode::ELEM) );
//    root->AddChildAttribute("id",idString, true );
//    
//    // Important: pass flag, if coil current should get normalized to one!
//    if( actPart.uniformCurrentDens_ ) {
//      root->AddChildAttribute("normalize","true", true );
//    } else {
//      root->AddChildAttribute("normalize","false", true );
//    }
//    
//    PtrParamNode rList (new ParamNode("regionList", ParamNode::ELEM ) );
//    root->AddChild(rList, true);
//    // add all regions
//    for( UInt i = 0; i < regions.GetSize(); ++i ) {
//      PtrParamNode rNode(new ParamNode("region", ParamNode::ELEM));
//      rNode->AddChildAttribute("name",ptGrid_->RegionIdToName(regions[i]), true );
//      rList->AddChild( rNode, false );
//    }
//    root->AddChildElem("dirName", "r", true );
//    
//    // WARNING: Currently, the free coordinate system only supports
//    // 1 surface for inflow/outflow. 
//    // Thus, we issue a warning, if the user has provided several interfaces
//    if( inSurfaces.GetSize() > 1 || outSurfaces.GetSize() > 1 ) {
//      EXCEPTION("Currently, NACS supports just 1 surfaces for in/out flow of "
//                "coil current");
//    }
//    root->AddChildElem("surf1", inSurfaces[0], true );
//    root->AddChildElem("surf2", outSurfaces[0], true );
//    
//    // debugging: dump tree
//    //root->Dump();
//    
//    
//    shared_ptr<FreeCoordSystem> cSys;
//    try {
//     cSys.reset(new FreeCoordSystem(idString, ptGrid_, root, false, false) );
//    } catch (Exception& e) {
//      RETHROW_EXCEPTION(e, "Can not determine current direction in coil '" << coilId_ << "'");
//    }
//    return cSys;
    return shared_ptr<CoordSystem>();
  }*/
  
  shared_ptr<EntityList> Coil::GetElems() {
    shared_ptr<ElemList> elems;
    elems.reset( new ElemList( ptGrid_ ) );
    std::map<RegionIdType, shared_ptr<Coil::Part> >::iterator partIt =
        parts_.begin();
    while( partIt != parts_.end() ){
      StdVector<Elem*> elemsReg;
      ptGrid_->GetElems( elemsReg, partIt->first );
      elems->AddElements( elemsReg );
      partIt++;
    }
    return elems;
  }
  
  Coil::Part::Part() {
    orientFlag = 0;
    resistance = "0.0";
    numTurns = 0;
    //coilCrossSect = 0.0;
    wireCrossSect = 0.0;
    //fillFactor = 1.0;
    //uniformCurrentDens_ = true;
    //hasAnalyticalDir_ = true;
  }
 
  Coil::Part::~Part() {
  }
  

}
