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
#include "Utils/mathParser/mathParser.hh"
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
    coilOptimization_ = false;

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
      } else if ( exType == "specialvoltage" ) {
        sourceType_ = SPECIALVOLTAGE;
      } else if ( exType == "specialcurrent" ) {
        sourceType_ = SPECIALCURRENT;
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

      std::string exType = myParam_->Get("sourceMultiharmonic")->Get("type")->As<std::string>();
      if ( exType == "current" ) {
        sourceType_ = CURRENT_MULTHARM;
      } else if ( exType == "specialvoltage" ) {
        sourceType_ = SPECIALVOLTAGE_MULTHARM;
      } else{ EXCEPTION("The kind of excitation does not exist in the coil excitation!"); }

      std::string value, phase;
      ParamNodeList harmonicList = myParam_->Get("sourceMultiharmonic")->GetList("harmonic");
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

    // ============================
    //  Loop over parts
    // ============================
    ParamNodeList parts = myParam_->GetList("part");
    for( UInt iPart = 0; iPart < parts.GetSize(); ++iPart ) {
     

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
      }

      // read direction
      PtrParamNode dirNode = actPartNode->Get("direction");

      actPart.orientFlag = dirNode->Get("orientation")->As<Integer>();

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

        // check if coil-part is used for special A-V formulation
        if( actPartNode->Has("gradV_in_gradV") ){
          PtrParamNode extNodegV = actPartNode->Get("gradV_in_gradV/external");
          partsExtIntgVgV_[partPtr] = extNodegV;
        }

      } else {
        EXCEPTION("Only analytic or automatic current directions are allowed");
      }


      // ========================================================
      //  Intiialize setup of coil (geometry, current direction)
      // ========================================================
      // initialize geometry setup (turns, cross section, orientation) 
      // a) wire cross section / kappa -> evaluate number of turns
      if(actPartNode->Has("wireCrossSection"))
      {
        std::string areaStr = actPartNode->Get("wireCrossSection/area")->As<std::string>();
        unsigned int handle = mParser_->GetNewHandle();
        mParser_->SetExpr(handle,areaStr);
        actPart.wireCrossSect = mParser_->Eval(handle);
        //actPart.fillFactor =
        //    actPartNode->Get("wireCrossSection")->Get("fillFactor")->As<Double>();
        
        //actPart.numTurns = UInt((actPart.coilCrossSect * actPart.fillFactor) 
        //                        / actPart.wireCrossSect );
        coilOptimization_ = actPartNode->Get("wireCrossSection/coil_top_opt")->As<bool>();
      } else 

      // b) turns / kappa given -> determine wire crossSection
      if( actPartNode->Has("windingTurns") ) {
        EXCEPTION("Currently the windingTurns can not be specified")
      } else {
        EXCEPTION( "Either the wireCrossSection or the number of turns has to be "
                   "specified." );
      }


      // read coil resistance (if given)
      if( actPartNode->Has("resistance") ) {
        actPartNode->Get("resistance")->GetValue("value",actPart.resistance);
      }

      // In the end, loop over all regions and add part to every region
      for( UInt i = 0; i < actPart.regions.GetSize(); ++i ) {
        parts_[actPart.regions[i]] = partPtr;
      }
    } // loop over parts
  }
  
  
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
    wireCrossSect = 0.0;
  }
 
  Coil::Part::~Part() {
  }
  

}
