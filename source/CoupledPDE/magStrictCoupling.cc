// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <map>
#include <string>
#include <utility>

#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/entityList.hh"
#include "Domain/resultInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Forms/baseForm.hh"
// integrator (bi-)linear forms
#include "Forms/linMagStrictInt.hh"
#include "Forms/linearForm.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/magneticScalarPDE.hh"
#include "PDE/mechPDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/result.hh"
#include "magStrictCoupling.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
class BiotSavart;
class BaseMaterial;

  MagStrictCoupling::MagStrictCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                        PtrParamNode paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode, "magnetostriction") {

    materialClass_ = MAGNETOSTRICTIVE;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );
  }


  // **************
  //   Destructor
  // **************
  MagStrictCoupling::~MagStrictCoupling() {
  }


  // ***************
  //   CalcResults
  // ***************
  void MagStrictCoupling::CalcResults( shared_ptr<BaseResult> result ) {

    switch (result->GetResultInfo()->resultType ) {
      case MAG_FLUX_DENSITY:
        if ( isComplex_ ) {
          CalcBField<Complex>( result );
        } else{
          CalcBField<Double>( result );
        }
        break;

      case MECH_STRESS:
        if ( isComplex_ ) {
          CalcStress<Complex>( result );
        } else {
          CalcStress<Double>( result );

        }
        break;

      default:
        WARN("Result type not computable by magnmechanic coupling");
    }
  }

  template<class TYPE>
  void MagStrictCoupling::CalcBField( shared_ptr<BaseResult> result ) {
    EXCEPTION("Not yet implemented. Merge from implementation of Ni Yunyun" );
  }
  
  template<class TYPE>
  void MagStrictCoupling::CalcStress( shared_ptr<BaseResult> result ) {
    EXCEPTION("Not yet implemented. Merge from implementation of Ni Yunyun") ;
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void MagStrictCoupling::DefineIntegrators() {

    
    // perform explicit cast of PDEs
    MechPDE *mechPDE = dynamic_cast<MechPDE*>(pde1_);
    MagScalarPDE * magPDE = dynamic_cast<MagScalarPDE*>(pde2_);
    
    Global::ComplexPart matType = Global::REAL;
    RegionIdType actRegion;
    // BaseMaterial * actSDMat = NULL; // TODO: Unused variable actSDMat

    // get material from magnetostatics
    std::map<RegionIdType, BaseMaterial*> magMat =
      pde2_->getPDEMaterialData();


    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      // actSDMat = it->second;
      matType = Global::REAL;

      // determine the tenstor type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // ==== Define Stiffness integrator K_U\Phi =======
      BiLinFormContext *actContextStiff = NULL;
      BaseForm *bilinearStiff = new LinMagStrictInt(materials_[actRegion], tensorType);
      bilinearStiff->SetMatDataType( matType );
      actContextStiff = new BiLinFormContext( bilinearStiff, STIFFNESS  );


      // We also need to set the transposed of the coupling
      // matrix to the lower diagonal side
      actContextStiff->SetCounterPart( true );
      actContextStiff->SetEntryType( matType );
      actContextStiff->SetPtPdes( mechPDE, magPDE );  // PDE1: MECH, PDE2: Mag)
      actContextStiff->SetResults( results1_[0], results2_[0],
                                   actSDList, actSDList );
      assemble_->AddBiLinearForm( actContextStiff );
      
      // ==== Add BiotSavart- Mech Coupling RHS === 
      // Only add integ
      shared_ptr<BiotSavart> bs =  magPDE->GetBiotSavart(); 
      if( bs ) {
        BiotSavartMechCouplingInt * bsInt = 
            new BiotSavartMechCouplingInt( materials_[actRegion], tensorType,
                                           bs, isaxi_ );
        LinearFormContext * linRhs = new LinearFormContext( bsInt);
        linRhs->SetPtPde( pde1_ );//mechanical
        linRhs->SetResult( results1_[0], actSDList );
        assemble_->AddLinearForm( linRhs );
      }
      
    }
  }


  void MagStrictCoupling::DefineAvailResults() {

   // Check for subType
    StdVector<std::string> stressDofNames, vecDofNames;
    if( subType_ == "3d" ) {
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
      vecDofNames = "x", "y", "z";

    } else if( subType_ == "planeStrain" ) {
      stressDofNames = "xx", "yy", "xy";
      vecDofNames = "x", "y";

    } else if( subType_ == "planeStress" ) {
      stressDofNames = "xx", "yy", "xy";
      vecDofNames = "x", "y";

    } else if( subType_ == "axi" ) {
      stressDofNames = "rr", "zz", "rz", "phiphi";
      vecDofNames = "r", "z";

    }

    // === MECHANIC STRESS ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressDofNames;
    stress->unit = "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    stress->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( stress );

    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX_DENSITY;
    flux->dofNames = vecDofNames;
    flux->unit = "Vs/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    flux->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( flux );
  }
}
