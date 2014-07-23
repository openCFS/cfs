// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MagneticPDE.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/Coil.hh"

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefXpr.hh"


#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

namespace CoupledField {

DECLARE_LOG(magpde)
DEFINE_LOG(magpde, "magpde")

MagneticPDE::MagneticPDE(Grid * aptgrid, PtrParamNode paramNode,
                         PtrParamNode infoNode,
                         shared_ptr<SimState> simState, Domain* domain)
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

  // =====================================================================
  // set solution information
  // =====================================================================
  pdename_          = "magnetic";
  pdematerialclass_ = ELECTROMAGNETIC;
  
  //! Always use updated Lagrangian formulation 
  updatedGeo_        = true;
  
  reluc_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_, dim_, isComplex_));

  // determine if there are coils excited by voltage
  hasVoltCoils_ = false;
  PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
  if ( coilNode ){
    ParamNodeList coilNodes = coilNode->GetChildren();
    for( UInt k = 0; k < coilNodes.GetSize(); k++ ){
      if( coilNodes[k]->Has("source") ){
        std::string exType = coilNodes[k]->Get("source")->Get("type")->As<std::string>();
        if( exType == "voltage" ){
          hasVoltCoils_ = true;
          break;
        }
      }
    }
  }

}

  MagneticPDE::~MagneticPDE()
  {

  }

  shared_ptr<Coil> MagneticPDE::GetCoilById(const Coil::IdType& id) {
    return coils_.at(id);
  }

  void MagneticPDE::InitNonLin()
  {

    SinglePDE::InitNonLin();
  }


  void MagneticPDE::DefineIntegrators() {

    // determine tensor representation of the material parameters needed
    SubTensorType tensorType;
    if ( dim_ == 3 ) {
      tensorType = FULL;
    } else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      } else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }
    
    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    
    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      
      // Set current region and material
      actRegion = it->first;
      actMat = it->second;

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list and add it fefunction
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      myFct->AddEntityList( actSDList );
      
      // --------------------------
      //  Set region approximation
      // --------------------------
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")
          ->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      
      //get possible nonlinearities defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
      if ( nonLinTypes.GetSize() > 0 ) {
        if ( nonLinTypes.Find(HYSTERESIS) != -1 ) {
          Exception("Hysteresis nonlinearity is not supported.");
        }else if (  nonLinTypes.Find(PERMEABILITY) != -1 ) {
          PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
          PtrCoefFct nuNl =
              actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY, Global::REAL, magFluxCoef);

          BaseBDBInt * stiffInt = NULL;
          if( dim_ == 2) {
            if( isaxi_ ) {
              // axisymmetric case
              stiffInt = new BBInt<>(new CurlOperatorAxi<Double>(), nuNl, 1.0, updatedGeo_);
            } else {
              // plane 2D case
              stiffInt = new BBInt<>(new CurlOperator<FeH1,2,Double>(), nuNl, 1.0, updatedGeo_);
            }
          } else {
            // 3D case
            stiffInt = new BBInt<>(new CurlOperator<FeH1,3,Double>(), nuNl, 1.0, updatedGeo_);

          }
          stiffInt->SetName("CurlCurlIntegrator-NL");
          BiLinFormContext * stiffContext =
              new BiLinFormContext(stiffInt, STIFFNESS );
          stiffContext->SetEntities( actSDList, actSDList );
          stiffContext->SetFeFunctions( myFct, myFct );
          assemble_->AddBiLinearForm( stiffContext );
          // Important: Add bdb-integrator to global list, as we need them later
          // for calculation of postprocessing results
          bdbInts_[actRegion] = stiffInt;
          // add also material to global, distributed reluctivity coefficient function
          //reluc_->AddRegion(actRegion, nuNl);

          // ================================================
          //  Nonlinear Stiffness Integrator (only Newton )
          // ================================================
          // Note: currently we set the nonlinear method hard-coded to NEWTON for
          // testing purpose
          if( nonLinMethod_ == NEWTON ) {
            PtrCoefFct nuDeriv = actMat->GetTensorCoefFncNonLin( MAG_RELUCTIVITY_DERIV, tensorType,
                                                                 Global::REAL, magFluxCoef );

            //create stiffness integrator
            BiLinearForm* stiff2 = NULL;
            //stiff2 = new BDBInt<>(new CurlOperator<FeHCurl,3, Double>(), nuDeriv, 1.0, updatedGeo_) ;
            if( dim_ == 2) {
              if( isaxi_ ) {
                // axisymmetric case
                stiff2 = new BDBInt<>(new CurlOperatorAxi<Double>(), nuDeriv, 1.0, updatedGeo_);
              } else {
                // plane 2D case
                stiff2 = new BDBInt<>(new CurlOperator<FeH1,2,Double>(), nuDeriv, 1.0, updatedGeo_);
              }
            } else {
              // 3D case
              stiff2 = new BDBInt<>(new CurlOperator<FeH1,3,Double>(), nuDeriv, 1.0, updatedGeo_);

            }

            stiff2->SetName("CurlCurlIntegrator-NL-Newton");
            //! mark the bi-linear form to be a Newton part
            stiff2->SetNewtonBilinearForm();

            BiLinFormContext * stiffContext2 =
                new BiLinFormContext(stiff2, STIFFNESS );
            stiffContext2->SetEntities( actSDList, actSDList );
            stiffContext2->SetFeFunctions( myFct, myFct );
            assemble_->AddBiLinearForm( stiffContext2 );
          }
        }
      }else{
        // ====================================================================
        //  Standard Linear CASE (2D AND 3D)
        // ====================================================================

          BaseBDBInt * stiffInt = NULL;
          PtrCoefFct curCoef =
              actMat->GetTensorCoefFnc( MAG_RELUCTIVITY, tensorType,
                                       Global::REAL );
          if( dim_ == 2) {
            if( isaxi_ ) {
              // axisymmetric case
              stiffInt = new BDBInt<>(new CurlOperatorAxi<Double>(), curCoef, 1.0, updatedGeo_);
            } else {
              // plane 2D case
              stiffInt = new BDBInt<>(new CurlOperator<FeH1,2,Double>(), curCoef, 1.0, updatedGeo_);
            }
          } else {
            // 3D case
            stiffInt = new BDBInt<>(new CurlOperator<FeH1,3,Double>(), curCoef, 1.0, updatedGeo_);
          }
          stiffInt->SetName("CurlCurlIntegrator");
          stiffInt->SetFeSpace( mySpace);
          BiLinFormContext * stiffIntDescr =
              new BiLinFormContext(stiffInt, STIFFNESS );
          stiffIntDescr->SetEntities( actSDList, actSDList );
          stiffIntDescr->SetFeFunctions( myFct, myFct );

          assemble_->AddBiLinearForm( stiffIntDescr );

          // Important: Add bdb-integrator to global list, as we need them later
          // for calculation of postprocessing results
          bdbInts_[actRegion] = stiffInt;

          // add also material to global, distributed reluctivity coefficient function
          reluc_->AddRegion(actRegion, curCoef);

          // ====================================================================
          //  3D CASE (additional stiffness integrator)
          // ====================================================================
          if( dim_ == 3 ) {
            BaseBDBInt * divInt =
                new BDBInt<>(new DivOperator<FeH1,3,Double>(), curCoef, 1.0, updatedGeo_);
            divInt->SetFeSpace( mySpace );
            divInt->SetName("DivDivIntegrator");
            BiLinFormContext * divIntDescr =
                        new BiLinFormContext(divInt, STIFFNESS );
            divIntDescr->SetEntities( actSDList, actSDList );
            divIntDescr->SetFeFunctions( myFct, myFct );
            assemble_->AddBiLinearForm( divIntDescr );
          }
      }
        
        
        // ====================================================================
        //  MASS - MATRIX (all dimensions)
        // ====================================================================
        // Note: here we have currently the problem, that we can only treat
        // constant / double valued conductivity values, as we must check for 
        // non-zero. However, how do we do this in case of non-constant parameters?
        Double conductivity;
        
        materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
        PtrCoefFct conducCoef =
            materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY,Global::REAL);
//                                         lexical_cast<std::string>(conductivity));
        {
          BaseBDBInt *massInt = NULL;
          if( dim_ == 2 ) {
            massInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), conducCoef,1.0, updatedGeo_);
          } else {
            massInt = new BBInt<>(new IdentityOperator<FeH1,3,3>(), conducCoef,1.0, updatedGeo_ );
          }
          massInt->SetName("MassIntegrator");
          BiLinFormContext * massContext = new BiLinFormContext(massInt, DAMPING );
          massContext->SetEntities( actSDList, actSDList );
          massContext->SetFeFunctions( myFct, myFct );
          assemble_->AddBiLinearForm( massContext );
          
          // insert mass integrator to list of defined mass integrators
          massInts_[actRegion] = massInt;
        }
        
        // ====================================================================
        //  3D MIXED CASE (coupling between vector and scalar potential)
        // ====================================================================
        if( isMixed_ && conductivity > EPS ) {
          shared_ptr<BaseFeFunction> potFct = feFunctions_[ELEC_POTENTIAL];
          shared_ptr<FeSpace> potSpace = potFct->GetFeSpace();
          // Important: set region approximation also at space for scalar 
          // potential
          potSpace->SetRegionApproximation(actRegion, polyId,integId);
          
          // add region to feFunction for scalar potential
          potFct->AddEntityList( actSDList );
          
          // ------------------------------
          // diagonal mass matrix M_phiphi
          // ------------------------------
          {
            BiLinearForm * phiDivInt = 
                new BBInt<>(new GradientOperator<FeH1,3,1,Double>(), conducCoef, 1.0, updatedGeo_);
            phiDivInt->SetName("MassIntegrator_PhiPhi");
            BiLinFormContext * massContext = 
                new BiLinFormContext(phiDivInt, DAMPING);
            massContext->SetEntities( actSDList, actSDList );
            massContext->SetFeFunctions( potFct, potFct );
            assemble_->AddBiLinearForm( massContext );
          }
          
          // ------------------------------
          // off-diagonal mass matrix M_A_phi
          // ------------------------------
          {
            BiLinearForm * cplInt = 
                new ABInt<>(new IdentityOperator<FeH1,3,3,Double>() ,
                            new GradientOperator<FeH1,3,1,Double>(), conducCoef, 1.0, updatedGeo_);
            cplInt->SetName("MassIntegrator_Coupling_Phi_A");
            BiLinFormContext * cplContext = 
                new BiLinFormContext(cplInt, DAMPING );
            cplContext->SetCounterPart(true);
            cplContext->SetEntities( actSDList, actSDList );
            cplContext->SetFeFunctions( myFct, potFct );
            assemble_->AddBiLinearForm( cplContext );
          }
        } // mixed
        
       // linear part
    } // regions
    
    // ============================
    // COIL INTEGRATORS
    // ============================
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    std::map<Coil::IdType, shared_ptr<Coil> >::iterator coilIt;
    coilIt = coils_.begin();
    for( ; coilIt != coils_.end(); coilIt++ ) {
      Coil& actCoil = *(coilIt->second);
      // run over all parts
      std::map<RegionIdType,shared_ptr<Coil::Part> >::iterator partIt;
      partIt = actCoil.parts_.begin();
      if( actCoil.sourceType_ == Coil::CURRENT ) {
        /*
        ============================================
         1) CURRENT driven coils

         Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
              Sens. and Act., 2nd edition, p. 131ff
        ============================================
        */
        for( partIt = actCoil.parts_.begin(); 
            partIt != actCoil.parts_.end(); 
            partIt++ ) {
          Coil::Part & actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
          actSDList->SetRegion( actRegion );
          LinearForm* curInt = NULL;

         
            // generate source current vector
            CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_,
                                                     CoefXpr::OP_MULT);
            PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

            CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), 
                                                     CoefXpr::OP_DIV);
            PtrCoefFct jFct = CoefFunction::Generate(mp_, part, jVec);
            
            CoefFunction * tmp = jFct.get();
            std::cerr << "jFct is " << tmp->ToString() << std::endl;
            // ===========
            //  3D CASE
            // ===========
            if( dim_ == 3 ) {
            if( isComplex_ ) {
              curInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,3,3,Complex>(),
                                                  1.0, jFct, updatedGeo_);
            }
            else {
              curInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3,Double>(),
                                                 1.0, jFct, updatedGeo_);
            }
            
            coilCurrentDens_[actRegion] = jFct;
          } else {
            // ===============
            //  2D / AXI CASE
            // ===============
//            CoefXprBinOp jVec = CoefXprBinOp(mp_,  actCoil.srcVal_,
//                                             boost::lexical_cast<std::string>(actPart.wireCrossSect), 
//                                             CoefXpr::OP_DIV);
//            PtrCoefFct jScal = CoefFunction::Generate(mp_, part, jVec);
//            
//            // Generate 1-component vector function
//            PtrCoefFct jFct = CoefFunction::Generate(mp_, part, )
            
            coilCurrentDens_[actRegion] = jFct;
            if( isComplex_ ) {
              curInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1>(),
                                                  1.0, jFct, updatedGeo_);
            } else {
              curInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1>(),
                                                 1.0, jFct, updatedGeo_);
            }
          }


          curInt->SetName("CoilIntegrator");
          LinearFormContext * coilContext =
              new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( myFct );
          assemble_->AddLinearForm( coilContext );
          // obtain coefficient function
        } // loop: parts
      } else {

        /*
        ============================================
         2) VOLTAGE driven coils

         Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
              Sens. and Act., 2nd edition, p. 211ff
        ============================================
         The coupled equation system in this case looks like
           ( M_A     0 ) ( A_dot ) + ( K_A -f_A  ) ( A ) = ( 0 )
           ( (f_A)^T 0 ) ( i_dot )   ( 0     R   ) ( i )   ( u )
        */

        std::string totRstr = "";
        shared_ptr<CoilList> singleCoilList( new CoilList( ptGrid_ ) );
        singleCoilList->AddCoil( coilIt->second );
        feFunctions_[COIL_CURRENT]->AddEntityList( singleCoilList );

        for( partIt = actCoil.parts_.begin();
             partIt != actCoil.parts_.end();
             partIt++ ) {

          Coil::Part & actPart = *(partIt->second);

          if( totRstr.empty() ){
            totRstr = actPart.resistance;
          } else {
            totRstr += " + " + actPart.resistance;
          }

          CoefXprVecScalOp eJscaledOp = CoefXprVecScalOp( mp_, actPart.jUnitVec,
              boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV );
          PtrCoefFct eJscaled = CoefFunction::Generate( mp_, part, eJscaledOp );

          shared_ptr<ElemList> actSDList( new ElemList( ptGrid_ ) );
          RegionIdType actRegion = partIt->first;
          actSDList->SetRegion( actRegion );

          // === -f_A ===
          LinearForm* psiDotInt;
          if( dim_ == 3 ) {
            if( isComplex_ ) {
              psiDotInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,3,3,Complex>(),
                  -1.0, eJscaled, updatedGeo_);
            } else {
              psiDotInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3,Double>(),
                  -1.0, eJscaled, updatedGeo_);
            }
          } else {
            if( isComplex_ ) {
              psiDotInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1,Complex>(),
                  -1.0, eJscaled, updatedGeo_);
            } else {
              psiDotInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1,Double>(),
                  -1.0, eJscaled, updatedGeo_);
            }
          }
          psiDotInt->SetName("CoilVoltCouplInt");

          bool assembleTransposed = false;
          BiLinearForm* pseudoBiLin = new BiLinWrappedLinForm( psiDotInt, assembleTransposed );
          BiLinFormContext* voltCoilContext = new BiLinFormContext( pseudoBiLin, STIFFNESS );
          voltCoilContext->SetEntities( actSDList, singleCoilList );
          voltCoilContext->SetFeFunctions( myFct, feFunctions_[COIL_CURRENT] );
          voltCoilContext->SetCounterPart(false);
          assemble_->AddBiLinearForm( voltCoilContext );

          // === (f_A)^T ===
          if( analysistype_ != STATIC ){
            LinearForm* psiDotIntT;
            if( dim_ == 3 ) {
              if( isComplex_ ) {
                psiDotIntT = new BUIntegrator<Complex>( new IdentityOperator<FeH1,3,3,Complex>(),
                    1.0, eJscaled, updatedGeo_);
              } else {
                psiDotIntT = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3,Double>(),
                    1.0, eJscaled, updatedGeo_);
              }
            } else {
              if( isComplex_ ) {
                psiDotIntT = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1,Complex>(),
                    1.0, eJscaled, updatedGeo_);
              } else {
                psiDotIntT = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1,Double>(),
                    1.0, eJscaled, updatedGeo_);
              }
            }
            psiDotIntT->SetName("CoilVoltCouplIntTransposed");

            assembleTransposed = true;
            BiLinearForm* pseudoBiLinT = new BiLinWrappedLinForm( psiDotIntT, assembleTransposed );
            BiLinFormContext* voltCoilContextT = new BiLinFormContext( pseudoBiLinT, DAMPING );
            voltCoilContextT->SetEntities( singleCoilList, actSDList );
            voltCoilContextT->SetFeFunctions( feFunctions_[COIL_CURRENT], myFct );
            voltCoilContextT->SetCounterPart(false);
            assemble_->AddBiLinearForm( voltCoilContextT );
          }

        } // loop: parts

        // === R ===
        PtrCoefFct totR = CoefFunction::Generate( mp_, part, totRstr, "0.0" );
        LinearForm* totRint = new SingleEntryInt( totR );
        totRint->SetName( "CoilResistanceInt" );
        BiLinearForm* totRBiLin = new BiLinWrappedLinForm( totRint, false );
        BiLinFormContext* totRcontext = new BiLinFormContext( totRBiLin, STIFFNESS );
        totRcontext->SetEntities( singleCoilList, singleCoilList );
        totRcontext->SetFeFunctions( feFunctions_[COIL_CURRENT], feFunctions_[COIL_CURRENT] );
        totRcontext->SetCounterPart(false);
        assemble_->AddBiLinearForm( totRcontext );

        // === u ===
        LinearForm* voltInt = new SingleEntryInt( actCoil.srcVal_ );
        voltInt->SetName( "CoilVoltageLoadInt" );
        LinearFormContext* voltContext = new LinearFormContext( voltInt );
        voltContext->SetEntities( singleCoilList );
        voltContext->SetFeFunction( feFunctions_[COIL_CURRENT] );
        assemble_->AddLinearForm( voltContext );
        
      } //if: current / voltage driven
    } // loop: coils
  }
  
  void MagneticPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        if (dim_ == 2) {
          DefineMortarCoupling<2,1>(MAG_POTENTIAL, *ncIt);
          if (isMixed_)
            DefineMortarCoupling<2,1>(ELEC_POTENTIAL, *ncIt);
        }
        else {
          DefineMortarCoupling<3,3>(MAG_POTENTIAL, *ncIt);
          if (isMixed_)
            DefineMortarCoupling<3,1>(ELEC_POTENTIAL, *ncIt);
        }
        break;
      case NC_NITSCHE:
        if (dim_ == 2)
          DefineNitscheCoupling<2,1>(MAG_POTENTIAL, *ncIt );
        else
          EXCEPTION("Currently Nitsche for 3D not implemented");
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

 
  void MagneticPDE::DefineRhsLoadIntegrators() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> mySpace = feFct->GetFeSpace();    
    
    // ==================
    //  COIL INTEGRATORS
    // ==================
    // Loop over all coils
//    for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
//
//      // Set current region and material
//      RegionIdType actRegion = coilRegionId_[coil];
//
//      // Get current region name
//      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
//
//      // create new entity list
//      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
//      actSDList->SetRegion( actRegion );
//
//      LinearForm * curInt = NULL;
//      std::string factor = coilDef_[coil]->value_ + "/" +
//          lexical_cast<std::string>(coilDef_[coil]->windingCrossSection_);
//      PtrCoefFct coef;
//      // ===========
//      //  3D CASE
//      // ===========
//      if( dim_ == 3 ) {
//        StdVector<std::string> currDensity(3);
//        currDensity[0] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[0]);
//        currDensity[1] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[1]);
//        currDensity[2] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[2]);
//
//        if( isComplex_ ) {
//          StdVector<std::string> phaseVec(3);
//          phaseVec.Init(coilDef_[coil]->phase_);
//          coef = CoefFunction::Generate(mp_, Global::COMPLEX, currDensity, phaseVec );
//          coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
//          curInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,3,3>(),
//                                              1.0, coef, updatedGeo_);
//        } else {
//          coef = CoefFunction::Generate(mp_, Global::REAL, currDensity);
//          coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
//          curInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3>(),
//                                             1.0, coef, updatedGeo_);
//        } // complex
//      } else {
//        // ===============
//        //  2D / AXI CASE
//        // ===============
//        StdVector<std::string> currDensity(1);
//        currDensity[0] = factor;
//
//        if( isComplex_ ) {
//          StdVector<std::string> phaseVec(1);
//          phaseVec.Init(coilDef_[coil]->phase_);
//          coef = CoefFunction::Generate(mp_, Global::COMPLEX, currDensity, phaseVec );
//          //coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_;
//          curInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1>(),
//                                              1.0, coef, updatedGeo_);
//        } else {
//          coef = CoefFunction::Generate(mp_, Global::REAL, currDensity);
//          //coef->SetCoordinateSystem(coilDef_[coil]->flowCoordSys_);
//          curInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1>(),
//                                             1.0, coef, updatedGeo_);
//
//        } // complex
//      } // dimension
//      
//      // remember coefficient for later use
//      coilCoefs_[actRegion] = coef;
//
//      LinearFormContext * coilContext =
//          new LinearFormContext( curInt );
//      coilContext->SetEntities( actSDList );
//      coilContext->SetFeFunction( feFct );
//      assemble_->AddLinearForm( coilContext );
//
//    }
    
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magpde) << "Reading prescribed flux density";
    StdVector<std::string> vecComponents;
    if( dim_ == 3 ) {
      vecComponents = "x", "y", "z";
    }
    else if( isaxi_ ) {
      vecComponents = "r", "z";
    } 
    else {
      vecComponents = "x", "y";
    }
    
    LinearForm * lin = NULL;
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    bool coefUpdateGeo = true;
    ReadRhsExcitation( "fluxDensity", vecComponents, ResultInfo::VECTOR, isComplex_, 
                       ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be defined im volume")
      }

      

      if(isComplex_) {
        if( dim_ == 2) {
          if( isaxi_ ) {
            // axisymmetric case
            lin = new BDUIntegrator<CurlOperatorAxi<Double>, Complex>(Complex(1.0), 
                                                                     coef[i], reluc_, coefUpdateGeo);
          } else {
            lin = new BDUIntegrator<CurlOperator<FeH1,2,Double>,Complex>(Complex(1.0), 
                                                                        coef[i], reluc_,coefUpdateGeo);
          }
        } else {
          // 3D case
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Double>,Complex>(Complex(1.0), 
                                                                      coef[i], reluc_,coefUpdateGeo);
        }
      } else {
        if( dim_ == 2) {
          if( isaxi_ ) {
            // axisymmetric case
            lin = new BDUIntegrator<CurlOperatorAxi<Double>, Double>(1.0,  coef[i], reluc_,coefUpdateGeo);
          } else {
            lin = new BDUIntegrator<CurlOperator<FeH1,2,Double>,Double>(1.0, coef[i], reluc_,coefUpdateGeo);
          }
        } else {
          // 3D case
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Double>,Double>(1.0, coef[i], reluc_,coefUpdateGeo);
        }  
      }
      
      lin->SetName("FluxIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(feFct);
      assemble_->AddLinearForm(ctx);
      feFct->AddEntityList(ent[i]);
    } // for
  }
  
  void MagneticPDE::ReadSpecialBCs() {


    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();
  }


 
  void MagneticPDE::DefineSolveStep()
  {
		  solveStep_ = new StdSolveStep(*this);
  }




  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MagneticPDE::InitTimeStepping() {
    // Use complete implicit scheme
    Double gamma = 1;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    if( hasVoltCoils_ ){
      // Important: Create a new time scheme just for the current unknowns, as otherwise the
      // size of the vectors does not match!
      GLMScheme * scheme2 = new Trapezoidal(gamma);
      shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(scheme2, 0, nlType) );

      feFunctions_[COIL_CURRENT]->SetTimeScheme(myScheme2);
    }
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagneticPDE::ReadCoils() {

    // Check if the element "coils" is present at all.
    // Otherwise leave
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    PtrParamNode coilInfoNode = myInfo_->Get( "coilList", ParamNode::PASS );
    if ( !coilNode )
      return;

    // Get single coil nodes
    ParamNodeList coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    if( coilNodes.GetSize() > 0 ) {
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {

        // get coil and id
        std::string coilId = coilNodes[i]->Get("id")->As<std::string>();

        // Check if coil with same ID already exists
        if( coils_.find(coilId) != coils_.end() ) {
          EXCEPTION("A coil with ID '" << coilId << "' was already defined.")
        }

        // Create new coil
        shared_ptr<Coil> actCoil( new Coil( coilNodes[i], coilInfoNode, 
                                            ptGrid_, mp_, part) );
        coils_[coilId] = actCoil;

        // Associate mapping of coil parts with regions
        std::map<RegionIdType, shared_ptr<Coil::Part> >::const_iterator it;
        for( it = actCoil->parts_.begin(); it != actCoil->parts_.end(); it++ ) {
          coilRegions_[it->first] = actCoil;
        }
      }

      // Adjust printing of coil information to info node
      WARN("Adapt printing of coils to InfoNode");
    }
  }

  void MagneticPDE::DefinePrimaryResults()
  {
    StdVector<std::string> vecComponents;
    if( dim_ == 3 ) {
      vecComponents = "x", "y", "z";
    }
    else if( isaxi_ ) {
      vecComponents = "phi";
    } 
    else {
      vecComponents = "z";
    }
    
    shared_ptr<BaseFeFunction> vecFct = feFunctions_[MAG_POTENTIAL];
    
    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = vecComponents;
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::VECTOR;
    res1->unit = "Vs/m";
    res1->SetFeFunction(vecFct);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    vecFct->SetResultInfo(res1);
    DefineFieldResult( vecFct, res1);
    
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";

    
    // === MAGNETIC SCALAR POTENTIAL ===
    if (isMixed_) {
      shared_ptr<BaseFeFunction> scalFct = feFunctions_[ELEC_POTENTIAL];
      shared_ptr<ResultInfo> res2(new ResultInfo);
      res2->resultType = ELEC_POTENTIAL;
      res2->dofNames = "";
      res2->unit = "V";
      res2->definedOn = ResultInfo::NODE;
      res2->entryType = ResultInfo::SCALAR;
      results_.Push_back( res2 );
      availResults_.insert( res2 );
      scalFct->SetResultInfo(res2);
      DefineFieldResult( scalFct, res2 );
      
      // -----------------------------------
      //  Define xml-names of Dirichlet BCs
      // -----------------------------------
      hdbcSolNameMap_[ELEC_POTENTIAL] = "ground";
      idbcSolNameMap_[ELEC_POTENTIAL] = "elecPotential";
    }

    // === COIL CURRENT ===
    if( hasVoltCoils_ ){
      shared_ptr<ResultInfo> currentInfo(new ResultInfo);
      currentInfo->resultType = COIL_CURRENT;
      currentInfo->dofNames = "";
      currentInfo->unit = "A";
      currentInfo->definedOn = ResultInfo::COIL;
      currentInfo->entryType = ResultInfo::SCALAR;

      feFunctions_[COIL_CURRENT]->SetResultInfo(currentInfo);
      DefineFieldResult( feFunctions_[COIL_CURRENT], currentInfo );
    }

  }
    
  void MagneticPDE::DefinePostProcResults() {
    StdVector<std::string> vecComponents, aVecComponents;
      if( dim_ == 3 ) {
        vecComponents = "x", "y", "z";
        aVecComponents = "x", "y", "z";
      }
      else if( isaxi_ ) {
        vecComponents = "r", "z";
        aVecComponents = "phi";
      } 
      else {
        vecComponents = "x", "y";
        aVecComponents = "z";
      }
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

      // === MAGNETIC VECTOR POTENTIAL - 1ST DERIVATIVE ===
      if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC ) {
        shared_ptr<ResultInfo> aDot(new ResultInfo);
        aDot->resultType = MAG_POTENTIAL_DERIV1;
        aDot->dofNames = aVecComponents;
        aDot->unit = "V/m";
        aDot->definedOn = ResultInfo::ELEMENT;
        aDot->entryType = ResultInfo::VECTOR;
        availResults_.insert( aDot );
        DefineTimeDerivResult( MAG_POTENTIAL_DERIV1, 1, MAG_POTENTIAL );

        // === COIL CURRENT OF COILS EXCITED BY VOLTAGE - 1ST DERIVATIVE ===
        if( hasVoltCoils_ ){
          shared_ptr<ResultInfo> iDot(new ResultInfo);
          iDot->resultType = COIL_CURRENT_DERIV1;
          iDot->dofNames = "";
          iDot->unit = "A/s";
          iDot->definedOn = ResultInfo::COIL;
          iDot->entryType = ResultInfo::SCALAR;
          availResults_.insert( iDot );
          DefineTimeDerivResult( COIL_CURRENT_DERIV1, 1, COIL_CURRENT );
        }
      }

      // === MAGNETIC RHS ===
      shared_ptr<ResultInfo> rhs(new ResultInfo);
      rhs->resultType = MAG_RHS_LOAD;
      rhs->dofNames = aVecComponents;
      rhs->unit = "";
      rhs->entryType = ResultInfo::VECTOR;
      rhs->definedOn = ResultInfo::NODE;
      rhsFeFunctions_[MAG_POTENTIAL]->SetResultInfo(rhs);
      DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL], rhs );

      // === MAGNETIC FLUX DENSITY ===
      shared_ptr<ResultInfo> flux(new ResultInfo);
      flux->resultType = MAG_FLUX_DENSITY;
      flux->dofNames = vecComponents;
      flux->unit = "Vs/m^2";
      flux->definedOn = ResultInfo::ELEMENT;
      flux->entryType = ResultInfo::VECTOR;
      availResults_.insert( flux );
      shared_ptr<CoefFunctionFormBased> bFunc;
      if( isComplex_ ) {
        bFunc.reset(new CoefFunctionBOp<Complex>(feFct, flux));
      } else {
        bFunc.reset(new CoefFunctionBOp<Double>(feFct, flux));
      }
      DefineFieldResult( bFunc, flux );
      stiffFormCoefs_.insert(bFunc);

      // === EDDY CURRENT DENSITY ===
      shared_ptr<CoefFunctionFormBased> jFunc;
      if( analysistype_ != STATIC ) {
        shared_ptr<BaseFeFunction> aDotFct = 
            timeDerivFeFunctions_[MAG_POTENTIAL_DERIV1];
        shared_ptr<ResultInfo> eddy(new ResultInfo);
        eddy->resultType = MAG_EDDY_CURRENT_DENSITY;
        eddy->dofNames = "";
        eddy->unit = "A/m^2";
        eddy->definedOn = ResultInfo::ELEMENT;
        eddy->entryType = ResultInfo::VECTOR;
        availResults_.insert( eddy );

        if( isMixed_) 
          WARN("Adjust eddy currents for mixed case");
        if( isComplex_ ) {
          jFunc.reset(new CoefFunctionFlux<Complex>(aDotFct, eddy, -1.0));
        } else {
          jFunc.reset(new CoefFunctionFlux<Double>(aDotFct, eddy, -1.0));
        }
        DefineFieldResult( jFunc, eddy );
        massFormCoefs_.insert(jFunc);

        // === EDDY POWER DENSITY ===
        shared_ptr<ResultInfo> epd(new ResultInfo());
        epd->resultType = MAG_EDDY_POWER_DENSITY;
        epd->dofNames = "";
        epd->unit = "W/m^3";
        epd->definedOn = ResultInfo::ELEMENT;
        epd->entryType = ResultInfo::SCALAR;
        shared_ptr<CoefFunctionFormBased> epdFunctor;
        if( isMixed_) 
          WARN("Adjust eddy power density for mixed case");
        if( isComplex_ ) { 
          epdFunctor.reset( new CoefFunctionBdBKernel<Complex>(aDotFct, 1.0));
        } else {
          epdFunctor.reset( new CoefFunctionBdBKernel<Double>(aDotFct, 1.0));
        }
        DefineFieldResult( epdFunctor, epd );
        massFormCoefs_.insert(epdFunctor);

        // === EDDY POWER ===
        shared_ptr<ResultInfo> ep(new ResultInfo());
        ep->resultType = MAG_EDDY_POWER;
        ep->dofNames = "";
        ep->unit = "W";
        ep->definedOn = ResultInfo::REGION;
        ep->entryType = ResultInfo::SCALAR;
        availResults_.insert( ep );
        if( isMixed_) 
          WARN("Adjust eddy power for mixed case");
        shared_ptr<ResultFunctor> epFunctor;
        if( isComplex_ ) {
          epFunctor.reset(new EnergyResultFunctor<Complex>(aDotFct, ep, 1.0));
        } else {
          epFunctor.reset(new EnergyResultFunctor<Double>(aDotFct, ep, 1.0));
        }
        resultFunctors_[MAG_EDDY_POWER] = epFunctor;
        massFormFunctors_.insert(epFunctor);
      }

      // determine dimensionality of current density
      UInt jDim = aVecComponents.GetSize();
      
      // === COIL CURRENT DENSITY ===
      shared_ptr<ResultInfo> ccd(new ResultInfo);
      ccd->resultType = MAG_COIL_CURRENT_DENSITY;
      ccd->dofNames = aVecComponents;
      ccd->unit = "A/m^2";
      ccd->definedOn = ResultInfo::ELEMENT;
      ccd->entryType = ResultInfo::VECTOR;
      availResults_.insert( ccd );
      shared_ptr<CoefFunctionMulti> ccdCoef(new CoefFunctionMulti(CoefFunction::VECTOR, jDim, 1,isComplex_));
      DefineFieldResult( ccdCoef, ccd );


      // === TOTAL CURRENT DENSITY ===
      shared_ptr<ResultInfo> tcd(new ResultInfo);
      tcd->resultType = MAG_TOTAL_CURRENT_DENSITY;
      tcd->dofNames = aVecComponents;
      tcd->unit = "A/m^2";
      tcd->definedOn = ResultInfo::ELEMENT;
      tcd->entryType = ResultInfo::VECTOR;
      availResults_.insert( tcd );
      shared_ptr<CoefFunctionMulti> tcdCoef(new CoefFunctionMulti(CoefFunction::VECTOR,jDim,1, 
                                                                  isComplex_));
      DefineFieldResult( tcdCoef, tcd );

      
      // === LORENTZ FORCE DENSITY ===
      shared_ptr<ResultInfo> lfd(new ResultInfo);
      lfd->resultType = MAG_FORCE_LORENTZ_DENSITY;
      lfd->dofNames = vecComponents;
      lfd->unit = "N/m^3";
      lfd->definedOn = ResultInfo::ELEMENT;
      lfd->entryType = ResultInfo::VECTOR;
      availResults_.insert( lfd );

      // assemble coefficient function F_L = J X B
      
      // switch type of cross-product depending on dimensionality
      CoefXpr::OpType op = isaxi_ ? CoefXpr::OP_CROSS_AXI : CoefXpr::OP_CROSS;
      PtrCoefFct lfdFunc = CoefFunction::Generate( mp_, part, 
                                                   CoefXprBinOp(mp_,  tcdCoef, bFunc, op) );
      DefineFieldResult( lfdFunc, lfd);

      // === LORENTZ FORCE (TOTAL) ===
      shared_ptr<ResultInfo> lf(new ResultInfo);
      lf->resultType = MAG_FORCE_LORENTZ;
      lf->dofNames = vecComponents;
      lf->unit = "N";
      lf->definedOn = ResultInfo::REGION;
      lf->entryType = ResultInfo::VECTOR;
      availResults_.insert( lf );

      // build result functor for integration
      shared_ptr<ResultFunctor> lfFunc;
      if( isComplex_ ) {
        lfFunc.reset(new ResultFunctorIntegrate<Complex>(lfdFunc, feFct, lf ) );
      } else {
        lfFunc.reset(new ResultFunctorIntegrate<Double>(lfdFunc, feFct, lf ) );
      }
      resultFunctors_[MAG_FORCE_LORENTZ] = lfFunc;
      
      // === MAGNETIC ENERGY ===
      if( isMixed_)
        WARN("Adjust energy for mixed case.");
      shared_ptr<ResultInfo> energy(new ResultInfo);
      energy->resultType = MAG_ENERGY;
      energy->dofNames = "";
      energy->unit = "Ws";
      energy->definedOn = ResultInfo::REGION;
      energy->entryType = ResultInfo::SCALAR;
      availResults_.insert( energy );
      shared_ptr<ResultFunctor> energyFunc;
      if( isComplex_ ) {
        energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy, 0.5));
      } else {
        energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy, 0.5));
      }
      resultFunctors_[MAG_ENERGY] = energyFunc;
      stiffFormFunctors_.insert(energyFunc);

  }
  
  void MagneticPDE::FinalizePostProcResults() {

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

     //=== COIL CURRENT DENSITY ===
    shared_ptr<CoefFunctionMulti> ccdCoef
    = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_COIL_CURRENT_DENSITY]);
    // loop over all coil coefficients and add contribution to coef 
    std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCurrentDens_.begin();
    for( ; coilIt != coilCurrentDens_.end(); ++coilIt ) {
      ccdCoef->AddRegion( coilIt->first, coilIt->second);
    }

    // === TOTAL CURRENT DENSITY ===
    PtrCoefFct jEddy = GetCoefFct(MAG_EDDY_CURRENT_DENSITY);
    shared_ptr<CoefFunctionMulti> tcdCoef 
    = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
    // loop over all regions and assemble total current density:
    //  - if region is coil -> take coil current
    //  - if region is no coil and analyis is transient/harmonic -> eddy

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ) {
      RegionIdType actRegion = *regIt;
      if( coilCurrentDens_.find(actRegion) != coilCurrentDens_.end() ) {
        // region is a coil
        tcdCoef->AddRegion( actRegion, coilCurrentDens_[actRegion] );
      } else {
        // region is no coil
        if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC ) {
          tcdCoef->AddRegion( actRegion, jEddy );
        }
      }
    }
  }
  
  
  std::map<SolutionType, shared_ptr<FeSpace> >
   MagneticPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
     
    // and 3D analysis
    isMixed_ = false;
    if( analysistype_ != BasePDE::STATIC && dim_ == 3 ) {
      isMixed_ = true;
    }
    
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    
    if( formulation == "default" || formulation == "H1" ){
      
      // 1) create space for magnetic vector potential
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
      
      // 1) create space for electric scalar potential
      if( isMixed_ ) {
        crSpaces[ELEC_POTENTIAL] =
            FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
        crSpaces[ELEC_POTENTIAL]->Init(solStrat_);
      }

    }else{
       EXCEPTION( "The formulation " << formulation 
                  << "of the magnetic PDE is not known!" );
    }
    // ===================================================================
    // Create FeSpaceConst for coil currents (of coils excited by voltage)
    // ===================================================================
    if( hasVoltCoils_ ){
      PtrParamNode voltSpaceNode = infoNode->Get("coilCurrent");
      crSpaces[COIL_CURRENT] =
          FeSpace::CreateInstance(myParam_, voltSpaceNode, FeSpace::CONSTANT, ptGrid_);
      crSpaces[COIL_CURRENT]->Init(solStrat_);
    }
    return crSpaces;
   }

} // end namespace CoupledField
