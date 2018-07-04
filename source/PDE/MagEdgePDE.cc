#include <fstream>

#include "MagEdgePDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DECLARE_LOG(magEdgePde)
DEFINE_LOG(magEdgePde, "magEdgePde")


  // **************
  //  Constructor
  // **************
  MagEdgePDE::MagEdgePDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdge";
    pdematerialclass_ = ELECTROMAGNETIC;

    //! Always use updated Lagrangian formulation 
    updatedGeo_        = true; //true;

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));

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


  // *************
  //  Destructor
  // *************
  MagEdgePDE::~MagEdgePDE() {
  }

  shared_ptr<Coil> MagEdgePDE::GetCoilById(const Coil::IdType& id) {
    return coils_.at(id);
  }

  void MagEdgePDE::ReadSpecialBCs() {


    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();

  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagEdgePDE::InitNonLin() {

    SinglePDE::InitNonLin();
  }


  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgePDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // flag for updatedLagrange formulation
    //bool upLagrangeForm = true;

    
    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);
    
    //==============================================================
    //begin new implementation
    //==============================================================
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      feSpace->SetRegionApproximation(actRegion,polyId,integId);
      
      //get possible nonlinearities defined in this region
      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // pass entitylist ot fespace / fefunction
      feFunc->AddEntityList( actSDList );
      
      // Switch, if region is linear / nonlinear
      if ( nonLinTypes.GetSize() > 0 ) {
        
        // ================================
        //  Nonlinear Stiffness Integrator 
        // =================================
        
        // create stiffness integrator
        //BaseBOperator* bOp = new CurlOperator<FeHCurl,3, Double>();
        PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
        PtrCoefFct nuNl = 
            actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY, Global::REAL, 
                                          magFluxCoef  );
        //compute permeability
        PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        PtrCoefFct permeability = CoefFunction::Generate( mp_,  Global::REAL,
         	            CoefXprBinOp(mp_, constOne, nuNl, CoefXpr::OP_DIV ) );
        matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

        BaseBDBInt* stiff1 = NULL;
        stiff1 = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), nuNl, 1.0, updatedGeo_) ;
        stiff1->SetName("CurlCurlIntegrator-NL");

       BiLinFormContext * stiffContext =
           new BiLinFormContext(stiff1, STIFFNESS );
       stiffContext->SetEntities( actSDList, actSDList );
       stiffContext->SetFeFunctions( feFunc, feFunc );
       assemble_->AddBiLinearForm( stiffContext );
       // Important: Add bdb-integrator to global list, as we need them later
       // for calculation of postprocessing results
       bdbInts_[actRegion] = stiff1;
       // add also material to global, distributed reluctivity coefficient function
       reluc_->AddRegion(actRegion, nuNl);

       // ================================================
       //  Nonlinear Stiffness Integrator (only Newton )
       // ================================================
       // Note: currently we set the nonlinear method hard-coded to NEWTON for
       // testing purpose
       if( nonLinMethod_ == NEWTON ) {

         //BaseBOperator* bOp = new CurlOperator<FeHCurl,3, Double>();
         PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
         PtrCoefFct nuDeriv = actMat->GetTensorCoefFncNonLin( MAG_RELUCTIVITY_DERIV, FULL,  
                                                            Global::REAL, magFluxCoef );
         
         //create stiffness integrator
         BiLinearForm* stiff2 = NULL;
         stiff2 = new BDBInt<>(new CurlOperator<FeHCurl,3, Double>(), nuDeriv, 1.0, updatedGeo_) ;
         stiff2->SetName("CurlCurlIntegrator-NL-Newton");
         stiff2->SetNewtonBilinearForm();

         BiLinFormContext * stiffContext2 =
             new BiLinFormContext(stiff2, STIFFNESS );
         stiffContext2->SetEntities( actSDList, actSDList );
         stiffContext2->SetFeFunctions( feFunc, feFunc );
         assemble_->AddBiLinearForm( stiffContext2 );
       }
        
       /*if ( analysistype_ == STATIC ) {
         // not needed for StdSolveStep
         // =================================
         //  Nonlinear RHS-integrator
         // =================================
         LinearForm * rhsNlinForm = new KXIntegrator<Double>(stiff1, -1.0, feFunc );
         rhsNlinForm->SetName("RHSNonLinForm");
         
         // Set explicitly the solution dependency
         rhsNlinForm->SetSolDependent();
         LinearFormContext * rhsNlinContext =
             new LinearFormContext( rhsNlinForm );
         rhsNlinContext->SetEntities( actSDList );
         rhsNlinContext->SetFeFunction( feFunc );
         assemble_->AddLinearForm( rhsNlinContext );
       }*/
      } else {

       // ***************************************
       // LINEAR PART
       // ***************************************

       // ===============================
       //  Standard Stiffness Integrator
       // ===============================
        PtrCoefFct curCoef = 
            //actMat->GetTensorCoefFnc(MAG_RELUCTIVITY,FULL,Global::REAL );
            actMat->GetScalCoefFnc(MAG_RELUCTIVITY,Global::REAL );
        //compute permeability
        PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        PtrCoefFct permeability = CoefFunction::Generate( mp_,  Global::REAL,
        					      	            CoefXprBinOp(mp_, constOne, curCoef, CoefXpr::OP_DIV ) );
        matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

        BaseBDBInt* curlcurl;
        //curlcurl = new BDBInt< CurlOperator<FeHCurl,3, Double> >(curCoef,1.0) ;
        curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), curCoef,1.0, updatedGeo_) ;
        curlcurl->SetName("CurlCurlIntegrator");

       BiLinFormContext * stiffContext =
           new BiLinFormContext(curlcurl, STIFFNESS );
       stiffContext->SetEntities( actSDList, actSDList );
       stiffContext->SetFeFunctions( feFunc, feFunc );
       assemble_->AddBiLinearForm( stiffContext );
       
       // Important: Add bdb-integrator to global list, as we need them later
       // for calculation of postprocessing results
       bdbInts_[actRegion] = curlcurl;
       
       // add also material to global, distributed reluctivity coefficient function
       reluc_->AddRegion(actRegion, curCoef);

     } // END OF NONLIN/LIN PART


      // ============================================================
      // Mass Matrix
      // ============================================================
      if ( matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1 && nonLinTypes.GetSize() == 0) {
          // =================================================
          // Pure temperature-dependent nonlinear Mass Matrix
          // =================================================
    	  // purely temperature dependent (no further non linearity)
    	  PtrCoefFct coef;
          shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
          StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
          shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );

          //get coeff-Fnc to evaluate the temperature
          ReadMaterialDependency( "electricConductivity", dispDofNames, ResultInfo::SCALAR, isComplex_,
                                  ent, coef, updatedGeo_ );
          //coef-Fnc for electric conductivity
          PtrCoefFct condNL = actMat->GetScalCoefFncNonLin( MAG_CONDUCTIVITY, Global::REAL, coef);
          if ( analysistype_ == STATIC ) {
        	  EXCEPTION("No temperature-dependent conductivity implemented for static magnetic analysis");
          }
          conduc_->AddRegion(actRegion, condNL);
          BaseBDBInt *massInt;
          BiLinFormContext * massContext;
          massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                                        condNL, 1.0, updatedGeo_);
          massContext = new BiLinFormContext(massInt, DAMPING );
          massInt->SetName("MassIntegrator-NL");
		  massContext->SetEntities( actSDList, actSDList );
		  massContext->SetFeFunctions( feFunc, feFunc );
		  assemble_->AddBiLinearForm( massContext );

	      // insert mass integrator to list of defined mass integrators
	  	  massInts_[actRegion] = massInt;

      }else{
          // =================================================
          // Standard linear Mass Matrix
          // =================================================
    	  Double conductivity = 0.0; //
    	        materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
    	        bool scaleByEdgeSize = false;
    	        // use gradient of shape functions?
    	        //bool useGrad = true;
    	        if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
    	          Matrix<Double> reluc;
    	          // get tensor of permeability and determine max. value
    	          materials_[actRegion]->GetTensor( reluc, MAG_RELUCTIVITY, Global::REAL );
    	          conductivity =  regularizationFactor * reluc[0][0];
    	          scaleByEdgeSize = true;
    	          // add region to set of "regularized" regions
    	          regularizedRegions_.insert(actRegion);
    	        }

    	        PtrCoefFct coeff =
    	            CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
    	        // add also material to global, distributed reluctivity coefficient function
    	        conduc_->AddRegion(actRegion, coeff);
    	        BaseBDBInt *massInt;

    	        BiLinFormContext * massContext;
    	        if ( analysistype_ == STATIC) {
    	          // we have to guarantee, that we add some mass to curl-curl integrator.
    	          // Additionally, the integrator gets scaled by the edge size for a uniform
    	          // conditioning
    	          massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),
    	                                        coeff,1.0);
    	          massInt->SetName("MassIntegrator");
    	          massContext =  new BiLinFormContext(massInt, STIFFNESS );
    	        } else {
    	          // here we add the "normal" mass integrator, which gets not scaled by the
    	          // edge size
    	          if( scaleByEdgeSize ) {
    	            massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),
    	                                          coeff,1.0, updatedGeo_);
    	            massContext = new BiLinFormContext(massInt, STIFFNESS );
    	          } else {
    	            massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(),
    	                                          coeff,1.0, updatedGeo_);
    	            massContext = new BiLinFormContext(massInt, DAMPING );
    	          }
    	          massInt->SetName("MassIntegrator");

    	        }
    	        massContext->SetEntities( actSDList, actSDList );
    	        massContext->SetFeFunctions( feFunc, feFunc );
    	        assemble_->AddBiLinearForm( massContext );

    	        // insert mass integrator to list of defined mass integrators
    	        massInts_[actRegion] = massInt;
      }// End of nonlin/lin mass matrix part

      

    } // end for regions
    
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
      if(( actCoil.sourceType_ == Coil::CURRENT )||
         ( actCoil.sourceType_ == Coil::EXTERNAL )) {

        /*
        =====================================================
         1) CURRENT driven coils OR EXTERNAL current density

         Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
              Sens. and Act., 2nd edition, p. 131ff
        =====================================================
        */

        for( partIt = actCoil.parts_.begin();
            partIt != actCoil.parts_.end();
            partIt++ ) {
          Coil::Part & actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
          actSDList->SetRegion( actRegion );

          // generate source current vector
          PtrCoefFct jFct;
          if( actCoil.sourceType_ == Coil::CURRENT ){
            CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_,
                                                   CoefXpr::OP_MULT);
            PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

            CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect),
                                                   CoefXpr::OP_DIV);
            jFct = CoefFunction::Generate(mp_, part, jVec);
          } else {
            jFct = coilPartsExtJ_[partIt->second];
          }
          coilCurrentDens_[actRegion] = jFct;
          LinearForm* curInt;
          if( isComplex_ ) {
            curInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                                                1.0, jFct, updatedGeo_);
          }
          else {
            curInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                                               1.0, jFct, updatedGeo_);
          }
          curInt->SetName("CoilIntegrator");
          LinearFormContext * coilContext =
              new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );
        } // loop: parts

      } else {

        /*
        ============================================
         2) VOLTAGE driven coils

         Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
              Sens. and Act., 2nd edition, p. 211ff
        ============================================
         The coupled equation system in this case looks like
           ( M_A     0 ) ( A_dot ) + ( K_A -f_A ) ( A ) = ( 0 )
           ( (f_A)^T 0 ) ( i_dot )   ( 0     R  ) ( i )   ( u )
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

          // implementation of coil current density is difficult because of FeSpaceConst;
          // it looks simple: J = I/Gamma_c, where Gamma_c is the coil cross section;
          // 1) but the FeSpaceConst does not have elements and the CoefFunction asks
          //    for elements in order to evaluate its expression (FeFunction::GetScalar);
          //    although the origin of this problem does not seem to be the FeFunction;
          //    the coil result coil current works properly
          // 2) the automatic calculation of the cross section of the coil
          //    must be implemented because we need the number of turns or the coil cross
          //    section (only the winding cross section is not enough!)
          //    or
          //    number of turns or coil cross section could be additionally specified in xml;
          // However, the effects are not severe: Current density results are not available
          // for coils with voltage sources, which is generally not interesting anyway.
          // With the 2 points resolved, the code could look like:
          /*CoefXprVecScalOp testOp = CoefXprVecScalOp( mp_, actPart.jUnitVec,
            GetCoefFct( COIL_CURRENT ), CoefXpr::OP_MULT );
          PtrCoefFct test = CoefFunction::Generate( mp_, part, testOp );*/
          // now the division by the cross section would be necessary
          // the unit vector of the current density is added as dummy so that the other
          // current density results do not get mixed up
          coilCurrentDens_[actRegion] = actPart.jUnitVec;

          // === -f_A ===
          LinearForm* psiDotInt;
          if( isComplex_ ) {
            psiDotInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                -1.0, eJscaled, updatedGeo_);
          } else {
            psiDotInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                -1.0, eJscaled, updatedGeo_);
          }
          psiDotInt->SetName("CoilVoltCouplInt");

          bool assembleTransposed = false;
          BiLinearForm* pseudoBiLin = new BiLinWrappedLinForm( psiDotInt, assembleTransposed );
          BiLinFormContext* voltCoilContext = new BiLinFormContext( pseudoBiLin, STIFFNESS );
          voltCoilContext->SetEntities( actSDList, singleCoilList );
          voltCoilContext->SetFeFunctions( feFunc, feFunctions_[COIL_CURRENT] );
          voltCoilContext->SetCounterPart(false);
          assemble_->AddBiLinearForm( voltCoilContext );

          // === (f_A)^T ===
          if( analysistype_ != STATIC ){
            LinearForm* psiDotIntT;
            if( isComplex_ ) {
              psiDotIntT = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                  1.0, eJscaled, updatedGeo_);
            } else {
              psiDotIntT = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                  1.0, eJscaled, updatedGeo_);
            }
            psiDotIntT->SetName("CoilVoltCouplIntTransposed");

            assembleTransposed = true;
            BiLinearForm* pseudoBiLinT = new BiLinWrappedLinForm( psiDotIntT, assembleTransposed );
            BiLinFormContext* voltCoilContextT = new BiLinFormContext( pseudoBiLinT, DAMPING );
            voltCoilContextT->SetEntities( singleCoilList, actSDList );
            voltCoilContextT->SetFeFunctions( feFunctions_[COIL_CURRENT], feFunc );
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

      } // if: current / voltage driven
    } // loop: coils

    /*
      // ============================
      // PERMANENT MAGNETS
      // ============================
      //
      // check, if this subdomain is a permanent magnet
      for ( UInt perm = 0; perm < magnetsDomain_.GetSize(); perm++ ) {
        if ( actRegion == magnetsDomain_[perm] ) {
          EXCEPTION("Currently magnetic 3D with edge elements do not support permanent magnets");

          Vector<Double> magnetization(dim_);
          magnetization[0] = magnetsOriX_[perm];
          magnetization[1] = magnetsOriY_[perm];
          magnetization[2] = magnetsOriZ_[perm];

          // Get reluctivity for this domain and perform consistency check
          Double reluctivity;
          actMat->GetScalar(reluctivity,MAG_RELUCTIVITY,Global::REAL);

          std::string fncname = "none";
          LinearForm *permSource =
            new MagPerm3DInt(magnetization, reluctivity,
                             isaxi_, upLagrangeForm );

          LinearFormContext * permContext =
            new LinearFormContext( permSource );
          permContext->SetPtPde( this );
          permContext->SetResult( results_[0], actSDList );
          assemble_->AddLinearForm( permContext );
        }
      }*/

  } // end DefineIntegrators
  
  
  void MagEdgePDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
            endIt = ncInterfaces_.End();

    for ( ; ncIt != endIt; ++ncIt ) {
      if( analysistype_ == STATIC ){
        EXCEPTION("Nitsche interface not yet tested for static analysis!\n"
                  "You only have to delete this exception and verify the results");
      }
      switch (ncIt->type) {
        case NC_MORTAR:
          EXCEPTION("No Mortar nonconforming interface for magnetic PDE with edge elements.\n"
                    "Try using H1 nodal elements in MagneticPDE")
          break;
        case NC_NITSCHE:
          if (dim_ == 2)
            EXCEPTION("MagEdgePDE only works for 3D geometry!")
          else
            DefineNitscheCoupling<3,1>(MAG_POTENTIAL, *ncIt );
          break;
        default:
          EXCEPTION("Unknown type of ncInterface");
          break;
      }
    }
  }

  void MagEdgePDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(magEdgePde) << "Defining rhs load integrators for magEdgePDE";
       
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> vecDofNames = myFct->GetResultInfo()->dofNames;
    bool coefUpdateGeo = true;
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgePde) << "Reading prescribed flux density";

    ReadRhsExcitation( "fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, 
                       ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be specified in a volume!")
      }
      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
      PtrCoefFct factor = CoefFunction::Generate(mp_, part,
                                                 CoefXprBinOp(mp_, reluc_, coef[i] , CoefXpr::OP_MULT ) );

      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new CurlOperator<FeHCurl,3, Complex>(),
                                         Complex(1.0), factor, coefUpdateGeo);
      } else {
        lin = new BUIntegrator<Double>( new CurlOperator<FeHCurl,3, Double>(),
                                        1.0, factor, coefUpdateGeo);
      }
      lin->SetName("FluxIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];
    } // for
    
    // ==================
    //  FIELD INTENSITY
    // ==================
    LOG_DBG(magEdgePde) << "Reading prescribed field intensity";

    ReadRhsExcitation( "fieldIntensity", vecDofNames, ResultInfo::VECTOR, isComplex_, 
                       ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic field intensity can only be specified in a volume!")
      }

      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new CurlOperator<FeHCurl,3, Complex>(),
                                         Complex(1.0), coef[i], coefUpdateGeo);
      } else {
        lin = new BUIntegrator<Double>( new CurlOperator<FeHCurl,3, Double>(),
                                        1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("FieldIntensityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

    } // for
  }
  

  void MagEdgePDE::DefineSolveStep()
  {
    /*SolveStepMagEdge *magSolveStep = new SolveStepMagEdge(*this);
    solveStep_ = magSolveStep;*/
    
    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgePDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
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
  void MagEdgePDE::ReadCoils() {
    // Check if the element "coils" is present at all.
    // Otherwise leave
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    PtrParamNode coilInfoNode = myInfo_->Get( "coilList", ParamNode::PASS );
    if ( !coilNode )
      return;

    // Get single coil nodes
    ParamNodeList coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    Global::ComplexPart cplx = isComplex_ ? Global::COMPLEX : Global::REAL;
    if( coilNodes.GetSize() > 0 ) {
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {
        
        // get coil and id
        std::string coilId = coilNodes[i]->Get("id")->As<std::string>();

        // Check if coil with same ID already exists
        if( coils_.find(coilId) != coils_.end() ) {
          EXCEPTION("A coil with ID '" << coilId << "' was already defined.")
        }

        // Create new coil
        shared_ptr<Coil> actCoil( new Coil( coilNodes[i], coilInfoNode, 
                                            ptGrid_, mp_, cplx ) );
        coils_[coilId] = actCoil;

        // Associate mapping of coil parts with regions
        std::map<RegionIdType, shared_ptr<Coil::Part> >::const_iterator it;
        for( it = actCoil->parts_.begin(); it != actCoil->parts_.end(); it++ ) {
          coilRegions_[it->first] = actCoil;
        }
      }

      // Insert the current densities which are defined externally (simulation or sequence step).
      // This is done here because it is impossible for the coil to use a PDE pointer.
      // We have to distinguish between external current density direction and external source.
      // External source includes the direction, but not vice versa. Therefore, the external source
      // must be stored per part anyway, although it counts for the whole coil. Additionally, the
      // parts need the regions and coef functions.
      std::map<Coil::IdType, shared_ptr<Coil> >::iterator coilIt;
      for( coilIt = coils_.begin(); coilIt != coils_.end(); ++coilIt ){
        std::map<shared_ptr<Coil::Part>, PtrParamNode >::iterator extPartIt;
        for( extPartIt = coilIt->second->partsExtJDir_.begin();
            extPartIt != coilIt->second->partsExtJDir_.end(); ++extPartIt ){
          PtrParamNode extNode = extPartIt->second;
          // determine if normalise is set
          bool normalise = true;
          if ( extNode->Has("normalise") ) {
              if ( extNode->Get("normalise")->As<std::string>() == "no" ) {
                  normalise = false;
              }
          }

          shared_ptr<CoefFunctionMulti> unitCurrDens(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1,
              isComplex_));
          shared_ptr<CoefFunctionMulti> currDens(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1,
                        isComplex_));
          for( UInt k_reg = 0; k_reg < extPartIt->first->regions.GetSize(); ++k_reg ){
            std::string regName = ptGrid_->regionData[extPartIt->first->regions[k_reg]].name;
            shared_ptr<EntityList> elems;
            elems = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regName );
            PtrCoefFct regCurrDens; // ReadUserFieldValues assigns a value to this
            StdVector<std::string> vecComponents;
            vecComponents = "x", "y", "z";
            std::set<UInt> definedDofs; // ReadUserFieldValues assigns a value to this
            bool updateGeo; // ReadUserFieldValues assigns a value to this
            ReadUserFieldValues(elems,extNode,vecComponents,
                ResultInfo::VECTOR,isComplex_,regCurrDens,
                definedDofs,updateGeo);
            // take the read values and normalise to a length of 1
            PtrCoefFct unitDir;
            if ( normalise ) {
                CoefXprUnaryOp dirAbsOp = CoefXprUnaryOp( mp_, regCurrDens, CoefXpr::OP_NORM );
                PtrCoefFct dirAbs = CoefFunction::Generate( mp_, cplx, dirAbsOp );
                CoefXprVecScalOp unitOp = CoefXprVecScalOp( mp_, regCurrDens, dirAbs, CoefXpr::OP_DIV );
                unitDir = CoefFunction::Generate( mp_, cplx, unitOp );
            }
            else {
                unitDir = regCurrDens;
            }
            unitCurrDens->AddRegion(extPartIt->first->regions[k_reg],unitDir);
            if( coilIt->second->sourceType_ == Coil::EXTERNAL ){
              currDens->AddRegion(extPartIt->first->regions[k_reg],regCurrDens);
            }
          }
          extPartIt->first->jUnitVec = unitCurrDens;
          if( coilIt->second->sourceType_ == Coil::EXTERNAL ){
            coilPartsExtJ_[extPartIt->first] = currDens;
          }
        }
      }

      // Adjust printing of coil information to info node
      WARN("Adapt printing of coils to InfoNode");
    }
  }


  void MagEdgePDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;

    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL], potInfo );

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

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";

    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability ( new ResultInfo );
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult(permFct, permeability);
  }
  
  void MagEdgePDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // === TIME DERIVATIVES OF PRIMARY RESULTS ===
    if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC ) {
      // === MAGNETIC VECTOR POTENTIAL - 1ST DERIVATIVE ===
      shared_ptr<ResultInfo> aDot(new ResultInfo);
      aDot->resultType = MAG_POTENTIAL_DERIV1;
      aDot->dofNames = vecComponents;
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


    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
    fluxDens->resultType = MAG_FLUX_DENSITY;
    fluxDens->dofNames = vecComponents;
    fluxDens->unit = "Vs/m^2";
    fluxDens->definedOn = ResultInfo::ELEMENT;
    fluxDens->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> bFunc;
    if( isComplex_ ) {
      bFunc.reset(new CoefFunctionBOp<Complex>(feFct, fluxDens));
    } else {
      bFunc.reset(new CoefFunctionBOp<Double>(feFct, fluxDens));
    }
    DefineFieldResult( bFunc, fluxDens );
    stiffFormCoefs_.insert(bFunc);


    // === MAGNETIC NORMAL FLUX DENSITY ===
    shared_ptr<ResultInfo> normFlux(new ResultInfo);
    normFlux->resultType = MAG_NORMAL_FLUX_DENSITY;
    normFlux->dofNames = "";
    normFlux->unit = "Vs/m^2";
    normFlux->entryType = ResultInfo::SCALAR;
    normFlux->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionSurf> sNormFDens;
    sNormFDens.reset(new CoefFunctionSurf(true, 1.0, normFlux));
    DefineFieldResult( sNormFDens, normFlux );
    surfCoefFcts_[sNormFDens] = bFunc;


    // === MAGNETIC_FLUX ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX;
    flux->dofNames = "";
    flux->unit = "Vs";
    flux->entryType = ResultInfo::SCALAR;
    flux->definedOn = ResultInfo::SURF_REGION;
    shared_ptr<ResultFunctor> fluxFct;
    if( isComplex_ ) {
      fluxFct.reset(new ResultFunctorIntegrate<Complex>(sNormFDens,
                                                          feFct, flux ) );
    } else {
      fluxFct.reset(new ResultFunctorIntegrate<Double>(sNormFDens,
                                                         feFct, flux ) );
    }
    resultFunctors_[MAG_FLUX] = fluxFct;
    availResults_.insert(flux);


    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD;
    rhs->dofNames = vecComponents;
    rhs->unit = "-";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::ELEMENT;
    rhsFeFunctions_[MAG_POTENTIAL]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL], rhs );


    // === RESULTS RELATED TO TIME DERIVATIVES ===
    shared_ptr<CoefFunctionFormBased> jFunc;
    shared_ptr<CoefFunction> jPowerDensFunc;
    if( analysistype_ != STATIC ) {
      // === EDDY CURRENT DENSITY ===
      shared_ptr<BaseFeFunction> aDotFct = 
          timeDerivFeFunctions_[MAG_POTENTIAL_DERIV1];
      shared_ptr<ResultInfo> eddy(new ResultInfo);
      eddy->resultType = MAG_EDDY_CURRENT_DENSITY;
      eddy->dofNames = vecComponents;
      eddy->unit = "A/m^2";
      eddy->definedOn = ResultInfo::ELEMENT;
      eddy->entryType = ResultInfo::VECTOR;

      if( isComplex_ ) {
        jFunc.reset(new CoefFunctionFlux<Complex>(aDotFct, eddy, -1.0));
      } else {
        jFunc.reset(new CoefFunctionFlux<Double>(aDotFct, eddy, -1.0));
      }
      DefineFieldResult( jFunc, eddy );

      // === EDDY CURRENT (SURFACE RESULT) ===
      shared_ptr<ResultInfo> ec(new ResultInfo());
      ec->resultType = MAG_EDDY_CURRENT;
      ec->dofNames = "";
      ec->unit = "A";
      ec->definedOn = ResultInfo::SURF_REGION;
      ec->entryType = ResultInfo::SCALAR;
      availResults_.insert( ec );

      // first, create normal mapping
      shared_ptr<CoefFunctionSurf> ncd(new CoefFunctionSurf(true, 1.0, ec));
      surfCoefFcts_[ncd] = jFunc;

      // then, integrate values
      shared_ptr<ResultFunctor> eddyCurrentFunc;
      if( isComplex_ ) {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Complex>(ncd, 
                                                                  feFct, ec ) );
      } else {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Double>(ncd, 
                                                           feFct, ec ) );
      }
      resultFunctors_[MAG_EDDY_CURRENT] = eddyCurrentFunc;

      // === EDDY POWER DENSITY ===
      shared_ptr<ResultInfo> epd(new ResultInfo());
      epd->resultType = MAG_EDDY_POWER_DENSITY;
      epd->dofNames = "";
      epd->unit = "W/m^3";
      epd->definedOn = ResultInfo::ELEMENT;
      epd->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionFormBased> epdFunctor;
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
      shared_ptr<ResultFunctor> epFunctor;
      if( isComplex_ ) {
        epFunctor.reset(new EnergyResultFunctor<Complex>(aDotFct, ep, 1.0));
      } else {
        epFunctor.reset(new EnergyResultFunctor<Double>(aDotFct, ep, 1.0));
      }
      resultFunctors_[MAG_EDDY_POWER] = epFunctor;
      massFormFunctors_.insert(epFunctor);

      // === ELECTRIC FIELD INTENSITY ===
      shared_ptr<ResultInfo> elecIntens(new ResultInfo);
      elecIntens->resultType = ELEC_FIELD_INTENSITY;
      elecIntens->SetVectorDOFs(dim_, isaxi_);
      elecIntens->dofNames = vecComponents;
      elecIntens->unit = "V/m";
      elecIntens->definedOn = ResultInfo::ELEMENT;
      elecIntens->entryType = ResultInfo::VECTOR;

      // assemble coefficient function E = - dA/dt
      PtrCoefFct eIFuncTmp = CoefFunction::Generate( mp_, part,
                             CoefXprBinOp(mp_,  aDotFct, aDotFct, CoefXpr::OP_MULT ) );
      PtrCoefFct eIFunc = CoefFunction::Generate( mp_, part,
                             CoefXprBinOp(mp_,  "-1.0", aDotFct, CoefXpr::OP_MULT ) );
      DefineFieldResult( eIFunc, elecIntens);



      // === COIL LINKED FLUX - 1ST DERIVATIVE, CORRESPONDS TO INDUCED VOLTAGE ===
      shared_ptr<ResultInfo> psiDotRes(new ResultInfo());
      psiDotRes->resultType = COIL_INDUCED_VOLTAGE;
      psiDotRes->dofNames = "";
      psiDotRes->unit = "V";
      psiDotRes->definedOn = ResultInfo::COIL;
      psiDotRes->entryType = ResultInfo::SCALAR;

      availResults_.insert( psiDotRes );
      shared_ptr<ResultFunctor> psiDotFunc;
      shared_ptr<CoefFunctionMulti> psiDotDens(new CoefFunctionMulti(CoefFunction::SCALAR,1,1,
                                                                      isComplex_));
      if( isComplex_ ){
        psiDotFunc.reset( new ResultFunctorIntegrate<Complex>(psiDotDens, aDotFct, psiDotRes) );
      } else {
        psiDotFunc.reset( new ResultFunctorIntegrate<Double>(psiDotDens, aDotFct, psiDotRes) );
      }
      resultFunctors_[COIL_INDUCED_VOLTAGE] = psiDotFunc;
      // it is an integrated result but we need to save the coef function
      // somewhere for the finalization
      fieldCoefs_[COIL_INDUCED_VOLTAGE] = psiDotDens;

    }

    // === COIL CURRENT DENSITY ===
    shared_ptr<ResultInfo> ccd(new ResultInfo);
    ccd->resultType = MAG_COIL_CURRENT_DENSITY;
    ccd->dofNames = vecComponents;
    ccd->unit = "A/m^2";
    ccd->definedOn = ResultInfo::ELEMENT;
    ccd->entryType = ResultInfo::VECTOR;
    availResults_.insert( ccd );
    shared_ptr<CoefFunctionMulti> ccdCoef(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, 
                                                                isComplex_));
    DefineFieldResult( ccdCoef, ccd );


    // === TOTAL CURRENT DENSITY ===
    shared_ptr<ResultInfo> tcd(new ResultInfo);
    tcd->resultType = MAG_TOTAL_CURRENT_DENSITY;
    tcd->dofNames = vecComponents;
    tcd->unit = "A/m^2";
    tcd->definedOn = ResultInfo::ELEMENT;
    tcd->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> tcdCoef(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,
                                                                isComplex_));
    DefineFieldResult( tcdCoef, tcd );


    // === LORENTZ FORCE DENSITY ===
    shared_ptr<ResultInfo> lfd(new ResultInfo);
    lfd->resultType = MAG_FORCE_LORENTZ_DENSITY;
    lfd->dofNames = vecComponents;
    lfd->unit = "N/m^3";
    lfd->definedOn = ResultInfo::ELEMENT;
    lfd->entryType = ResultInfo::VECTOR;

    // assemble coefficient function F_L = J X B
    PtrCoefFct lfdFunc = CoefFunction::Generate( mp_, part, 
                                                 CoefXprBinOp(mp_,  tcdCoef, bFunc, CoefXpr::OP_CROSS ) );
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


    // === PERMEABILITY  ===
    shared_ptr<ResultInfo> perm(new ResultInfo);
    perm->resultType = MAG_ELEM_PERMEABILITY;
    perm->dofNames = "";
    perm->unit = "Vs/Am";
    perm->definedOn = ResultInfo::ELEMENT;
    perm->entryType = ResultInfo::SCALAR;

    PtrCoefFct muFunc = CoefFunction::Generate( mp_, part,
        CoefXprUnaryOp( mp_, reluc_, CoefXpr::OP_INV ) );
    DefineFieldResult( muFunc, perm );


    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> magIntens(new ResultInfo);
    magIntens->resultType = MAG_FIELD_INTENSITY;
    magIntens->SetVectorDOFs(dim_, isaxi_);
    magIntens->dofNames = vecComponents;
    magIntens->unit = "A/m";
    magIntens->definedOn = ResultInfo::ELEMENT;
    magIntens->entryType = ResultInfo::VECTOR;

    // assemble coefficient function hFunc = reluctivity * (flux density - remanence)
    // has to be done in finalize because RHS is not known now
    shared_ptr<CoefFunctionMulti> hFunc(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1,
                                                                isComplex_));
    DefineFieldResult( hFunc, magIntens );

    if( analysistype_ != HARMONIC ) {
    	// === MAXWELL FORCE DENSITY ===
    	shared_ptr<ResultInfo> mfd(new ResultInfo);
    	mfd->resultType = MAG_FORCE_MAXWELL_DENSITY;
    	mfd->dofNames = vecComponents;
    	mfd->unit = "N/m^3";
    	mfd->definedOn = ResultInfo::SURF_ELEM;
    	mfd->entryType = ResultInfo::VECTOR;
    	availResults_.insert( mfd );
    	// Note: The positive normal direction in this case is defined as the
    	//       inward facing one.
    	shared_ptr<CoefFunctionSurfMaxwell> maxForceDens(new CoefFunctionSurfMaxwell(false, matCoefs_, ptGrid_, -1.0, mfd));
    	DefineFieldResult( maxForceDens, mfd);
    	surfCoefFcts_[maxForceDens] = bFunc;

    	// === MAXWELL FORCE (TOTAL) ===
    	shared_ptr<ResultInfo> mf(new ResultInfo);
    	mf->resultType = MAG_FORCE_MAXWELL;
    	mf->dofNames = vecComponents;
    	mf->unit = "N";
    	mf->definedOn = ResultInfo::SURF_REGION;
    	mf->entryType = ResultInfo::VECTOR;
    	availResults_.insert( mf );

    	// build result functor for integration
    	shared_ptr<ResultFunctor> mfFunc;
    	if( isComplex_ ) {
    		mfFunc.reset(new ResultFunctorIntegrate<Complex>(maxForceDens, feFct, mf ) );
    	} else {
    		mfFunc.reset(new ResultFunctorIntegrate<Double>(maxForceDens, feFct, mf ) );
    	}
    	resultFunctors_[MAG_FORCE_MAXWELL] = mfFunc;

    	// === VIRTUAL WORK PRINCIPLE FORCE (TOTAL) ===
    	shared_ptr<ResultInfo> vwp(new ResultInfo);
    	vwp->resultType = MAG_FORCE_VWP;
    	vwp->dofNames = vecComponents;
    	vwp->unit = "N";
    	vwp->definedOn = ResultInfo::SURF_REGION;
    	vwp->entryType = ResultInfo::VECTOR;
    	availResults_.insert( vwp );

    	// define and save coefFunction
    	shared_ptr<CoefFunctionSurfVWP> vwpForce(new CoefFunctionSurfVWP(false, matCoefs_, 1.0, vwp));
    	surfCoefFcts_[vwpForce] = bFunc;

    	// build result functor for integration
    	shared_ptr<ResultFunctor> vwpFunc;
    	if( isComplex_ ) {
    	  vwpFunc.reset(new ResultFunctorVWP<Complex>(vwpForce, feFct, vwp, ptGrid_ ) );
    	} else {
    	  vwpFunc.reset(new ResultFunctorVWP<Double>(vwpForce, feFct, vwp, ptGrid_ ) );
    	}
    	resultFunctors_[MAG_FORCE_VWP] = vwpFunc;
    }


    // === MAGNETIC ENERGY ===
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


    // === COIL LINKED FLUX  ===
    shared_ptr<ResultInfo> psiRes(new ResultInfo());
    psiRes->resultType = COIL_LINKED_FLUX;
    psiRes->dofNames = "";
    psiRes->unit = "Vs/m^2";
    psiRes->definedOn = ResultInfo::COIL;
    psiRes->entryType = ResultInfo::SCALAR;

    availResults_.insert( psiRes );
    shared_ptr<ResultFunctor> psiFunc;
    shared_ptr<CoefFunctionMulti> psiDens(new CoefFunctionMulti(CoefFunction::SCALAR,1,1,
                                                                    isComplex_));
    if( isComplex_ ){
      psiFunc.reset( new ResultFunctorIntegrate<Complex>(psiDens, feFct, psiRes) );
    } else {
      psiFunc.reset( new ResultFunctorIntegrate<Double>(psiDens, feFct, psiRes) );
    }
    resultFunctors_[COIL_LINKED_FLUX] = psiFunc;
    // it is an integrated result but we need to save the coef function
    // somewhere for the finalization
    fieldCoefs_[COIL_LINKED_FLUX] = psiDens;


    // === CORE LOSS DENSITY ===
    shared_ptr<ResultInfo> cldRes(new ResultInfo());
    cldRes->resultType = MAG_CORE_LOSS_DENSITY;
    cldRes->dofNames = "";
    cldRes->unit = "W/kg";
    cldRes->definedOn = ResultInfo::ELEMENT;
    cldRes->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> coreLossDensCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    DefineFieldResult( coreLossDensCoef, cldRes );

    // === CORE LOSS ===
    shared_ptr<ResultInfo> clRes(new ResultInfo());
    clRes->resultType = MAG_CORE_LOSS;
    clRes->dofNames = "";
    clRes->unit = "W";
    clRes->definedOn = ResultInfo::REGION;
    clRes->entryType = ResultInfo::SCALAR;
    //DefineFieldResult( coreLossCoef, clRes );
    availResults_.insert( clRes );
    shared_ptr<ResultFunctor> coreLossFunc;
    shared_ptr<CoefFunctionMulti> coreLossCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    if( isComplex_ ){
      coreLossFunc.reset( new ResultFunctorIntegrate<Complex>(coreLossCoef, feFct, clRes) );
    } else {
      coreLossFunc.reset( new ResultFunctorIntegrate<Double>(coreLossCoef, feFct, clRes) );
    }
    resultFunctors_[MAG_CORE_LOSS] = coreLossFunc;
    // it is an integrated result but we need to save the coef function
    // somewhere for the finalization
    fieldCoefs_[MAG_CORE_LOSS] = coreLossCoef;


    // === JOULE LOSS Power DENSITY INTEGRATED OVER PERIOD	 ===
    if( analysistype_ != STATIC ){
      shared_ptr<ResultInfo> jld(new ResultInfo);
      jld->resultType = MAG_JOULE_LOSS_POWER_DENSITY;
      jld->dofNames = "";
      jld->unit = "W/m^3";
      jld->definedOn = ResultInfo::ELEMENT;
      jld->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionMulti> jldCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
      DefineFieldResult( jldCoef, jld );
    }
  }

  void MagEdgePDE::FinalizePostProcResults() {

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // Initialized results involving coils

    // === EDDY CURRENT DENSITY ===
    shared_ptr<CoefFunctionFormBased> jEddyCoef =
            dynamic_pointer_cast<CoefFunctionFormBased>(fieldCoefs_[MAG_EDDY_CURRENT_DENSITY]);
    std::map<RegionIdType, BaseBDBInt*>::iterator massIt = massInts_.begin();
    for( ; massIt != massInts_.end(); ++massIt ) {
      RegionIdType region = massIt->first;
      BaseBDBInt* massInt = massIt->second;
      
      // only assign region to jEddy, if it not a "regularized" region, i.e. 
      // only regions with "physical" conductivity get assigned
      if( regularizedRegions_.find(region) == regularizedRegions_.end()) {
        jEddyCoef->AddIntegrator( massInt, region);
      }
    }
    
    // === COIL CURRENT DENSITY ===
    shared_ptr<CoefFunctionMulti> ccdCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_COIL_CURRENT_DENSITY]);
    // loop over all coil coefficients and add contribution to coef 
    std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCurrentDens_.begin();
    for( ; coilIt != coilCurrentDens_.end(); ++coilIt ) {
      ccdCoef->AddRegion( coilIt->first, coilIt->second);
    }

    // === TOTAL CURRENT DENSITY ===
    PtrCoefFct jEddy = GetCoefFct(MAG_EDDY_CURRENT_DENSITY);
    shared_ptr<CoefFunctionMulti> tcdCoef = 
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
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

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // === H field ===
    // assemble coefficient function field intensity = reluctivity * (flux density - remanence)
    // the remanence is the RHS load flux density
    shared_ptr<CoefFunctionMulti> hCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
    regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ){
      if( bRHSRegions_.find( *regIt ) != bRHSRegions_.end() ){
        PtrCoefFct bPM = CoefFunction::Generate( mp_, part,
            CoefXprBinOp( mp_, GetCoefFct(MAG_FLUX_DENSITY),
                bRHSRegions_[*regIt], CoefXpr::OP_SUB ) );
        PtrCoefFct hPM = CoefFunction::Generate( mp_, part,
            CoefXprVecScalOp( mp_, bPM, reluc_, CoefXpr::OP_MULT ) );
        hCoef->AddRegion(*regIt,hPM);
      } else {
        PtrCoefFct h = CoefFunction::Generate( mp_, part,
            CoefXprBinOp( mp_, reluc_, GetCoefFct(MAG_FLUX_DENSITY), CoefXpr::OP_MULT ) );
        hCoef->AddRegion(*regIt,h);
      }
    }

    // === INDUCED VOLTAGE ===
    // This code assembles a CoefFunctionMulti containing all coil regions, which
    // is integrated by a ResultFunctorIntegrate. Ref: M. Kaltenbacher, Numer. Sim.
    // of. Mech. Sens. and Act., 2nd edition, p. 212, Eq. (7.19)
    // Every coil part has an own CoefFunction for the current density, which is
    // why the coil regions are searched to find the corresponding part.
    // It is assumed implicitly that a part cannot contain a region contained by
    // another part. If that happens, the CoefFunctionMulti throws.
    if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC ) {
      PtrCoefFct temp = GetCoefFct(COIL_INDUCED_VOLTAGE);
      shared_ptr<CoefFunctionMulti> psiDotDens =
          dynamic_pointer_cast<CoefFunctionMulti>(temp);
      CoilRegionMap::iterator coilRegsIt = coilRegions_.begin();
      for( ; coilRegsIt != coilRegions_.end(); ++coilRegsIt ){
        std::map<RegionIdType,shared_ptr<Coil::Part> >::iterator partIt =
            coilRegsIt->second->parts_.find( coilRegsIt->first );
        CoefXprVecScalOp eJscaledOp = CoefXprVecScalOp( mp_, partIt->second->jUnitVec,
            boost::lexical_cast<std::string>(partIt->second->wireCrossSect), CoefXpr::OP_DIV );
        PtrCoefFct eJscaled = CoefFunction::Generate( mp_, part, eJscaledOp );
        CoefXprBinOp integrandOp = CoefXprBinOp( mp_, eJscaled,
            GetCoefFct( MAG_POTENTIAL_DERIV1 ), CoefXpr::OP_MULT );
        PtrCoefFct integrand = CoefFunction::Generate( mp_, part, integrandOp );
        psiDotDens->AddRegion( coilRegsIt->first, integrand );
      }
    }

    // === COIL LINKED FLUX ===
    // same as for induced voltage, but with the vector potential instead
    // of the first time derivative of the vector potential
    PtrCoefFct temp = GetCoefFct(COIL_LINKED_FLUX);
    shared_ptr<CoefFunctionMulti> psiDotDens =
        dynamic_pointer_cast<CoefFunctionMulti>(temp);
    CoilRegionMap::iterator coilRegsIt = coilRegions_.begin();
    for( ; coilRegsIt != coilRegions_.end(); ++coilRegsIt ){
      std::map<RegionIdType,shared_ptr<Coil::Part> >::iterator partIt =
          coilRegsIt->second->parts_.find( coilRegsIt->first );
      CoefXprVecScalOp eJscaledOp = CoefXprVecScalOp( mp_, partIt->second->jUnitVec,
          boost::lexical_cast<std::string>(partIt->second->wireCrossSect), CoefXpr::OP_DIV );
      PtrCoefFct eJscaled = CoefFunction::Generate( mp_, part, eJscaledOp );
      CoefXprBinOp integrandOp = CoefXprBinOp( mp_, eJscaled,
          GetCoefFct( MAG_POTENTIAL ), CoefXpr::OP_MULT );
      PtrCoefFct integrand = CoefFunction::Generate( mp_, part, integrandOp );
      psiDotDens->AddRegion( coilRegsIt->first, integrand );
    }

    // === CORE LOSS DENSITY ===
    // This is a "density" per kg, not per m^3. That's why we cannot integrate it over
    // the volume to get the core loss (see below).

//TODO Core loss computation is disabled temporarily for harmonic analysis because otherwise we
// get an error in AddRegion, due to the multiplication of a real- and a complex-valued
// coefficient function. CHECK THIS !!!
    if( analysistype_ != HARMONIC ) {
		shared_ptr<CoefFunctionMulti> cldCoef =
			dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_CORE_LOSS_DENSITY]);
		BaseMaterial* actMat = NULL;
		regIt = regions_.Begin();
		for( ; regIt != regions_.End(); ++regIt ) {
		  RegionIdType actRegion = *regIt;
		  actMat = materials_[actRegion];
		  PtrCoefFct coreLossFct = actMat->GetScalCoefFncNonLin( CORE_LOSS, Global::REAL,
			  GetCoefFct(MAG_FLUX_DENSITY) );
		  cldCoef->AddRegion(actRegion, coreLossFct);
		}

		// === CORE LOSS ===
		// The core loss per kg has to be multiplied by the density to get it per m^3.
		// Then we can integrate over the volume in order to get the core loss.
		shared_ptr<CoefFunctionMulti> clCoef =
			dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_CORE_LOSS]);
		actMat = NULL;
		regIt = regions_.Begin();
		for( ; regIt != regions_.End(); ++regIt ) {
		  RegionIdType actRegion = *regIt;
		  actMat = materials_[actRegion];
		  // core loss density in W/kg
		  PtrCoefFct coreLossFct = actMat->GetScalCoefFncNonLin( CORE_LOSS, Global::REAL,
			  GetCoefFct(MAG_FLUX_DENSITY) );
		  // core loss density in W/m^3, which deserves to be called density
		  PtrCoefFct densFct = actMat->GetScalCoefFnc( DENSITY, Global::REAL );
		  PtrCoefFct coreLossDensCoef = CoefFunction::Generate( mp_, part,
			  CoefXprBinOp( mp_, coreLossFct, densFct, CoefXpr::OP_MULT ));
		  clCoef->AddRegion( actRegion, coreLossDensCoef );
		}
    }


    // === EDDY CURRENT (JOULE) LOSS DENSITY INTEGRATED===
    /* Already integrated over one period of excitation "signal"
     * Q = \gamma \omega^2 ||\hat{A}||^2 + \omega Im{ \hat{J}_i \hat{\overbar{A}} }
     * \hat{\overbar{A}} is the conjugate complex of \hat{A}
     */
    if( analysistype_ == HARMONIC){
		shared_ptr<CoefFunctionMulti> eddyLossCoef =
				dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY]);
		regIt = regions_.Begin();
		for( ; regIt != regions_.End(); ++regIt ) {
		  RegionIdType actRegion = *regIt;
		  // first part : \gamma \omega^2 ||\hat{A}||^2
		  PtrCoefFct normA = CoefFunction::Generate( mp_, part,
				  CoefXprUnaryOp( mp_, GetCoefFct(MAG_POTENTIAL), CoefXpr::OP_NORM) );
		  PtrCoefFct squareNormA = CoefFunction::Generate( mp_, part,
				  CoefXprBinOp( mp_, normA, normA, CoefXpr::OP_MULT ) );
		  PtrCoefFct tmp = CoefFunction::Generate( mp_, part,
				  CoefXprBinOp( mp_, squareNormA, "0.5*2*pi*f*2*pi*f", CoefXpr::OP_MULT ) );
		  PtrCoefFct part1 = CoefFunction::Generate( mp_, part,
		                      CoefXprBinOp( mp_, materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY,Global::REAL),
		                    		  tmp, CoefXpr::OP_MULT ) );

		  PtrCoefFct  part2 = NULL;
	      if( coilCurrentDens_.find(actRegion) != coilCurrentDens_.end() ) {
	          // region is a coil
			  // second part : \omega Im{ \hat{J}_i \hat{\overbar{A}} }
			  PtrCoefFct part2Tmp = CoefFunction::Generate( mp_, part,
					  CoefXprBinOp( mp_, GetCoefFct(MAG_POTENTIAL), GetCoefFct(MAG_COIL_CURRENT_DENSITY), CoefXpr::OP_MULT_CONJ ) );
			  PtrCoefFct part2TmpIM = CoefFunction::Generate( mp_, part,
					  CoefXprUnaryOp( mp_, part2Tmp, CoefXpr::OP_IM ) );
			  part2 = CoefFunction::Generate( mp_, part,
					  CoefXprBinOp( mp_, part2TmpIM, "pi*f", CoefXpr::OP_MULT ) );

	      }else{
	    	  part2 = CoefFunction::Generate( mp_, part, "0.0");
	      }

		  // add both parts
		  PtrCoefFct  partAdd =  CoefFunction::Generate( mp_, part,
				  CoefXprBinOp( mp_, part1, part2, CoefXpr::OP_ADD) );
		  // the division by the period length is already incorporated
  		  //PtrCoefFct  partAdd =  CoefFunction::Generate( mp_, part,
		  //				  CoefXprBinOp( mp_, partAddTmp, "1.0", CoefXpr::OP_DIV) );
		  eddyLossCoef->AddRegion(actRegion, partAdd);
		}
    }

    if( analysistype_ == TRANSIENT){
		shared_ptr<CoefFunctionMulti> eddyLossCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY]);
		regIt = regions_.Begin();
		for( ; regIt != regions_.End(); ++regIt ) {
		  RegionIdType actRegion = *regIt;
		  // J_total * E = J_total * dA/dt
		  PtrCoefFct partTmp = CoefFunction::Generate( mp_, part,
				  CoefXprBinOp( mp_, GetCoefFct(ELEC_FIELD_INTENSITY), GetCoefFct(MAG_TOTAL_CURRENT_DENSITY), CoefXpr::OP_MULT ) );
		  //PtrCoefFct partT = CoefFunction::Generate( mp_, part, CoefXprBinOp( mp_, "t", partTmp, CoefXpr::OP_MULT ) );

		  eddyLossCoef->AddRegion(actRegion, partTmp);
		}
    }

  }

  std::map<SolutionType, shared_ptr<FeSpace> > 
  MagEdgePDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H_CURL"){
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] = 
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_ );
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation 
                << "of magnetic edge PDE is not known!");
    }

    // in addition query, if special treatment of anisotropic elements
    // is activated
    if( myParam_->Has("thinElements") ) {
      Double aspectRatio = 0.0;
      aspectRatio = myParam_->Get("thinElements")->Get("maxAspectRatio")->As<Double>();
      dynamic_pointer_cast<FeSpaceHCurlHi>(crSpaces[MAG_POTENTIAL])
          ->TreatThinElements(aspectRatio);
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

} // end of namespace

