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
#include "Forms/Operators/SurfaceNormalFluxDensityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "Domain/CoefFunction/CoefFunctionHyst.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Domain/CoefFunction/CoefFunctionDiagTensorFromScalar.hh"
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"

#include "Driver/SolveSteps/SolveStepHyst.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

namespace CoupledField {

  DEFINE_LOG(magpde, "magpde")

  MagneticPDE::MagneticPDE(Grid * aptgrid, PtrParamNode paramNode,
          PtrParamNode infoNode,
          shared_ptr<SimState> simState, Domain* domain)
  :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magnetic";
    pdematerialclass_ = ELECTROMAGNETIC;
    formulation_ = MagBasePDE::NODAL;

    updatedGeo_        = !true; //! true needed for NC interface to give correct results while moving; in master false is used however; needs more tests

    // for magnetostrictive coupling
    isMagnetoStrictCoupled_ = false;
    mechanicPDE_ = NULL;

    anyRegionHasConductivity_ = false;

//    isMixed_ = false;
//    regionApproxSet_ = false;
//    anyRegionHasConductivity_ = false;
//    // can the reluctivity be complex? before the change it had the same type as the PDE
//    relucTensor_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_, dim_, false));
//    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
//
//    // determine if there are coils excited by voltage
//    hasVoltCoils_ = false;
//    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
//    if ( coilNode ){
//      ParamNodeList coilNodes = coilNode->GetChildren();
//      for( UInt k = 0; k < coilNodes.GetSize(); k++ ){
//        if( coilNodes[k]->Has("source") ){
//          std::string exType = coilNodes[k]->Get("source")->Get("type")->As<std::string>();
//          if( exType == "voltage" ){
//            hasVoltCoils_ = true;
//            break;
//          }
//        }
//      }
//    }

  }

  MagneticPDE::~MagneticPDE()
  {

  }

  void MagneticPDE::SetMagnetoStrictCoupling(SinglePDE *mechanicPDE)
  {
    mechanicPDE_ = mechanicPDE;
    isMagnetoStrictCoupled_ = true;

  }

//  shared_ptr<Coil> MagneticPDE::GetCoilById(const Coil::IdType& id) {
//    return coils_.at(id);
//  }

//  void MagneticPDE::InitNonLin() {
//    MagBasePDE::InitNonLin();
//  }

  void MagneticPDE::CheckForConductivity() {
    /*
     * Helper function that iterates over all regions and checks whether at least one
     * has a non-zero conductivity; if this is not the case, we do not have to consider
     * the electric potential in the 3d case
     */

    anyRegionHasConductivity_ = false;
    RegionIdType actRegion;
    PtrCoefFct conducCoef;
	  BaseMaterial * actMat = NULL;

    //  Loop over all regions
	  std::map<RegionIdType, BaseMaterial*>::iterator it;
	  //hysteresisCoefs_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));

	  for ( it = materials_.begin(); it != materials_.end(); it++ ) {

		  // Set current region and material
		  actRegion = it->first;
		  actMat = it->second;

		  conducCoef = actMat->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL);

      // use this change to directly store the conductivity
      conduc_->AddRegion(actRegion, conducCoef);

      if(!conducCoef->IsZero()){
//        std::cout << "Region with id=" << actRegion << " has non-zero conductivity" << std::endl;
        anyRegionHasConductivity_ = true;
      }
    }

//    if(!anyRegionHasConductivity_){
//      std::cout << "No region with non-zero conductivity found" << std::endl;
//    }
//
  }


  void MagneticPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

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

    bool isPermFrozen = false;

	  shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
	  shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();


    double factor = 1.0;
    if ( isMagnetoStrictCoupled_ == true ){
      // similar to the piezoelectric case we have to multiply the magnetic pde in the magnetostrictive case with -1 to
      // get a symmetric equation system (in mechanics we have -div(sigma) in magnetics +rot(H) -> multiply magnetics with -1)
      factor = -1.0;
    }

    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    //hysteresisCoefs_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];

		  // Get current region name
		  std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

		  shared_ptr<ElemList> actSDList;
      if(!regionApproxSet_){
        // create new entity list
        shared_ptr<ElemList> newSDList( new ElemList(ptGrid_ ) );
        actSDList = newSDList;
        actSDList->SetRegion( actRegion );
        SDLists_[actRegion] = actSDList;
      } else {
        actSDList = SDLists_[actRegion];
      }

		  // --------------------------
		  //  Set region approximation
		  // --------------------------
		  // --- Set the approximation for the current region ---
		  PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
		  std::string polyId = curRegNode->Get("polyId")->As<std::string>();
		  std::string integId = curRegNode->Get("integId")->As<std::string>();

      if(!regionApproxSet_){
        // --- Set the approximation for the current region ---
        mySpace->SetRegionApproximation(actRegion, polyId, integId);
      }

      myFct->AddEntityList( actSDList );

		  //get possible nonlinearities defined in this region
		  const StdVector<NonLinType>& nonLinTypes = regionNonLinTypes_[actRegion];

		  // ====================================================================
		  //  NONLINEAR BH RELATION (NON-HYSTERETIC)
		  // ====================================================================
		  if (  nonLinTypes.Find(PERMEABILITY) != -1 ) {
		    CoefFunctionOpt* cfo = NULL; // we might do optimization and then we have such a thing
			  PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
			  PtrCoefFct nuNl = actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY_SCALAR, Global::REAL, magFluxCoef);

        if(domain->HasDesign())
        {
          cfo = new CoefFunctionOpt(domain->GetDesign(), nuNl, MAG_RELUCTIVITY_SCALAR,  this);
          nuNl.reset(cfo);
        }
        
			  PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
			  PtrCoefFct permeability = CoefFunction::Generate( mp_,  Global::REAL,
                CoefXprBinOp(mp_, constOne, nuNl, CoefXpr::OP_DIV ) );
			  matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

			  BaseBDBInt * stiffInt = NULL;
			  if( dim_ == 2) {
				  if( isaxi_ ) {
					  // axisymmetric case
					  stiffInt = new BBInt<>(new CurlOperatorAxi<Double>(), nuNl,factor, updatedGeo_);
				  } else {
					  // plane 2D case
					  stiffInt = new BBInt<>(new CurlOperator<FeH1,2,Double>(), nuNl, factor, updatedGeo_);
				  }
			  } else {
				  // 3D case
				  stiffInt = new BBInt<>(new CurlOperator<FeH1,3,Double>(), nuNl, factor, updatedGeo_);
			  }
			  stiffInt->SetName("CurlCurlIntegrator-NL");
			  BiLinFormContext * stiffContext =
                new BiLinFormContext(stiffInt, STIFFNESS );
			  stiffContext->SetEntities( actSDList, actSDList );
			  stiffContext->SetFeFunctions( myFct, myFct );
			  assemble_->AddBiLinearForm( stiffContext );

        // when we have a CoefFunctionOpt, we tell it the proper form, which we only have now
        if(cfo)
          cfo->SetForm(stiffInt);

			  // Important: Add bdb-integrator to global list, as we need them later
			  // for calculation of postprocessing results
			  bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );
			  // add also material to global, distributed reluctivity coefficient function
			  //relucTensor_->AddRegion(actRegion, nuNl);

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
						  stiff2 = new BDBInt<>(new CurlOperatorAxi<Double>(), nuDeriv, factor, updatedGeo_);
					  } else {
						  // plane 2D case
						  stiff2 = new BDBInt<>(new CurlOperator<FeH1,2,Double>(), nuDeriv, factor, updatedGeo_);
					  }
				  } else {
					  // 3D case
					  stiff2 = new BDBInt<>(new CurlOperator<FeH1,3,Double>(), nuDeriv, factor, updatedGeo_);

				  }

				  stiff2->SetName("CurlCurlIntegrator-NL-Newton");
				  //! mark the bi-linear form to be a Newton part
				  stiff2->SetNewtonBiLinearForm();

				  BiLinFormContext * stiffContext2 =
                  new BiLinFormContext(stiff2, STIFFNESS );
				  stiffContext2->SetEntities( actSDList, actSDList );
				  stiffContext2->SetFeFunctions( myFct, myFct );
				  assemble_->AddBiLinearForm( stiffContext2 );
			  }
		  } else {
			  // ====================================================================
			  //  HYSTERESIS + LINEAR CASE
			  // ====================================================================
			  shared_ptr<CoefFunction > curCoef;
        
			  CoefFunctionOpt* cfo = NULL; // we might do optimization and then we have such a thing
        BaseBDBInt * stiffInt = NULL;

        if ( nonLinTypes.Find(HYSTERESIS) != -1 ){
          /* for both the delta material method as well as the std fixpoint method we have to know
           * which regions are affected by hysteresis
           */
          // NEW: coefFncHyst should already be created during DefinePostProcResults!
          curCoef = GenerateHystereticCoefFunctions(actRegion);

//          PtrCoefFct hystPol = hysteresisCoefs_->GetRegionCoef(actRegion);
//          curCoef = hystPol->GenerateMatCoefFnc("Reluctivity");
//
//          PtrCoefFct hystOutput = hystPol->GenerateOutputCoefFnc("MagPolarization");
//          polarization_->AddRegion( actRegion, hystOutput);
//
//          PtrCoefFct hystOutput2 = hystPol->GenerateOutputCoefFnc("MagMagnetization");
//          magnetization_->AddRegion( actRegion, hystOutput2);
//
//          PtrCoefFct hystOutput3 = hystPol->GenerateOutputCoefFnc("MagFieldIntensityHyst");
//          fieldIntensity_->AddRegion( actRegion, hystOutput3);

          stiffInt = GeHystStiffInt( factor, curCoef );
			  } else {
				  // ====================================================================
				  //  Standard Linear CASE (2D AND 3D)
				  // ====================================================================

          //check for frozen permeability
          isPermFrozen = false;
          const StdVector<NonLinType>& matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
          if ( matDepenTypes.Find(PERMEABILITY_FROZEN) != -1 ) {
            // Here, the magnetic permeability is not used from the material file but rather from a previous simulation.
            // Hence, we freeze the permeability from this last simulation step and use it for further computations, much
            // like freezing the permeability for an operating point.

            // we read the computed permeability from file
            shared_ptr<ResultInfo> resultInfo = GetResultInfo(MAG_ELEM_PERMEABILITY); 
            std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
            shared_ptr<EntityList> entity = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );

            PtrCoefFct permeability;
            // get coeff-Fnc for the magnetic permeability
            ReadMaterialDependency( "permeabilityFrozen", resultInfo->dofNames, resultInfo->entryType, false,
                                    entity, permeability, updatedGeo_ );
            // compute the reluctivity
            PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
            curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, permeability,
                                             CoefXpr::OP_DIV));      
            // For postprocessing                        
            matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);    
            isPermFrozen = true;                    
          } 
          else {
				    curCoef = actMat->GetTensorCoefFnc( MAG_RELUCTIVITY_TENSOR, tensorType, Global::REAL );

				    // for postprocessing
            if (materials_[actRegion]->GetSymmetryType(MAG_PERMEABILITY_TENSOR) == BaseMaterial::ISOTROPIC) {
              PtrCoefFct permeability = materials_[actRegion]->GetScalCoefFnc( MAG_PERMEABILITY_SCALAR, Global::REAL);
              matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);
            }
          }

			    if(domain->HasDesign())
			    {
			      cfo = new CoefFunctionOpt(domain->GetDesign(), curCoef, MAG_RELUCTIVITY_TENSOR, this);
			      curCoef.reset(cfo);
			    }

          if( dim_ == 2) {
            if( isaxi_ ) {
              // axisymmetric case
              if ( isPermFrozen ) {
                stiffInt = new BBInt<>(new CurlOperatorAxi<Double>(), curCoef, factor, updatedGeo_);
              }
              else {
              stiffInt = new BDBInt<>(new CurlOperatorAxi<Double>(), curCoef, factor, updatedGeo_);
              }
            } else {
              // plane 2D case
              if ( isPermFrozen ) {
                stiffInt = new BBInt<>(new CurlOperator<FeH1,2,Double>(), curCoef,factor, updatedGeo_);
              }
              else {
                stiffInt = new BDBInt<>(new CurlOperator<FeH1,2,Double>(), curCoef,factor, updatedGeo_);
              }
            }
          } else {
            // 3D case
            if ( isPermFrozen ) {
              stiffInt = new BBInt<>(new CurlOperator<FeH1,3,Double>(), curCoef, factor, updatedGeo_);
            }
            else {
              stiffInt = new BDBInt<>(new CurlOperator<FeH1,3,Double>(), curCoef, factor, updatedGeo_);
            }
          }
          
          // shifted to DefineRHSLoadIntegrator!
//          // check for fixed magnetization > defined in mat.xml
//          int hasFixedMagnetization = 0;
//          materials_[actRegion]->GetScalar(hasFixedMagnetization,PRESCRIBED_MAGNETIZATION);
//
//          // mRHSRegions_ gets filled in MagBasePDE.cc > InitMagnetization
//          PtrCoefFct magFncSharedPtr = mRHSRegions_[actRegion];
//            
//          if(hasFixedMagnetization == 1){
//            // add RHS integrator (compare rhsIntegrator flux)
//            LinearForm * lin = NULL;
//            bool coefUpdateGeo = true;
//            
//            if( dim_ == 2) {
//              if( isaxi_ ) {
//                // axisymmetric case
//                lin = new BUIntegrator<Double>(new CurlOperatorAxi<Double>(), factor,  magFncSharedPtr,coefUpdateGeo);
//              } else {
//                lin = new BUIntegrator<Double>(new CurlOperator<FeH1,2,Double>(),factor, magFncSharedPtr,coefUpdateGeo);
//              }
//            } else {
//              // 3D case
//              lin = new BUIntegrator<Double>(new CurlOperator<FeH1,3,Double>(),factor, magFncSharedPtr,coefUpdateGeo);
//            }
//      
//            lin->SetName("FixedMagnetizationIntegrator");
//            LinearFormContext *ctx = new LinearFormContext( lin );
//            ctx->SetEntities( actSDList );
//            ctx->SetFeFunction(myFct);
//            assemble_->AddLinearForm(ctx);
//            //feFct->AddEntityList(ent[i]);
//          }        
			  }

			  stiffInt->SetName("CurlCurlIntegrator");
			  stiffInt->SetFeSpace( mySpace);
			  BiLinFormContext* stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );
			  stiffIntDescr->SetEntities( actSDList, actSDList );
			  stiffIntDescr->SetFeFunctions( myFct, myFct );

			  assemble_->AddBiLinearForm( stiffIntDescr );
        
			  // when we have a CoefFunctionOpt, we tell it the proper form, which we only have now
        if(cfo)
          cfo->SetForm(stiffInt);

			  // Important: Add bdb-integrator to global list, as we need them later
			  // for calculation of postprocessing results
			  bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

        // NEW: curCoef is a specialized coef function responsible
        // for delivering the mat tensor > type is tensor
        //			  // add also material to global, distributed reluctivity coefficient function
        //			  if ( nonLinTypes.Find(HYSTERESIS) != -1){// || nonLinTypes.Find(HYSTERESIS_FIXPOINT) != -1 ){
        //			    /*
        //			     * we cannot directly add coefFunctionHyst to relucTensor_ as relucTensor_ expects tensorial coefFncs
        //			     * but coefFunctionHyst has to be a vector coefFnc
        //			     */
        //			  } else {
			  if ( isPermFrozen == false ) 
          relucTensor_->AddRegion(actRegion, curCoef);
        //			  }

			  // ====================================================================
			  //  3D CASE (additional stiffness integrator)
			  // ====================================================================
			  if( dim_ == 3 ) {
//          /*
//           * This term is required to satisfy the coloumb gauge in 3d
//           * According to "Finite Element Method in Magnetics" (Kuczmann, Ivanyi) this
//           * term gets assembled with nu_0, i.e. in the hysteresis case, we have to return only nu_0 (no deltaMAT!)
//           * > create another coeffunction for reluctivity here
//           */
//          if ( nonLinTypes.Find(HYSTERESIS) != -1 ){
//            WARN("Hysteresis in 3d not tested yet");
//            PtrCoefFct hystPol = hysteresisCoefs_->GetRegionCoef(actRegion);
//            curCoef = hystPol->GenerateMatCoefFnc("LinearReluctivity");
//          }

          std::string nu0String = "795774.7155";
          PtrCoefFct nu0Scal = CoefFunction::Generate(mp_, Global::REAL, nu0String);
          StdVector<PtrCoefFct> diagValues = StdVector<PtrCoefFct>(dim_);
          for(UInt i = 0; i < dim_; i++){
            diagValues[i] = nu0Scal;
          }

          shared_ptr<CoefFunction> nu0 (new CoefFunctionDiagTensorFromScalar(diagValues));

          BaseBDBInt * divInt = new BDBInt<>(new DivOperator<FeH1,3,Double>(), nu0, factor, updatedGeo_);
          // NACS version:
//				  BaseBDBInt * divInt = new BDBInt<>(new DivOperator<FeH1,3,Double>(), curCoef, factor, updatedGeo_);
				  divInt->SetFeSpace( mySpace );
				  divInt->SetName("DivDivIntegrator");
				  BiLinFormContext * divIntDescr =
                  new BiLinFormContext(divInt, STIFFNESS );
				  divIntDescr->SetEntities( actSDList, actSDList );
				  divIntDescr->SetFeFunctions( myFct, myFct );
				  assemble_->AddBiLinearForm( divIntDescr );
			  }

		  }  //end: nonlinear, hysteresis, linear

		  // ====================================================================
		  //  MASS - MATRIX (all dimensions)
		  // ====================================================================
		  // Note: here we have currently the problem, that we can only treat
		  // constant / double valued conductivity values, as we must check for
		  // non-zero. However, how do we do this in case of non-constant parameters?
		  Double conductivity;

		  materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
		  PtrCoefFct conducCoef =
              materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL);
		  //                                         lexical_cast<std::string>(conductivity));
		  // FROM NACS
//      conduc_->AddRegion(actRegion, conducCoef);

//      if (!conducCoef->IsZero()) { // does not work if coupled to mechanics > test case failing due to EddyPower
      // > possible reason: EddyPower is calculated from/via mass integrators; if conductivity is zero, no mass integrator
      //    is defined on this region and thus the EddyPower will not be computed (see SinglePDE::FinalizePostProcResults);
      //    in the coupled mech-mag case, a mass integrator is defined from mechanics; thus CFS might try to compute the
      //    EddyPower for regions without conductivity and then fail?
      {
			  BaseBDBInt *massInt = NULL;
			  if( dim_ == 2 ) {
				  massInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), conducCoef,factor, updatedGeo_);
			  } else {
				  massInt = new BBInt<>(new IdentityOperator<FeH1,3,3>(), conducCoef,factor, updatedGeo_ );
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
      //  3D MIXED CASE (coupling between magnetic vector and electric scalar potential) (A-V)-Form
      // ====================================================================
      /*
       *
       *
       */
      // FROM NACS
      if( isMixed_ && !conducCoef->IsZero() ) {
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
          BiLinearForm * phiDivInt = new BBInt<>(new GradientOperator<FeH1,3,1,Double>(), conducCoef, factor, updatedGeo_);
          phiDivInt->SetName("MassIntegrator_PhiPhi");
          BiLinFormContext * massContext = new BiLinFormContext(phiDivInt, DAMPING);
          massContext->SetEntities( actSDList, actSDList );
          massContext->SetFeFunctions( potFct, potFct );
          assemble_->AddBiLinearForm( massContext );
        }

        // ------------------------------
        // off-diagonal mass matrix M_A_phi
        // ------------------------------
        {
          //BiLinearForm * cplInt =
          BaseBDBInt * cplInt = new ABInt<>(new IdentityOperator<FeH1,3,3,Double>() , new GradientOperator<FeH1,3,1,Double>(), conducCoef, factor, updatedGeo_);
          cplInt->SetName("MassIntegrator_Coupling_Phi_A");
          BiLinFormContext * cplContext = new BiLinFormContext(cplInt, DAMPING );
          cplContext->SetCounterPart(true);
          cplContext->SetEntities( actSDList, actSDList );
          cplContext->SetFeFunctions( myFct, potFct );
          assemble_->AddBiLinearForm( cplContext );

          // store bilinear form for later use in current density calculation
          // FROM NACS
          std::set<shared_ptr<CoefFunctionFormBased> >::iterator formIt;
          for( formIt = mixedFormFunctor_.begin();
                  formIt != mixedFormFunctor_.end(); ++formIt ) {
            (*formIt)->AddIntegrator(cplInt, actRegion);
          }

        }
      } // mixed

      std::string velocityId = curRegNode->Get("velocityId")->As<std::string>();
      if(velocityId != "")
      {
        bool coefUpdateGeo;
        ReadRegionVelocityField(velocityId,actSDList,actRegion,coefUpdateGeo);

        //coef-Fnc for electric conductivity
        Matrix<Double> reluc;
        Double conductivity = 0.0;

        // get conductivity
        materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
        assert(conductivity != 0.0);
        //PtrCoefFct coeff = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
        PtrCoefFct coeff = materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL);

        // Create the integrators
        BaseBDBInt   *velocityStiff = NULL;

        if( isComplex_ )
        {
//          if(dim_ == 2)
//          {
//            velocityStiff  = new ABInt<>(new IdentityOperator<FeH1,2,1>(),new CurlOperatorMag<FeH1,2,Double,Complex>(), coeff, 1.0, coefUpdateGeo);
//          }
//          else
//          {
//            velocityStiff  = new ABInt<>(new IdentityOperator<FeH1,3,1>(),new CurlOperatorMag<FeH1,3,Double,Complex>(), coeff, 1.0, coefUpdateGeo);
//          }
          EXCEPTION("VelocityStiffInt is not implemented for complex case!");
        }
        else
        {
          if(dim_ == 2)
          {
            velocityStiff  = new ABInt<>(new IdentityOperator<FeH1,2,1>(),new CurlOperatorMag<FeH1,2,Double>(),coeff, -1.0, coefUpdateGeo);
          }
          else
          {
            velocityStiff  = new ABInt<>(new IdentityOperator<FeH1,3,1>(),new CurlOperatorMag<FeH1,3,Double>(),coeff, -1.0, coefUpdateGeo);
          }
        }
        assert(velocityStiff != NULL);

        velocityStiff->SetBCoefFunctionOpB(VelocityCoef_);
        velocityStiff->SetName("VelocityStiffInt");
        velocityInts_[actRegion] = velocityStiff;

        BiLinFormContext *VelocityContextStiff =  new BiLinFormContext(velocityStiff, STIFFNESS );
        VelocityContextStiff->SetEntities( actSDList, actSDList );
        VelocityContextStiff->SetFeFunctions( feFunctions_[MAG_POTENTIAL],myFct);
        assemble_->AddBiLinearForm( VelocityContextStiff );
      }

	  } // regions

	  // ============================
	  // COIL INTEGRATORS
	  // ============================
    DefineCoilIntegrators(factor);
  }


  void MagneticPDE::DefineSurfaceIntegrators() {
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    this->ptGrid_->MapEdges();
    this->ptGrid_->MapFaces();
    if (!bcNode)
      return;
  
    ParamNodeList blochNodesList = bcNode->GetList("blochPeriodic");
    for (UInt i = 0; i < blochNodesList.GetSize(); i++) {
      std::string str_value = blochNodesList[i]->Get("factor_value")->As<std::string>();
      std::string formulation = blochNodesList[i]->Get("formulation")->As<std::string>();
      
      PtrCoefFct factor = CoefFunction::Generate(mp_, Global::REAL, str_value);
      PtrCoefFct factorSqr = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, factor, 
                                                    factor, CoefXpr::OP_MULT));      
      PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");

      ParamNodeList regionsList = blochNodesList[i]->GetList("region");
      for (UInt j = 0; j < regionsList.GetSize(); j++) {
        std::string ncRegionName = regionsList[j]->Get("name")->As<std::string>();
        shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ptGrid_->GetNcInterfaceId(ncRegionName));
        if (!ncIf) {
          EXCEPTION("No interface with the name '" << ncRegionName << "' found!");
        }
        shared_ptr<MortarInterface> mortarIf = dynamic_pointer_cast<MortarInterface>(ncIf);
        assert(mortarIf);
        
        if (formulation == "Nitsche") {          
          // obtain a proper scaling for the penalty terms
          std::string nitFac = blochNodesList[i]->Get("nitscheFactor")->As<std::string>();
          Double nitscheFactor = boost::lexical_cast<Double>(nitFac);

          RegionIdType volPrimaryId = mortarIf->GetPrimaryVolRegion();
          RegionIdType volSecondaryId = mortarIf->GetSecondaryVolRegion();
          
          PtrCoefFct permeabilityS, permeabilityP;
          permeabilityS = materials_[volSecondaryId]->GetScalCoefFnc( MAG_PERMEABILITY_SCALAR, 
                                                                                      Global::REAL );                                                                                        
          permeabilityP = materials_[volPrimaryId]->GetScalCoefFnc( MAG_PERMEABILITY_SCALAR, 
                                                                                       Global::REAL ); 
          StdVector<Vector<Double> > points(1);
          Vector<Double> p1(dim_);
          p1.Init();
          points[0]= p1;

          StdVector<Double> valsM, valsS;
          permeabilityS->GetScalarValuesAtCoords(points, valsS, this->ptGrid_);
          permeabilityP->GetScalarValuesAtCoords(points, valsM, this->ptGrid_);
          Double beta = 0.5*( 1.0/valsS[0] + 1.0/valsM[0]) * nitscheFactor;

          // primary & secondary penalty integrals
          BiLinearForm *pnlt_PhiP_PsiP = NULL;
          BiLinearForm *pnlt_PhiP_PsiS = NULL;
          BiLinearForm *pnlt_PhiS_PsiP = NULL;
          BiLinearForm *pnlt_PhiS_PsiS = NULL;
          // primary & secondary integrals with normal derivatives
          BiLinearForm *flux_DPhiP_PsiP = NULL;
          BiLinearForm *flux_DPhiP_PsiS = NULL;
          BiLinearForm *flux_PhiP_DPsiP = NULL;
          BiLinearForm *flux_PhiS_DPsiP = NULL;
          
          shared_ptr<ElemList> actSDList = ncIf->GetElemList();          
          // define bilinear forms for Nitsche coupling
          // penalty integrators   
          if (dim_ == 2) {
            pnlt_PhiP_PsiP = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 2, 1>(),
              new SurfaceIdentityOperator<FeH1, 2, 1>(),
              factor, beta, BiLinearForm::PRIM_PRIM, updatedGeo_, false, true);
            pnlt_PhiP_PsiS = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 2, 1>(),
              new SurfaceIdentityOperator<FeH1, 2, 1>(),
              factorSqr, -beta, BiLinearForm::PRIM_SEC, updatedGeo_, false, true); 
            pnlt_PhiS_PsiP = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 2, 1>(),
              new SurfaceIdentityOperator<FeH1, 2, 1>(),
              constOne, -beta, BiLinearForm::SEC_PRIM, updatedGeo_, false, true);   
            pnlt_PhiS_PsiS = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 2, 1>(),
              new SurfaceIdentityOperator<FeH1, 2, 1>(),
              factor, beta, BiLinearForm::SEC_SEC, updatedGeo_, false, true);   
          }
          else {
            pnlt_PhiP_PsiP = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 3, 1>(),
              new SurfaceIdentityOperator<FeH1, 3, 1>(),
              factor, beta, BiLinearForm::PRIM_PRIM, updatedGeo_, false, true);
            pnlt_PhiP_PsiS = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 3, 1>(),
              new SurfaceIdentityOperator<FeH1, 3, 1>(),
              factorSqr, -beta, BiLinearForm::PRIM_SEC, updatedGeo_, false, true); 
            pnlt_PhiS_PsiP = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 3, 1>(),
              new SurfaceIdentityOperator<FeH1, 3, 1>(),
              constOne, -beta, BiLinearForm::SEC_PRIM, updatedGeo_, false, true);   
            pnlt_PhiS_PsiS = new SurfaceNitscheABInt<Double, Double>(new SurfaceIdentityOperator<FeH1, 3, 1>(),
              new SurfaceIdentityOperator<FeH1, 3, 1>(),
              factor, beta, BiLinearForm::SEC_SEC, updatedGeo_, false, true);     
          }   

          //flux integrators
          PtrCoefFct reluctivityS;
          reluctivityS = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, constOne, permeabilityS, 
                                                 CoefXpr::OP_DIV));                                                          
          
          if (dim_ == 2) {
            flux_DPhiP_PsiP = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceNormalDerivOperator<FeH1,2,1>(),
                                   new SurfaceIdentityOperator<FeH1,2,1>(),
                                   constOne, -1.0, BiLinearForm::PRIM_PRIM, updatedGeo_);
            
            flux_PhiP_DPsiP = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceIdentityOperator<FeH1,2,1>(),
                                   new SurfaceNormalDerivOperator<FeH1,2,1>(),                                   
                                   factor, -1.0, BiLinearForm::PRIM_PRIM, updatedGeo_);

            flux_DPhiP_PsiS = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceNormalDerivOperator<FeH1,2,1>(),
                                   new SurfaceIdentityOperator<FeH1,2,1>(),
                                   factor, 1.0, BiLinearForm::PRIM_SEC, updatedGeo_);

            flux_PhiS_DPsiP = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceIdentityOperator<FeH1,2,1>(),
                                   new SurfaceNormalDerivOperator<FeH1,2,1>(),
                                   constOne, 1.0, BiLinearForm::SEC_PRIM, updatedGeo_);                        
          } else {
            flux_DPhiP_PsiP = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceNormalDerivOperator<FeH1,3,1>(),
                                   new SurfaceIdentityOperator<FeH1,3,1>(),
                                   constOne, -1.0, BiLinearForm::PRIM_PRIM, updatedGeo_);
            
            flux_PhiP_DPsiP = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceIdentityOperator<FeH1,3,1>(),
                                   new SurfaceNormalDerivOperator<FeH1,3,1>(),                                   
                                   factor, -1.0, BiLinearForm::PRIM_PRIM, updatedGeo_);

            flux_DPhiP_PsiS = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceNormalDerivOperator<FeH1,3,1>(),
                                   new SurfaceIdentityOperator<FeH1,3,1>(),
                                   factor, 1.0, BiLinearForm::PRIM_SEC, updatedGeo_);

            flux_PhiS_DPsiP = new SurfaceNitscheABInt<Double, Double>(
                                   new SurfaceIdentityOperator<FeH1,3,1>(),
                                   new SurfaceNormalDerivOperator<FeH1,3,1>(),
                                   constOne, 1.0, BiLinearForm::SEC_PRIM, updatedGeo_);  
          }
         
          // pass material data to the flux operators
          flux_DPhiP_PsiP->SetBCoefFunctionOpA(reluctivityS);
          flux_PhiP_DPsiP->SetBCoefFunctionOpB(reluctivityS);
          flux_DPhiP_PsiS->SetBCoefFunctionOpA(reluctivityS);
          flux_PhiS_DPsiP->SetBCoefFunctionOpB(reluctivityS);
          
          // primary-primary
          pnlt_PhiP_PsiP->SetName("pnlt_PhiP_PsiP");
          flux_DPhiP_PsiP->SetName("flux_DPhiP_PsiP");
          flux_PhiP_DPsiP->SetName("flux_PhiP_DPsiP");
          //primary-secondary
          pnlt_PhiP_PsiS->SetName("pnlt_PhiP_PsiS");
          flux_DPhiP_PsiS->SetName("flux_DPhiP_PsiS");
          //secondary-primary
          pnlt_PhiS_PsiP->SetName("pnlt_PhiS_PsiP");
          flux_PhiS_DPsiP->SetName("flux_PhiS_DPsiP");
          //secondary-secondary
          pnlt_PhiS_PsiS->SetName("pnlt_PhiS_PsiS");
          
          // BiLinearForm::PRIM_PRIM
          SurfaceBiLinFormContext *pnlt_PhiP_PsiP_cont = new SurfaceBiLinFormContext(pnlt_PhiP_PsiP, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext *flux_DPhiP_PsiP_cont = new SurfaceBiLinFormContext(flux_DPhiP_PsiP, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext *flux_PhiP_DPsiP_cont = new SurfaceBiLinFormContext(flux_PhiP_DPsiP, STIFFNESS, BiLinearForm::PRIM_PRIM);
          // BiLinearForm::PRIM_SEC
          SurfaceBiLinFormContext *pnlt_PhiP_PsiS_cont = new SurfaceBiLinFormContext(pnlt_PhiP_PsiS, STIFFNESS, BiLinearForm::PRIM_SEC);
          SurfaceBiLinFormContext *flux_DPhiP_PsiS_cont = new SurfaceBiLinFormContext(flux_DPhiP_PsiS, STIFFNESS, BiLinearForm::PRIM_SEC);
          //BiLinearForm::SEC_PRIM
          SurfaceBiLinFormContext *pnlt_PhiS_PsiP_cont = new SurfaceBiLinFormContext(pnlt_PhiS_PsiP, STIFFNESS, BiLinearForm::SEC_PRIM);
          SurfaceBiLinFormContext *flux_PhiS_DPsiP_cont = new SurfaceBiLinFormContext(flux_PhiS_DPsiP, STIFFNESS, BiLinearForm::SEC_PRIM);
          // BiLinearForm::SEC_SEC
          SurfaceBiLinFormContext *pnlt_PhiS_PsiS_cont = new SurfaceBiLinFormContext(pnlt_PhiS_PsiS, STIFFNESS, BiLinearForm::SEC_SEC);
          
          pnlt_PhiP_PsiP_cont->SetEntities(actSDList, actSDList);
          flux_DPhiP_PsiP_cont->SetEntities(actSDList, actSDList);
          flux_PhiP_DPsiP_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiP_PsiS_cont->SetEntities(actSDList, actSDList);
          flux_DPhiP_PsiS_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiS_PsiP_cont->SetEntities(actSDList, actSDList);
          flux_PhiS_DPsiP_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiS_PsiS_cont->SetEntities(actSDList, actSDList);
          
          pnlt_PhiP_PsiP_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          flux_DPhiP_PsiP_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          flux_PhiP_DPsiP_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          pnlt_PhiP_PsiS_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          flux_DPhiP_PsiS_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          pnlt_PhiS_PsiP_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          flux_PhiS_DPsiP_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          pnlt_PhiS_PsiS_cont->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunctions_[MAG_POTENTIAL]);
          
          assemble_->AddBiLinearForm(pnlt_PhiP_PsiP_cont);
          assemble_->AddBiLinearForm(flux_DPhiP_PsiP_cont);
          assemble_->AddBiLinearForm(flux_PhiP_DPsiP_cont);
          assemble_->AddBiLinearForm(pnlt_PhiP_PsiS_cont);
          assemble_->AddBiLinearForm(flux_DPhiP_PsiS_cont);
          assemble_->AddBiLinearForm(pnlt_PhiS_PsiP_cont);
          assemble_->AddBiLinearForm(flux_PhiS_DPsiP_cont);
          assemble_->AddBiLinearForm(pnlt_PhiS_PsiS_cont);
        } // end Nitsche
        else 
          EXCEPTION("Periodic boundary condition: just Nitsche formulation possoible");
      }
    }  
  }


  LinearForm* MagneticPDE::GetCurrentDensityInt( Double factor,
                                                PtrCoefFct coef, std::string coilId) {
    LinearForm * ret = NULL;
    // -------------------
    //  NODAL FORMULATION
    // -------------------
    // ===  3D CASE ===
    if( dim_ == 3 ) {
      if( isComplex_ ) {
        ret = new BUIntegrator<Complex>( new IdentityOperator<FeH1,3,3,Complex>(),
                                         factor, coef, updatedGeo_, true, false, coilId);
      }
      else {
        ret = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3,Double>(),
                                        factor, coef, updatedGeo_, true, false, coilId);
      }
    } else {
      // ===  2D / AXI CASE ===
      if( isComplex_ ) {
        ret = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1>(),
                                         factor, coef, updatedGeo_, true, false, coilId);
      } else {
        ret = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1>(),
                                        factor, coef, updatedGeo_, true, false, coilId);
      }
    }
    return ret;
  }

  LinearForm* MagneticPDE::GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation ){
    LinearForm * lin = NULL;
    bool coefUpdateGeo = true;

    if(isComplex_) {
      if( dim_ == 2 ) {
        if(isaxi_){
          lin = new BUIntegrator<Complex>( new CurlOperatorAxi<Double>(),
                  Complex(factor),rhsMag,  coefUpdateGeo, fullEvaluation);
        } else {
          lin = new BUIntegrator<Complex>( new CurlOperator<FeH1,2,Double>(),
                  Complex(factor),rhsMag,  coefUpdateGeo, fullEvaluation);
        }
      } else {
        lin = new BUIntegrator<Complex>( new CurlOperator<FeH1,3,Double>(),
                Complex(factor),rhsMag,  coefUpdateGeo, fullEvaluation);
      }
    } else  {
      if( dim_ == 2 ) {
        if(isaxi_){
          // we need +factor as we put -rotM to the rhs
          lin = new BUIntegrator<Double>( new CurlOperatorAxi<Double>(),
                  factor,rhsMag,  coefUpdateGeo, fullEvaluation);
        } else {
          lin = new BUIntegrator<Double>( new CurlOperator<FeH1,2,Double>(),
                  factor,rhsMag,  coefUpdateGeo, fullEvaluation);
        }
      } else {
        lin = new BUIntegrator<Double>( new CurlOperator<FeH1,3,Double>(),
                factor,rhsMag,  coefUpdateGeo, fullEvaluation);
      }
    }

    return lin;
  }

  BaseBDBInt* MagneticPDE::GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity ){
    BaseBDBInt* stiffInt = NULL;

    if( dim_ == 2) {
      if( isaxi_ ) {
        // axisymmetric case
        stiffInt = new BDBInt<>(new CurlOperatorAxi<Double>(), tensorReluctivity, factor, updatedGeo_);
      } else {
        // plane 2D case
        stiffInt = new BDBInt<>(new CurlOperator<FeH1,2,Double>(), tensorReluctivity,factor, updatedGeo_);
      }
    } else {
      // 3D case
      stiffInt = new BDBInt<>(new CurlOperator<FeH1,3,Double>(), tensorReluctivity, factor, updatedGeo_);
    }

    return stiffInt;
  }

  void MagneticPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
            endIt = ncInterfaces_.End();

    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
        case NC_MORTAR:
          if (dim_ == 2) {
            DefineMortarCoupling<2,1>(MAG_POTENTIAL, *ncIt);
            if (isMixed_)// mixed in 2d?
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

    double factor = 1.0;
    if ( isMagnetoStrictCoupled_ == true ){
      factor = -1.0;
    }

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // ==================
    //  COIL INTEGRATORS
    // ==================
	// > see MagBasePDE.cc
    
	
    LinearForm * lin = NULL;
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    bool coefUpdateGeo = true;

    // =================================
    //  Magnetization
    // =================================
    // magnetization should contain ALL regions; hysteresisCoefs only those where hysteresis is defined
    // and mRHSRegions_ only those where constant magnetization is prescribed
    std::map<RegionIdType,PtrCoefFct > regionCoefs = magnetization_->GetRegionCoefs();
    std::map<RegionIdType,PtrCoefFct > regionHystCoefs = hysteresisCoefs_->GetRegionCoefs();
    std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
    for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {

      // get regionIdType
      RegionIdType curReg = it->first;
      PtrCoefFct curHystCoef = regionHystCoefs[curReg];
      PtrCoefFct curFixedMagCoef = mRHSRegions_[curReg];

      // get SDList
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( curReg );

      shared_ptr<CoefFunction> rhsMag;

      // set fullevaluation to trigger evaluation at each integration point
      // the nonlinear parameter "evaluation depth" determines if each
      // integration point gets mapped to midpoint (> fullevaluation = false)
      // or if hyst operator really is evaluated at the actual int. point
      bool fullevaluation = true;
      bool forceSolDependent = false;
      
      if(curHystCoef != NULL){
//        std::cout << "Hysteresis region found" << std::endl;
        // NEW: we do not pass the hysteresis coefficient function
        // directly but instead a special class that returns the
        // correctly weighted term
        // even though, we have a similar function for output,
        // we need a separate coefFunction here as rhsMag might
        // be evaluated at another timestep/interation step as outputMag
        //shared_ptr<CoefFunction> rhsMag = it->second->GenerateRHSCoefFnc("MagMagnetization");
        rhsMag = curHystCoef->GenerateRHSCoefFnc("MagMagnetization");
        mRHSRegions_[curReg] = rhsMag;
        forceSolDependent = true;
      } else if(curFixedMagCoef != NULL) {
        if(magnetizationSet_ == false){
          EXCEPTION("Magnetization has not been initialized yet!");
        }
//        std::cout << "Use constant magnetization " << std::endl;
        rhsMag = mRHSRegions_[curReg];
      } else {
//        std::cout << "Neither hysteresis nor constant magnetization prescribed " << std::endl;
        continue;
      }

      lin = GetRHSMagnetizationInt( factor, rhsMag, fullevaluation );

      lin->SetName("rhs_magnetization");
      if(forceSolDependent){
        lin->SetSolDependent();
      }
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( actSDList );
      ctx->SetFeFunction(feFct);
      assemble_->AddLinearForm(ctx);
      // Add entity list will add nothing, if entities were already assigned
      feFct->AddEntityList(actSDList);

    }
   
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

    ReadRhsExcitation( "fluxDensity", vecComponents, ResultInfo::VECTOR, isComplex_,
            ent, coef, coefUpdateGeo );

    // Problems encountered with fluxDensity:
    // a)
    // if fluxDensity is prescribed it actually is not fixed to that value but it acts solely as excitation
    // > if fluxDensity is set to 0 in input.xml it will not force the flux to be 0 but only apply no excitation
    //    > flux will be different if other sources occur
    // b)
    // if fluxDensity is prescribed along any direction (lets say x) and the same excited block is bounded on top
    // and bottom (assuming 2d here) by fluxParallel boundary conditions in x direction (what should work if the
    // flux is excited in x); the flux inside the block will not at all be parallel to x but circulate; if the boundaries
    // are left default (which should force flux-normal) the flux inside the block look correct (but outside there are
    // additional issues ...)
    // c)
    // if fluxDensity is excited and a hysteretic material is specified, you would expect that H is developing such
    // that mu H + P(H) = B
    // however, as stated above, B is not fixed; instead B acts as excitation from which H arises; H excites P which
    // amplifies the excitation; strangely, this process converges; if the excitation flux is around saturation, the
    // resulting flux will be this excitation flux + saturation
    // > what we can do here is to reduce the rhs-load term by the value of current polarization
    // > we try this via coefFunctionHyst via an additional term which receives the rhsFlux coef as input
    // > first tests seems to suggest that this actually works
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
              ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be defined in a volume.")
      }
      bool forceSolDependency = false;
      shared_ptr<CoefFunction> rhsCoefFunctionToBeUsed;

      // check for hysteresis on current entity
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();

      if(regionCoefs.count(ent[i]->GetRegion())){
        // GenerateRHSCoefFnc("MagPolarization",coef[i]) returns coef[i] - MagPolarization = B - J
        rhsCoefFunctionToBeUsed = regionCoefs[ent[i]->GetRegion()]->GenerateRHSCoefFnc("MagPolarization",coef[i]);
        forceSolDependency = true;
      } else {
        // just take B
        rhsCoefFunctionToBeUsed = coef[i];
      }
      if(isComplex_) {
        if( dim_ == 2) {
          if( isaxi_ ) {
            // axisymmetric case
            lin = new BDUIntegrator<CurlOperatorAxi<Double>, Complex>(Complex(factor),
                    rhsCoefFunctionToBeUsed, relucTensor_, coefUpdateGeo);
          } else {
            lin = new BDUIntegrator<CurlOperator<FeH1,2,Double>,Complex>(Complex(factor),
                    rhsCoefFunctionToBeUsed, relucTensor_,coefUpdateGeo);
          }
        } else {
          // 3D case
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Double>,Complex>(Complex(factor),
                  rhsCoefFunctionToBeUsed, relucTensor_,coefUpdateGeo);
        }
      } else {
        if( dim_ == 2) {
          if( isaxi_ ) {
            // axisymmetric case
            lin = new BDUIntegrator<CurlOperatorAxi<Double>, Double>(factor,  rhsCoefFunctionToBeUsed, relucTensor_,coefUpdateGeo);
          } else {
            lin = new BDUIntegrator<CurlOperator<FeH1,2,Double>,Double>(factor, rhsCoefFunctionToBeUsed, relucTensor_,coefUpdateGeo);
          }
        } else {
          // 3D case
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Double>,Double>(factor, rhsCoefFunctionToBeUsed, relucTensor_,coefUpdateGeo);
        }
      }

      lin->SetName("FluxIntegrator");

      if(forceSolDependency){
        lin->SetSolDependent();
      }
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(feFct);
      assemble_->AddLinearForm(ctx);
      feFct->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];
    } // for


    // ===============
    //  mechStrain (couples in case of magnetostrictive material)
    // ===============
    LOG_DBG(magpde) << "Reading mechanical strain";

    // as we do not have access to the mechanical subtype we have to get the information about
    // the number opf dof from the overall dimensions; luckily planeStrain and planeStress do not differ
    // in case of strainDofNames and tensorType
    SubTensorType tensorType;
    StdVector<std::string> stressDofNames;
    if ( dim_ == 3 ) {
      tensorType = FULL;
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
    } else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
        stressDofNames = "rr", "zz", "rz", "phiphi";
      } else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
        stressDofNames = "xx", "yy", "xy";
      }
    }

    ReadRhsExcitation( "mechStrain", stressDofNames, ResultInfo::TENSOR, isComplex_,
            ent, coef, coefUpdateGeo );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {

      // get region and material from entity
      RegionIdType curRegionId = ent[i]->GetRegion();
      BaseMaterial* curMaterial = NULL;
      bool complexMat = complexMatData_[curRegionId];
      curMaterial = materials_[curRegionId];

      shared_ptr<CoefFunction > curCoef;


      if( complexMat ) {
        //MAGNETOSTRICTION_TENSOR_h_mag is already in a 3 x 6 tensor so we read it in without further transposing (last flag = false)
        curCoef = curMaterial->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h_mag, tensorType,
                Global::COMPLEX, false  );
      } else {
        curCoef = curMaterial->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h_mag, tensorType,
                Global::REAL, false );
      }


      if ( complexMat) {
        // NEW: 3.7.15
        // explanation for factor 1.0 (instead of -1):
        // the used coupling equations for magnotostriction are:
        // H = -hS + nuB
        // sigma = cS - hB
        // Bringing the term -hS to the right hand side would lead to +hS
        // Different from the weak formulation in mechanics, we do NOT loose a -1 factor
        // Reason: Greens Integral theorem for rot(rot H) has a shifted sign compared to the Integral theorem for div(div D)
        // (see "Numerical simulations of electromechanical transducers" Appendix B B47 and B50

        // Apprently, the results show, that it is basically the other way round -> check what the hell is going on here

        Complex couple_factor = Complex(-1.0);

        if( isaxi_ ) {

          lin = new BDUIntegrator<CurlOperatorAxi<Complex>, Complex> (factor*couple_factor,coef[i],curCoef,coefUpdateGeo);

        } else if( dim_ == 2 ) {
          lin = new BDUIntegrator<CurlOperator<FeH1,2,Complex>, Complex> (factor*couple_factor,coef[i],curCoef,coefUpdateGeo);

        } else if( dim_ == 3 ) {
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Complex>, Complex> (factor*couple_factor,coef[i],curCoef,coefUpdateGeo);

        } else {
          EXCEPTION( "unknown dimension in magnetics" );
        }
      }
      else {

        //Double couple_factor = -1.0;
        Double couple_factor = -1.0;

        if( isaxi_ ) {

          lin = new BDUIntegrator<CurlOperatorAxi<Double>, Double> (factor*couple_factor,coef[i],curCoef,coefUpdateGeo);

        } else if( dim_ == 2 ) {
          lin = new BDUIntegrator<CurlOperator<FeH1,2,Double>, Double> (factor*couple_factor,coef[i],curCoef,coefUpdateGeo);

        } else if( dim_ == 3 ) {
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Double>, Double> (factor*couple_factor,coef[i],curCoef,coefUpdateGeo);

        } else {
          EXCEPTION( "unknown dimension in magnetics" );
        }
      }

      lin->SetName("mechStrainInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(feFct);
      assemble_->AddLinearForm(ctx);
      feFct->AddEntityList(ent[i]);
    } // for

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
    res1->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_POTENTIAL);
    res1->entryType = ResultInfo::VECTOR;
    res1->unit = MapSolTypeToUnit(MAG_POTENTIAL);
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


    // === MAGNETIC(?, should be ELECTRIC) SCALAR POTENTIAL ===
    if (isMixed_) {
      shared_ptr<BaseFeFunction> scalFct = feFunctions_[ELEC_POTENTIAL];
      shared_ptr<ResultInfo> res2(new ResultInfo);
      res2->resultType = ELEC_POTENTIAL;
      res2->dofNames = "";
      res2->unit = "V";
      res2->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_POTENTIAL);
      res2->entryType = ResultInfo::SCALAR;
//      results_.Push_back( res2 );
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
      currentInfo->definedOn = ResultInfo::MapSolTypeToDefinedOn(COIL_CURRENT);
      currentInfo->entryType = ResultInfo::SCALAR;

      feFunctions_[COIL_CURRENT]->SetResultInfo(currentInfo);
      DefineFieldResult( feFunctions_[COIL_CURRENT], currentInfo );
    }

    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability ( new ResultInfo );
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_ELEM_PERMEABILITY);
    permeability->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult(permFct, permeability);
    availResults_.insert(permeability);

	  //creates the velocity
    StdVector<std::string> vecDofNames;
    if( ptGrid_->GetDim() == 3 ) {
      vecDofNames = "x", "y", "z";
    } else {
      if( ptGrid_->IsAxi() ) {
        vecDofNames = "r", "z";
      } else {
        vecDofNames = "x", "y";
      }
    }

    //// === VELOCITY ===
    shared_ptr<ResultInfo> velocity( new ResultInfo);
    velocity->resultType = MEAN_FLUIDMECH_VELOCITY;
    velocity->dofNames = vecDofNames;
    velocity->unit = "m/s";

    velocity->definedOn = ResultInfo::MapSolTypeToDefinedOn(MEAN_FLUIDMECH_VELOCITY);
    velocity->entryType = ResultInfo::VECTOR;

    VelocityCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    DefineFieldResult( VelocityCoef_, velocity );

    results_.Push_back( velocity );
    availResults_.insert( velocity );
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
      aDot->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_POTENTIAL_DERIV1);
      aDot->entryType = ResultInfo::VECTOR;
      availResults_.insert( aDot );
      DefineTimeDerivResult( MAG_POTENTIAL_DERIV1, 1, MAG_POTENTIAL );

      // === COIL CURRENT OF COILS EXCITED BY VOLTAGE - 1ST DERIVATIVE ===
      if( hasVoltCoils_ ){
        shared_ptr<ResultInfo> iDot(new ResultInfo);
        iDot->resultType = COIL_CURRENT_DERIV1;
        iDot->dofNames = "";
        iDot->unit = "A/s";
        iDot->definedOn = ResultInfo::MapSolTypeToDefinedOn(COIL_CURRENT_DERIV1);
        iDot->entryType = ResultInfo::SCALAR;
        availResults_.insert( iDot );
        DefineTimeDerivResult( COIL_CURRENT_DERIV1, 1, COIL_CURRENT );
      }
    }

    // from NACS
    // === ELECTRIC SCALAR POTENTIAL - 1ST DERIVATIVE===
    if (isMixed_) {
      shared_ptr<ResultInfo> phiDot(new ResultInfo);
      phiDot->resultType = ELEC_POTENTIAL_DERIV_1;
      phiDot->dofNames = "";
      phiDot->unit = "V/s";
      phiDot->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_POTENTIAL_DERIV_1);
      phiDot->entryType = ResultInfo::SCALAR;
      availResults_.insert( phiDot );
      DefineTimeDerivResult( ELEC_POTENTIAL_DERIV_1, 1, ELEC_POTENTIAL );
    }

    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD;
    rhs->dofNames = aVecComponents;
    rhs->unit = "";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_RHS_LOAD);
    rhsFeFunctions_[MAG_POTENTIAL]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL], rhs );

    if(!fluxDensityDefined_){
      DefineMagFluxDensity();
    }
//    shared_ptr<CoefFunctionFormBased>
    PtrCoefFct bFunc = this->GetCoefFct(MAG_FLUX_DENSITY);

//    // === MAGNETIC FLUX DENSITY ===
//    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
//    fluxDens->resultType = MAG_FLUX_DENSITY;
//    fluxDens->dofNames = vecComponents;
//    fluxDens->unit = "Vs/m^2";
//    fluxDens->definedOn = ResultInfo::ELEMENT;
//    fluxDens->entryType = ResultInfo::VECTOR;
//    availResults_.insert( fluxDens );
//    shared_ptr<CoefFunctionFormBased> bFunc;
//    if( isComplex_ ) {
//      bFunc.reset(new CoefFunctionBOp<Complex>(feFct, fluxDens));
//    } else {
//      bFunc.reset(new CoefFunctionBOp<Double>(feFct, fluxDens));
//    }
//    DefineFieldResult( bFunc, fluxDens );
//    stiffFormCoefs_.insert(bFunc);

    // === MAGNETIC NORMAL FLUX DENSITY ===
    shared_ptr<ResultInfo> normFlux(new ResultInfo);
    normFlux->resultType = MAG_NORMAL_FLUX_DENSITY;
    normFlux->dofNames = "";
    normFlux->unit = "Vs/m^2";
    normFlux->entryType = ResultInfo::SCALAR;
    normFlux->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_NORMAL_FLUX_DENSITY);
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
    flux->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_FLUX);
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

    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> magIntens ( new ResultInfo );
    magIntens->resultType = MAG_FIELD_INTENSITY;
    magIntens->dofNames = vecComponents;
    //magIntens->SetVectorDOFs(dim_, isaxi_);
    magIntens->unit = "A/m";
    magIntens->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_FIELD_INTENSITY);
    magIntens->entryType = ResultInfo::VECTOR;

    /*
     * TODO: check what happens if we have hysteretic and non-hysteretic regions
     * > in this case isHysteresis_ should be true, so that in the non-hysteretic regions
     *   no H-Field will be output?! > actually yes
     * > solution idea: add belows definition for magIntensFunc in MagBasePDE where fieldIntensity_ is filled
     *   there we check each region for hysteresis and currently add a coefFunction only on hysteretic regions
     *   
     */
//    if ( !isHysteresis_){
////      shared_ptr<CoefFunctionFormBased> magIntensFunc;
////        if( isComplex_ ) {
////          magIntensFunc.reset(new CoefFunctionFlux<Complex>(feFct, magIntens));
////        } else {
////          magIntensFunc.reset(new CoefFunctionFlux<Double>(feFct, magIntens));
////        }
////        
////        // subtract fixed magnetization!
////      // cannot be done here as we do not have access to different regions
////        CoefXprBinOp sub(mp_, magIntensFunc, mRHSRegions_[actRegion], CoefXpr::OP_SUB);
////        PtrCoefFct magIntensFuncCorrected(CoefFunction::Generate(mp_, part, sub));
////
////        fieldIntensity_->AddRegion( actRegion, magIntensFuncCorrected);
////        std::string warnmsg = "NON-Hysteretic H set on region " + regionName;
////        std::cout << warnmsg << std::endl;
////        //DefineFieldResult( magIntensFunc, magIntens );
////        stiffFormCoefs_.insert(magIntensFunc);
//
//      shared_ptr<CoefFunctionFormBased> magIntensFunc;
//      if( isComplex_ ) {
//        magIntensFunc.reset(new CoefFunctionFlux<Complex>(feFct, magIntens));
//      } else {
//        magIntensFunc.reset(new CoefFunctionFlux<Double>(feFct, magIntens));
//      }
//      DefineFieldResult( magIntensFunc, magIntens );
//      stiffFormCoefs_.insert(magIntensFunc);
//    }

    // Magnetization, Polarization and Field Intensity
//    if ( isHysteresis_){
      shared_ptr<ResultInfo> magJ ( new ResultInfo );
      magJ->resultType = MAG_POLARIZATION;
      magJ->SetVectorDOFs(dim_, isaxi_);
      magJ->unit = "Vs/m^2";
      magJ->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_POLARIZATION);
      magJ->entryType = ResultInfo::VECTOR;

      DefineFieldResult( polarization_, magJ );
      availResults_.insert( magJ );

      shared_ptr<ResultInfo> magM ( new ResultInfo );
      magM->resultType = MAG_MAGNETIZATION;
      magM->SetVectorDOFs(dim_, isaxi_);
      magM->unit = "A/m";
      magM->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_MAGNETIZATION);
      magM->entryType = ResultInfo::VECTOR;

      DefineFieldResult( magnetization_, magM );
      availResults_.insert( magM );

      DefineFieldResult( fieldIntensity_, magIntens );
      availResults_.insert( magIntens );
//    }

    // for both BdBKernel and EnergyResultFunctor, we need to apply the -1 factor
    // to get right sign in the results (even though the energy results are not really usable in the coupled case as they neglect the influnce of the coupled pde)
    Double factor = 1.0;
    if ( isMagnetoStrictCoupled_ ){
      factor = -1.0;
    }

    // === EDDY CURRENT DENSITY ===
    shared_ptr<CoefFunctionFormBased> jFunc;
    if( analysistype_ != STATIC ) {
      shared_ptr<BaseFeFunction> aDotFct =
              timeDerivFeFunctions_[MAG_POTENTIAL_DERIV1];
      shared_ptr<ResultInfo> eddy(new ResultInfo);
      eddy->resultType = MAG_EDDY_CURRENT_DENSITY;
      eddy->dofNames = aVecComponents;
      eddy->unit = "A/m^2";
      eddy->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_EDDY_CURRENT_DENSITY);
      eddy->entryType = ResultInfo::VECTOR;
      availResults_.insert( eddy );

      // from NACS
      if( isComplex_ ) {
        jFunc.reset(new CoefFunctionFlux<Complex>(aDotFct, eddy, -factor));
      } else {
        jFunc.reset(new CoefFunctionFlux<Double>(aDotFct, eddy, -factor));
      }
      massFormCoefs_.insert(jFunc);


      // === EDDY CURRENT (SURFACE RESULT) ===
      shared_ptr<ResultInfo> ec(new ResultInfo());
      ec->resultType = MAG_EDDY_CURRENT;
      ec->dofNames = "";
      ec->unit = "A";
      ec->definedOn = ResultInfo::REGION;
      ec->entryType = ResultInfo::SCALAR;
      availResults_.insert( ec );
      // first, create normal mapping
      //shared_ptr<CoefFunctionSurf> ncd(new CoefFunctionSurf(true, 1.0, ec));
      //surfCoefFcts_[ncd] = jFunc;
      // then, integrate values
      shared_ptr<ResultFunctor> eddyCurrentFunc;
      if( isComplex_ ) {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Complex>(jFunc,
                                                                  feFct, ec ) );
      } else {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Double>(jFunc,
                                                           feFct, ec ) );
      }
      resultFunctors_[MAG_EDDY_CURRENT] = eddyCurrentFunc;


      /*
       * currently, the following part leads to an Exception (when hex elements are used in 3d)
       * > ElemType is not supported by FESpace
       * CoupledField::BaseFE* CoupledField::FeSpaceH1Nodal::GetFe(CoupledField::UInt): requested FEType (11) is not supported by space
       * In file '/home/userHome/staff/mnierla/CFS/GIT/CFS/source/FeBasis/H1/FeSpaceH1Nodal.cc' at line 145
       * > error occurs only if result is to be written out;
       */
      if( isMixed_)  {
        shared_ptr<BaseFeFunction> phiDotFct =
                    timeDerivFeFunctions_[ELEC_POTENTIAL_DERIV_1];
        // 2nd part: - gamma * d(grad phi)/dt
        shared_ptr<CoefFunctionFormBased> part2;
        if(isComplex_) {
          part2.reset(new CoefFunctionFlux<Complex>(phiDotFct, eddy, -factor));
        } else {
          part2.reset(new CoefFunctionFlux<Double>(phiDotFct, eddy, -factor));
        }
        mixedFormFunctor_.insert(part2);

        // sum up both parts
        CoefXprBinOp sum(mp_, jFunc, part2, CoefXpr::OP_ADD);
        PtrCoefFct jFinalFunc(CoefFunction::Generate(mp_, part, sum));
        DefineFieldResult( jFinalFunc, eddy );
      } else
      {
        DefineFieldResult( jFunc, eddy );
      }

      // === EDDY POWER DENSITY ===
      shared_ptr<ResultInfo> epd(new ResultInfo());
      epd->resultType = MAG_EDDY_POWER_DENSITY;
      epd->dofNames = "";
      epd->unit = "W/m^3";
      epd->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_EDDY_POWER_DENSITY);
      epd->entryType = ResultInfo::SCALAR;

      // from NACS
      if( isMixed_) {
       // in this case we really have to compute the power density
       //    w_eddy = J_eddy^2 / gamma
       // as element result, as J_eddy is composed of of the A- and
       // phi-related part.
        PtrCoefFct jCoef = GetCoefFct(MAG_EDDY_CURRENT_DENSITY);
        CoefXprBinOp mult(mp_, jCoef, jCoef, CoefXpr::OP_MULT_CONJ);
        PtrCoefFct jTmp = CoefFunction::Generate(mp_, part, mult);
        CoefXprBinOp epdExpr(mp_, jTmp, conduc_, CoefXpr::OP_DIV);
        PtrCoefFct epdFct= CoefFunction::Generate(mp_, part, epdExpr);
        DefineFieldResult( epdFct, epd );
      } else {
        shared_ptr<CoefFunctionFormBased> epdFunctor;
        if( isComplex_ ) {
          epdFunctor.reset( new CoefFunctionBdBKernel<Complex>(aDotFct, factor));
        } else {
          epdFunctor.reset( new CoefFunctionBdBKernel<Double>(aDotFct, factor));
        }
        DefineFieldResult( epdFunctor, epd );
        massFormCoefs_.insert(epdFunctor);
      }

//      // === EDDY POWER ===
//      shared_ptr<ResultInfo> ep(new ResultInfo());
//      ep->resultType = MAG_EDDY_POWER;
//      ep->dofNames = "";
//      ep->unit = "W";
//      ep->definedOn = ResultInfo::REGION;
//      ep->entryType = ResultInfo::SCALAR;
//      availResults_.insert( ep );
//      if( isMixed_)
//        WARN("Adjust eddy power for mixed case");
//      shared_ptr<ResultFunctor> epFunctor;
//      if( isComplex_ ) {
//        epFunctor.reset(new EnergyResultFunctor<Complex>(aDotFct, ep, factor));
//      } else {
//        epFunctor.reset(new EnergyResultFunctor<Double>(aDotFct, ep, factor));
//      }
//      resultFunctors_[MAG_EDDY_POWER] = epFunctor;
//      massFormFunctors_.insert(epFunctor);

      // === EDDY POWER ===
      shared_ptr<ResultInfo> ep(new ResultInfo());
      ep->resultType = MAG_EDDY_POWER;
      ep->dofNames = "";
      ep->unit = "W";
      ep->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_EDDY_POWER);
      ep->entryType = ResultInfo::SCALAR;
      availResults_.insert( ep );
      shared_ptr<ResultFunctor> epFunctor;

      // from NACS
      if( isMixed_) {
        // Here we have to integrate the eddy current density
        PtrCoefFct eddyDens = GetCoefFct(MAG_EDDY_POWER_DENSITY);
        if( isComplex_ ) {
          epFunctor.reset(new ResultFunctorIntegrate<Complex>(eddyDens, aDotFct,
                                                              ep));
        } else {
          epFunctor.reset(new ResultFunctorIntegrate<Double>(eddyDens, aDotFct, ep));
        }
      } else {

        if( isComplex_ ) {
          epFunctor.reset(new EnergyResultFunctor<Complex>(aDotFct, ep, factor));
        } else {
          epFunctor.reset(new EnergyResultFunctor<Double>(aDotFct, ep, factor));
        }
        massFormFunctors_.insert(epFunctor);
      }
      resultFunctors_[MAG_EDDY_POWER] = epFunctor;

      // === COIL LINKED FLUX - 1ST DERIVATIVE, CORRESPONDS TO INDUCED VOLTAGE ===
      shared_ptr<ResultInfo> psiDotRes(new ResultInfo());
      psiDotRes->resultType = COIL_INDUCED_VOLTAGE;
      psiDotRes->dofNames = "";
      psiDotRes->unit = "V";
      psiDotRes->definedOn = ResultInfo::MapSolTypeToDefinedOn(COIL_INDUCED_VOLTAGE);
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


      // === ELECTRIC FIELD INTENSITY ===
      shared_ptr<ResultInfo> elecIntens(new ResultInfo);
      elecIntens->resultType = ELEC_FIELD_INTENSITY;
      elecIntens->SetVectorDOFs(dim_, isaxi_);
      elecIntens->dofNames = vecComponents;
      elecIntens->unit = "V/m";
      elecIntens->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_FIELD_INTENSITY);
      elecIntens->entryType = ResultInfo::VECTOR;

      // assemble coefficient function E = - dA/dt
      PtrCoefFct eIFuncTmp = CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,  aDotFct, aDotFct, CoefXpr::OP_MULT ) ); // not used?!
      PtrCoefFct eIFunc = CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,  "-1.0", aDotFct, CoefXpr::OP_MULT ) );
      DefineFieldResult( eIFunc, elecIntens);
    }

    // determine dimensionality of current density
    UInt jDim = aVecComponents.GetSize();

    // === COIL CURRENT DENSITY ===
    shared_ptr<ResultInfo> ccd(new ResultInfo);
    ccd->resultType = MAG_COIL_CURRENT_DENSITY;
    ccd->dofNames = aVecComponents;
    ccd->unit = "A/m^2";
    ccd->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_COIL_CURRENT_DENSITY);
    ccd->entryType = ResultInfo::VECTOR;
    availResults_.insert( ccd );
    shared_ptr<CoefFunctionMulti> ccdCoef(new CoefFunctionMulti(CoefFunction::VECTOR, jDim, 1,isComplex_));
    DefineFieldResult( ccdCoef, ccd );

    // === TOTAL CURRENT DENSITY ===
    shared_ptr<ResultInfo> tcd(new ResultInfo);
    tcd->resultType = MAG_TOTAL_CURRENT_DENSITY;
    tcd->dofNames = aVecComponents;
    tcd->unit = "A/m^2";
    tcd->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_TOTAL_CURRENT_DENSITY);
    tcd->entryType = ResultInfo::VECTOR;
    availResults_.insert( tcd );
    shared_ptr<CoefFunctionMulti> tcdCoef(new CoefFunctionMulti(CoefFunction::VECTOR,jDim,1,
            isComplex_));
    DefineFieldResult( tcdCoef, tcd );

    GenerateLorentzForceResults(vecComponents, tcdCoef, bFunc, part, feFct);

 	// === MAXWELL FORCE DENSITY ===
	shared_ptr<ResultInfo> mfd(new ResultInfo);
	mfd->resultType = MAG_FORCE_MAXWELL_DENSITY;
	mfd->dofNames = vecComponents;
	mfd->unit = MapSolTypeToUnit(mfd->resultType);
	mfd->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_FORCE_MAXWELL_DENSITY);
	mfd->entryType = ResultInfo::VECTOR;
	availResults_.insert( mfd );

	// Note: The positive normal direction in this case is defined as the
	//       inward facing one.
	shared_ptr<CoefFunctionSurfMaxwell> maxForceDens(new CoefFunctionSurfMaxwell(false, matCoefs_, ptGrid_, -1.0, mfd));
	DefineFieldResult(maxForceDens, mfd);
	surfCoefFcts_[maxForceDens] = bFunc;

  // === MAXWELL NORMAL_FORCE DENSITY ===
  shared_ptr<ResultInfo> mfdN(new ResultInfo);
  mfdN->resultType = MAG_NORMALFORCE_MAXWELL_DENSITY;
  mfdN->dofNames = vecComponents;
  mfdN->unit = MapSolTypeToUnit(mfdN->resultType);
  mfdN->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_NORMALFORCE_MAXWELL_DENSITY);
  mfdN->entryType = ResultInfo::VECTOR;
  availResults_.insert( mfdN );
  shared_ptr<CoefFunctionSurfMaxwell> maxNormalForceDens(new CoefFunctionSurfMaxwell(false,
                                                         matCoefs_, ptGrid_, -1.0, mfdN));
  DefineFieldResult(maxNormalForceDens, mfdN);
  surfCoefFcts_[maxNormalForceDens] = bFunc;

  // === MAXWELL TANGENTIAL_FORCE DENSITY ===
  shared_ptr<ResultInfo> mfdT(new ResultInfo);
  mfdT->resultType = MAG_TANGENTIALFORCE_MAXWELL_DENSITY;
  mfdT->dofNames = vecComponents;
  mfdT->unit = MapSolTypeToUnit(mfdT->resultType);
  mfdT->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_TANGENTIALFORCE_MAXWELL_DENSITY);
  mfdT->entryType = ResultInfo::VECTOR;
  availResults_.insert( mfdT );
  shared_ptr<CoefFunctionSurfMaxwell> maxTangentialForceDens(new CoefFunctionSurfMaxwell(false,
                                                             matCoefs_, ptGrid_, -1.0, mfdT));
  DefineFieldResult(maxTangentialForceDens, mfdT);
  surfCoefFcts_[maxTangentialForceDens] = bFunc;


  if( analysistype_ != HARMONIC ) {
	// === MAXWELL FORCE (TOTAL) ===
	shared_ptr<ResultInfo> mf(new ResultInfo);
	mf->resultType = MAG_FORCE_MAXWELL;
	mf->dofNames = vecComponents;
	mf->unit = "N";
	mf->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_FORCE_MAXWELL);
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
	vwp->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_FORCE_VWP);
	vwp->entryType = ResultInfo::VECTOR;
	availResults_.insert( vwp );

	// define and save coefFunction
	shared_ptr<CoefFunctionSurfVWP> vwpForce(new CoefFunctionSurfVWP(false, matCoefs_,
              1.0, vwp));
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
//    if( isMixed_)
//      WARN("Adjust energy for mixed case.");
    shared_ptr<ResultInfo> energy(new ResultInfo);
    energy->resultType = MAG_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_ENERGY);
    energy->entryType = ResultInfo::SCALAR;
    availResults_.insert( energy );
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy, 0.5*factor));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy, 0.5*factor));
    }
    resultFunctors_[MAG_ENERGY] = energyFunc;
    stiffFormFunctors_.insert(energyFunc);


    // === COIL LINKED FLUX  ===
    shared_ptr<ResultInfo> psiRes(new ResultInfo());
    psiRes->resultType = COIL_LINKED_FLUX;
    psiRes->dofNames = "";
    psiRes->unit = "Vs/m^2";
    psiRes->definedOn = ResultInfo::MapSolTypeToDefinedOn(COIL_LINKED_FLUX);
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

    // === COIL INDUCTANCE  ===
    // more or less the same as COIL LINKED FLUX with the only difference
    // that during the integration the factor 1/coil_current is applied
    shared_ptr<ResultInfo> indRes(new ResultInfo());
    indRes->resultType = COIL_INDUCTANCE;
    indRes->dofNames = "";
    indRes->unit = "Vs/A";
    indRes->definedOn = ResultInfo::MapSolTypeToDefinedOn(COIL_INDUCTANCE);
    indRes->entryType = ResultInfo::SCALAR;
    availResults_.insert( indRes );
    shared_ptr<ResultFunctor> indFunc;
    shared_ptr<CoefFunctionMulti> indDens(new CoefFunctionMulti(CoefFunction::SCALAR,1,1,
            isComplex_));
    if( isComplex_ ){
      indFunc.reset( new ResultFunctorIntegrate<Complex>(indDens, feFct, indRes) );
    } else {
      indFunc.reset( new ResultFunctorIntegrate<Double>(indDens, feFct, indRes) );
    }
    resultFunctors_[COIL_INDUCTANCE] = indFunc;
    // it is an integrated result but we need to save the coef function
    // somewhere for the finalization
    fieldCoefs_[COIL_INDUCTANCE] = indDens;

    // === CORE LOSS DENSITY ===
    shared_ptr<ResultInfo> cldRes(new ResultInfo());
    cldRes->resultType = MAG_CORE_LOSS_DENSITY;
    cldRes->dofNames = "";
    cldRes->unit = "W/kg";
    cldRes->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_CORE_LOSS_DENSITY);
    cldRes->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> coreLossDensCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
    DefineFieldResult( coreLossDensCoef, cldRes );

    // === CORE LOSS ===
    shared_ptr<ResultInfo> clRes(new ResultInfo());
    clRes->resultType = MAG_CORE_LOSS;
    clRes->dofNames = "";
    clRes->unit = "W";
    clRes->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_CORE_LOSS);
    clRes->entryType = ResultInfo::SCALAR;
    //DefineFieldResult( coreLossCoef, clRes );
    availResults_.insert( clRes );
    shared_ptr<ResultFunctor> coreLossFunc;
    shared_ptr<CoefFunctionMulti> coreLossCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
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
    shared_ptr<ResultInfo> jld(new ResultInfo);
    jld->resultType = MAG_JOULE_LOSS_POWER_DENSITY;
    jld->dofNames = "";
    jld->unit = "W/m^3";
    jld->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_JOULE_LOSS_POWER_DENSITY);
    jld->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> jldCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
    DefineFieldResult( jldCoef, jld );

    // optimization results are provided in DesignSpace::ExtractResults()
    DefineFieldResult(PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);
    DefineFieldResult(PHYSICAL_PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);
    DefineFieldResult(RHS_PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);
    DefineFieldResult(PHYSICAL_RHS_PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);

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

    StdVector<RegionIdType>::iterator regIt;

    // === TOTAL CURRENT DENSITY ===
    PtrCoefFct jEddy = GetCoefFct(MAG_EDDY_CURRENT_DENSITY);
    shared_ptr<CoefFunctionMulti> tcdCoef
            = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
    // loop over all regions and assemble total current density:
    //  - if region is coil -> take coil current
    //  - if region is no coil and analyis is transient/harmonic -> eddy

    regIt = regions_.Begin();

//    if( isMixed_){
//      WARN("Debug TOTAL CURRENT DENSITY for mixed case");
//    } else {
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
//      }
    }

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // === CORE LOSS DENSITY ===
    // This is a "density" per kg, not per m^3. That's why we cannot integrate it over
    // the volume to get the core loss (see below).
    shared_ptr<CoefFunctionMulti> cldCoef =
            dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_CORE_LOSS_DENSITY]);
    BaseMaterial* actMat = NULL;
    regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ) {
      RegionIdType actRegion = *regIt;
      actMat = materials_[actRegion];
      PtrCoefFct coreLossFct = actMat->GetScalCoefFncNonLin( MAG_CORE_LOSS_PER_MASS, Global::REAL, GetCoefFct(MAG_FLUX_DENSITY) );
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
      PtrCoefFct coreLossFct = actMat->GetScalCoefFncNonLin( MAG_CORE_LOSS_PER_MASS, Global::REAL,
              GetCoefFct(MAG_FLUX_DENSITY) );
      // core loss density in W/m^3, which deserves to be called density
      PtrCoefFct densFct = actMat->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct coreLossDensCoef = CoefFunction::Generate( mp_, Global::REAL,
              CoefXprBinOp( mp_, coreLossFct, densFct, CoefXpr::OP_MULT ));
      clCoef->AddRegion( actRegion, coreLossDensCoef );
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
    //
    // The coil inductance is COIL_LINKED_FLUX / COIL_CURRENT
    // integrate this at the same time
    // Note that this is only the static inductance, which is quite useless for
    // circuits with non-linear material or permanent magnets!
    PtrCoefFct temp = GetCoefFct(COIL_LINKED_FLUX);
    PtrCoefFct temp_ind = GetCoefFct(COIL_INDUCTANCE);

    shared_ptr<CoefFunctionMulti> psiDotDens =
            dynamic_pointer_cast<CoefFunctionMulti>(temp);

    shared_ptr<CoefFunctionMulti> indDens =
            dynamic_pointer_cast<CoefFunctionMulti>(temp_ind);

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

      // once again with 1/COIL_CURRENT
      if( coilRegsIt->second->sourceType_ == Coil::CURRENT ){
        // it won't work for voltage coils until the coilCurrentDens_ is fixed (see DefineIntegrators)

        PtrCoefFct coilCurrentDensity;
        PtrCoefFct coilCurrent;

        // this gives us just the current density! we have to multiply with the winding cross section!
        CoefXprUnaryOp coilCurrentDensityOp = CoefXprUnaryOp(mp_, coilCurrentDens_[coilRegsIt->first],CoefXpr::OP_NORM);
        coilCurrentDensity = CoefFunction::Generate( mp_, part, coilCurrentDensityOp );

        CoefXprBinOp coilCurrentOp = CoefXprBinOp(mp_, coilCurrentDensity,
                boost::lexical_cast<std::string>(partIt->second->wireCrossSect), CoefXpr::OP_MULT);

        coilCurrent = CoefFunction::Generate( mp_, part, coilCurrentOp );

        CoefXprBinOp indIntegrandOp = CoefXprBinOp( mp_, integrand,
                coilCurrent, CoefXpr::OP_DIV );

        PtrCoefFct indIntegrand = CoefFunction::Generate( mp_, part, indIntegrandOp );

        indDens->AddRegion( coilRegsIt->first, indIntegrand );
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
                CoefXprBinOp( mp_, materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL),
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

  void MagneticPDE::FinalizeAfterTimeStep() {}


  std::map<SolutionType, shared_ptr<FeSpace> >
  MagneticPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    // and 3D analysis
    isMixed_ = false;
    CheckForConductivity();

    if( analysistype_ != BasePDE::STATIC && dim_ == 3 ) {
      // Mixed formulation ONLY required if at least one region has conductivity; otherwise
      // algsys fails due to empty blocks
      if(anyRegionHasConductivity_)
        isMixed_ = true;
    }

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;

    if( formulation == "default" || formulation == "H1" ){

      // 1) create space for magnetic vector potential
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] =
              FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);

      // 2) create space for electric scalar potential
      if( isMixed_ ) {
        crSpaces[ELEC_POTENTIAL] = FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
        crSpaces[ELEC_POTENTIAL]->Init(solStrat_);
      }

    } else {
      EXCEPTION( "The formulation " << formulation << "of the magnetic PDE is not known!" );
    }

    // ===================================================================
    // Create FeSpaceConst for coil currents (of coils excited by voltage)
    // ===================================================================
    if( hasVoltCoils_ ){
      PtrParamNode voltSpaceNode = infoNode->Get("coilCurrent");
      crSpaces[COIL_CURRENT] = FeSpace::CreateInstance(myParam_, voltSpaceNode, FeSpace::CONSTANT, ptGrid_);
      crSpaces[COIL_CURRENT]->Init(solStrat_);
    }
    return crSpaces;
  }

} // end namespace CoupledField
