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
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
//#include "Materials/Models/VectorPreisach.hh"
#include "Domain/CoefFunction/CoefFunctionHyst.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"

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
    
    updatedGeo_        = !true; //! true needed for NC interface to give correct results while moving
    isMagnetoStrictCoupled_ = false;
    mechanicPDE_ = NULL;
    
    regionApproxSet_ = false;
    
    // can the reluctivity be complex? before the change it had the same type as the PDE
    reluc_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_, dim_, false));
    
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
  
  void MagneticPDE::SetMagnetoStrictCoupling(SinglePDE *mechanicPDE)
  {
    mechanicPDE_ = mechanicPDE;
    isMagnetoStrictCoupled_ = true;
    
  }
  
  shared_ptr<Coil> MagneticPDE::GetCoilById(const Coil::IdType& id) {
    return coils_.at(id);
  }
  
  void MagneticPDE::InitNonLin()
  {
    
    SinglePDE::InitNonLin();
  }
  
  void MagneticPDE::InitHystCoefs() {
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    SubTensorType tensorType;
    
    if ( dim_ == 3  || subType_ == "2.5d" ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }
    
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[MAG_POTENTIAL]->GetFeSpace();
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;
      
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // ==========================================================
      // New implementation
      // ==========================================================
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId, integId);
      
      SDLists_[actRegion] = actSDList;
      
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
      if ( nonLinTypes.Find(HYSTERESIS) != -1 ){
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( actRegion );
        
        // we set the dependency of the hyst operator to the actual
        // quantity that is computed (B) although the standard Preisach
        // as well as the VectorPreisach model need H as input
        PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FLUX_DENSITY);

        // please note: 
        //  in magnetics, the hysteresis model is supposed to return the
        //  magnetic polarizaiton J_P = mu*M
        //  (in older versions, the magnetization M was returned!)
        PtrCoefFct hystPol(new CoefFunctionHyst( actSDMat, actSDList,
                magFieldCoef,tensorType,MAG_RELUCTIVITY,mySpace));
        
        hysteresisCoefs_->AddRegion( actRegion, hystPol);
        
      }
    }
    regionApproxSet_ = true;
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
    
	  for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
		  // Set current region and material
		  actRegion = it->first;
		  actMat = it->second;
      
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
		  StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
      
		  // ====================================================================
		  //  NONLINEAR BH RELATION (NON-HYSTERETIC)
		  // ====================================================================
		  if (  nonLinTypes.Find(PERMEABILITY) != -1 ) {
		    CoefFunctionOpt* cfo = NULL; // we might do optimization and then we have such a thing
			  PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
			  PtrCoefFct nuNl = actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY, Global::REAL, magFluxCoef);

        if(domain->HasDesign())
        {
          cfo = new CoefFunctionOpt(domain->GetDesign(), nuNl, this);
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
				  stiff2->SetNewtonBilinearForm();
          
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

        if ( nonLinTypes.Find(HYSTERESIS) != -1 ){
          /* for both the delta material method as well as the std fixpoint method we have to know
           * which regions are affected by hysteresis
           */
          // NEW: coefFncHyst should already be created during DefinePostProcResults!      
          PtrCoefFct hystPol = hysteresisCoefs_->GetRegionCoef(actRegion);
          curCoef = hystPol->GenerateMatCoefFnc("Reluctivity");
          
          PtrCoefFct hystOutput = hystPol->GenerateOutputCoefFnc("MagPolarization");
          hysteresisPolarization_->AddRegion( actRegion, hystOutput);
          
          PtrCoefFct hystOutput2 = hystPol->GenerateOutputCoefFnc("MagMagnetization");
          hysteresisMagnetization_->AddRegion( actRegion, hystOutput2);
          
			  } else{
				  // ====================================================================
				  //  Standard Linear CASE (2D AND 3D)
				  // ====================================================================
				  curCoef = actMat->GetTensorCoefFnc( MAG_RELUCTIVITY, tensorType, Global::REAL );

				  // for postprocessing
				  PtrCoefFct permeability = materials_[actRegion]->GetScalCoefFnc( MAG_PERMEABILITY, Global::REAL);
			    if(domain->HasDesign())
			    {
			      cfo = new CoefFunctionOpt(domain->GetDesign(), curCoef, this);
			      curCoef.reset(cfo);
			    }
				  matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);
			  }
        
			  BaseBDBInt * stiffInt = NULL;
        
			  if( dim_ == 2) {
				  if( isaxi_ ) {
					  // axisymmetric case
					  stiffInt = new BDBInt<>(new CurlOperatorAxi<Double>(), curCoef, factor, updatedGeo_);
				  } else {
					  // plane 2D case
					  stiffInt = new BDBInt<>(new CurlOperator<FeH1,2,Double>(), curCoef,factor, updatedGeo_);
				  }
			  } else {
				  // 3D case
				  stiffInt = new BDBInt<>(new CurlOperator<FeH1,3,Double>(), curCoef, factor, updatedGeo_);
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
			  bdbInts_[actRegion] = stiffInt;
        
        // NEW: curCoef is a specialized coef function responsible
        // for delivering the mat tensor > type is tensor
        //			  // add also material to global, distributed reluctivity coefficient function
        //			  if ( nonLinTypes.Find(HYSTERESIS) != -1){// || nonLinTypes.Find(HYSTERESIS_FIXPOINT) != -1 ){
        //			    /*
        //			     * we cannot directly add coefFunctionHyst to reluc_ as reluc_ expects tensorial coefFncs
        //			     * but coefFunctionHyst has to be a vector coefFnc
        //			     */
        //			  } else {
			  reluc_->AddRegion(actRegion, curCoef);
        //			  }
        
			  // ====================================================================
			  //  3D CASE (additional stiffness integrator)
			  // ====================================================================
			  if( dim_ == 3 ) {
				  BaseBDBInt * divInt =
                  new BDBInt<>(new DivOperator<FeH1,3,Double>(), curCoef, factor, updatedGeo_);
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
      
		  materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
		  PtrCoefFct conducCoef =
              materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY,Global::REAL);
		  //                                         lexical_cast<std::string>(conductivity));
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
                          new BBInt<>(new GradientOperator<FeH1,3,1,Double>(), conducCoef, factor, updatedGeo_);
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
                          new GradientOperator<FeH1,3,1,Double>(), conducCoef, factor, updatedGeo_);
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
				  LinearForm* curInt = NULL;
          
				  // generate source current vector
				  CoefFunctionOpt* cfoc = NULL; // we might do optimization and then we have such a thing
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

          /*if(domain->HasDesign())
          {
            cfoc = new CoefFunctionOpt(domain->GetDesign(), jFct, this);
            jFct.reset(cfoc);
          }*/
          
				  if( dim_ == 3 ) {
					  // ===========
            //  3D CASE
            // ===========
            if( isComplex_ ) {
              curInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,3,3,Complex>(),
                      factor, jFct, updatedGeo_);
            }
            else {
              curInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3,Double>(),
                      factor, jFct, updatedGeo_);
            }
            
				  } else {
					  // ===============
					  //  2D / AXI CASE
					  // ===============
            
					  if( isComplex_ ) {
						  curInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1>(),
                      factor, jFct, updatedGeo_);
					  } else {
						  curInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1>(),
                      factor, jFct, updatedGeo_);
					  }
				  }
          
				  curInt->SetName("CoilIntegrator");
				  LinearFormContext * coilContext =
                  new LinearFormContext( curInt );
				  coilContext->SetEntities( actSDList );
				  coilContext->SetFeFunction( myFct );
				  assemble_->AddLinearForm( coilContext );

	        // when we have a CoefFunctionOpt, we tell it the proper form, which we only have now
	        if(cfoc)
	          cfoc->SetForm(curInt);
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
				  if( dim_ == 3 ) {
					  if( isComplex_ ) {
						  psiDotInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,3,3,Complex>(),
                      -1.0*factor, eJscaled, updatedGeo_);
					  } else {
						  psiDotInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3,Double>(),
                      -1.0*factor, eJscaled, updatedGeo_);
					  }
				  } else {
					  if( isComplex_ ) {
						  psiDotInt = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1,Complex>(),
                      -1.0*factor, eJscaled, updatedGeo_);
					  } else {
						  psiDotInt = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1,Double>(),
                      -1.0*factor, eJscaled, updatedGeo_);
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
                        1.0*factor, eJscaled, updatedGeo_);
						  } else {
							  psiDotIntT = new BUIntegrator<Double>( new IdentityOperator<FeH1,3,3,Double>(),
                        1.0*factor, eJscaled, updatedGeo_);
						  }
					  } else {
						  if( isComplex_ ) {
							  psiDotIntT = new BUIntegrator<Complex>( new IdentityOperator<FeH1,2,1,Complex>(),
                        1.0*factor, eJscaled, updatedGeo_);
						  } else {
							  psiDotIntT = new BUIntegrator<Double>( new IdentityOperator<FeH1,2,1,Double>(),
                        1.0*factor, eJscaled, updatedGeo_);
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
    
    double factor = 1.0;
    if ( isMagnetoStrictCoupled_ == true ){
      factor = -1.0;  
    }
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
    
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
    
    LinearForm * lin = NULL;
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    bool coefUpdateGeo = true;
    
    // =================================
    //  Magnetization -> from hysteresis (VOLUME)
    // =================================
    //check for hysteresis
    if ( isHysteresis_ ){
      //    std::cout << "Putting magnetization to rhs" << std::endl;
      LOG_DBG(magpde) << "Putting magnetization to rhs";
      
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
        
        // get regionIdType
        RegionIdType curReg = it->first;
        
        // get SDList
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( curReg );
        
        // set fullevaluation to trigger evaluation at each integration point
        // the nonlinear parameter "evaluation depth" determines if each
        // integration point gets mapped to midpoint (> fullevaluation = false)
        // or if hyst operator really is evaluated at the actual int. point
        bool fullevaluation = true;
        
        // NEW: we do not pass the hysteresis coefficient function
        // directly but instead a special class that returns the
        // correctly weighted term
        // even though, we have a similar function for output,
        // we need a separate coefFunction here as rhsMag might
        // be evaluated at another timestep/interation step as outputMag
        shared_ptr<CoefFunction> rhsMag = it->second->GenerateRHSCoefFnc("MagMagnetization");
        
        //std::cout << "Put magnetization on rhs" << std::endl;
        //factor = factor*(-1.0);
        if(isComplex_) {
          if( dim_ == 2 ) {
            if(isaxi_){
              lin = new BUIntegrator<Complex>( new CurlOperatorAxi<Double>(),
                      Complex(factor),rhsMag,  coefUpdateGeo, fullevaluation);
            } else {
              lin = new BUIntegrator<Complex>( new CurlOperator<FeH1,2,Double>(),
                      Complex(factor),rhsMag,  coefUpdateGeo, fullevaluation);
            }
          } else {
            lin = new BUIntegrator<Complex>( new CurlOperator<FeH1,3,Double>(),
                    Complex(factor),rhsMag,  coefUpdateGeo, fullevaluation);
          }
        } else  {
          if( dim_ == 2 ) {
            if(isaxi_){
              // we need +factor as we put -rotH to the rhs
              lin = new BUIntegrator<Double>( new CurlOperatorAxi<Double>(),
                      factor,rhsMag,  coefUpdateGeo, fullevaluation);
            } else {
              lin = new BUIntegrator<Double>( new CurlOperator<FeH1,2,Double>(),
                      factor,rhsMag,  coefUpdateGeo, fullevaluation);
            }
          } else {
            lin = new BUIntegrator<Double>( new CurlOperator<FeH1,3,Double>(),
                    factor,rhsMag,  coefUpdateGeo, fullevaluation);
          }
        }
        
        lin->SetName("rhs_magnetization");
        lin->SetSolDependent();
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( actSDList );
        ctx->SetFeFunction(feFct);
        assemble_->AddLinearForm(ctx);
        // Add entity list will add nothing, if entities were already assigned
        feFct->AddEntityList(actSDList);
      }
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
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
              ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be defined in a volume.")
      }
      
      if(isComplex_) {
        if( dim_ == 2) {
          if( isaxi_ ) {
            // axisymmetric case
            lin = new BDUIntegrator<CurlOperatorAxi<Double>, Complex>(Complex(factor), 
                    coef[i], reluc_, coefUpdateGeo);
          } else {
            lin = new BDUIntegrator<CurlOperator<FeH1,2,Double>,Complex>(Complex(factor), 
                    coef[i], reluc_,coefUpdateGeo);
          }
        } else {
          // 3D case
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Double>,Complex>(Complex(factor), 
                  coef[i], reluc_,coefUpdateGeo);
        }
      } else {
        if( dim_ == 2) {
          if( isaxi_ ) {
            // axisymmetric case
            lin = new BDUIntegrator<CurlOperatorAxi<Double>, Double>(factor,  coef[i], reluc_,coefUpdateGeo);
          } else {
            lin = new BDUIntegrator<CurlOperator<FeH1,2,Double>,Double>(factor, coef[i], reluc_,coefUpdateGeo);
          }
        } else {
          // 3D case
          lin = new BDUIntegrator<CurlOperator<FeH1,3,Double>,Double>(factor, coef[i], reluc_,coefUpdateGeo);
        }  
      }
      
      lin->SetName("FluxIntegrator");
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
  
  void MagneticPDE::ReadSpecialBCs() {
    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();
  }
  
  
  
  void MagneticPDE::DefineSolveStep()
  {
		if(isHysteresis_){
			solveStep_ = new SolveStepHyst(*this);
		} else {
			solveStep_ = new StdSolveStep(*this);
		}
    //solveStep_ = new StdSolveStep(*this);
  }
  
  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MagneticPDE::InitTimeStepping() {
    // Use complete implicit scheme
    Double gamma = 1;
    GLMScheme * scheme = new Trapezoidal(gamma);
    
    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
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
      // WARN("Adapt printing of coils to InfoNode");
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
    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
    fluxDens->resultType = MAG_FLUX_DENSITY;
    fluxDens->dofNames = vecComponents;
    fluxDens->unit = "Vs/m^2";
    fluxDens->definedOn = ResultInfo::ELEMENT;
    fluxDens->entryType = ResultInfo::VECTOR;
    availResults_.insert( fluxDens );
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

    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> magIntens ( new ResultInfo );
    magIntens->resultType = MAG_FIELD_INTENSITY;
    magIntens->SetVectorDOFs(dim_, isaxi_);
    magIntens->unit = "A/m";
    magIntens->definedOn = ResultInfo::ELEMENT;
    magIntens->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> magIntensFunc;
    if( isComplex_ ) {
      magIntensFunc.reset(new CoefFunctionFlux<Complex>(feFct, magIntens));
    } else {
      magIntensFunc.reset(new CoefFunctionFlux<Double>(feFct, magIntens));
    }
    stiffFormCoefs_.insert(magIntensFunc);
    
    // Magnetization and Polarization and Field Intensity
    if ( isHysteresis_){
      // this is the point where the actual storage is created
      hysteresisCoefs_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      hysteresisPolarization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      hysteresisMagnetization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      
      shared_ptr<ResultInfo> magJ ( new ResultInfo );
      magJ->resultType = MAG_POLARIZATION;
      magJ->SetVectorDOFs(dim_, isaxi_);
      magJ->unit = "Vs/m^2";
      magJ->definedOn = ResultInfo::ELEMENT;
      magJ->entryType = ResultInfo::VECTOR;
      
      DefineFieldResult( hysteresisPolarization_, magJ );
      availResults_.insert( magJ );
      
      shared_ptr<ResultInfo> magM ( new ResultInfo );
      magM->resultType = MAG_MAGNETIZATION;
      magM->SetVectorDOFs(dim_, isaxi_);
      magM->unit = "A/m";
      magM->definedOn = ResultInfo::ELEMENT;
      magM->entryType = ResultInfo::VECTOR;

      DefineFieldResult( hysteresisMagnetization_, magM );
      availResults_.insert( magM );
      
      // new: init hysteresis right now (so that it is available
      // in case of magnetostrictive coupling for the definition of
      // the integrators
      InitHystCoefs();
      
      
//      // we should use the same nu here that is used for inversion etc
//      // otherwise results do not match!
//      PtrCoefFct magMag;
//      PtrCoefFct nu0 = CoefFunction::Generate( mp_, Global::REAL, "795774.7155");
//
//      magMag = CoefFunction::Generate(mp_,Global::REAL,CoefXprBinOp(mp_,nu0,hysteresisPolarization_,CoefXpr::OP_MULT));
//      DefineFieldResult( magMag, magM );
//      availResults_.insert( magM );
      
      // Field intensity
      PtrCoefFct magIntensFuncMod;
      /*
       * in case of hysteresis, H has is computed from B using nu0;
       * however, B = nu0^-1 H + nu0^-1 M
       * i.e.
       * H = nu0*B - M
       * i.e.
       * subtract M
       */
      if( isComplex_ ) {
        magIntensFuncMod = CoefFunction::Generate(mp_,Global::COMPLEX,CoefXprBinOp(mp_,magIntensFunc,hysteresisMagnetization_,CoefXpr::OP_SUB));
      } else {
        magIntensFuncMod = CoefFunction::Generate(mp_,Global::REAL,CoefXprBinOp(mp_,magIntensFunc,hysteresisMagnetization_,CoefXpr::OP_SUB));
      }
      DefineFieldResult( magIntensFuncMod, magIntens );
      
    } else {
      // use std computation of mag intensity, i.e. H = B/mu
      DefineFieldResult( magIntensFunc, magIntens );
    }
    
    // for both BdBKernel and EnergyResultFunctor, we need to apply the -1 factor
    // to get right sign in the results (even though the energy results are not really usable in the coupled case as they neglect the influnce of the coupled pde)
    Double factor = 1.0;
    if ( isMagnetoStrictCoupled_ ){
      factor = -1.0;
    }
    
    if( analysistype_ != STATIC ) {
      // === EDDY CURRENT DENSITY ===
      shared_ptr<BaseFeFunction> aDotFct =
              timeDerivFeFunctions_[MAG_POTENTIAL_DERIV1];
      shared_ptr<ResultInfo> eddy(new ResultInfo);
      eddy->resultType = MAG_EDDY_CURRENT_DENSITY;
      eddy->dofNames = "";
      eddy->unit = "A/m^2";
      eddy->definedOn = ResultInfo::ELEMENT;
      eddy->entryType = ResultInfo::VECTOR;
      availResults_.insert( eddy );
      shared_ptr<CoefFunctionFormBased> jFunc;
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
        epdFunctor.reset( new CoefFunctionBdBKernel<Complex>(aDotFct, factor));
      } else {
        epdFunctor.reset( new CoefFunctionBdBKernel<Double>(aDotFct, factor));
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
        epFunctor.reset(new EnergyResultFunctor<Complex>(aDotFct, ep, factor));
      } else {
        epFunctor.reset(new EnergyResultFunctor<Double>(aDotFct, ep, factor));
      }
      resultFunctors_[MAG_EDDY_POWER] = epFunctor;
      massFormFunctors_.insert(epFunctor);
      
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
    
    // === COIL INDUCTANCE  ===
    // more or less the same as COIL LINKED FLUX with the only difference
    // that during the integration the factor 1/coil_current is applied
    shared_ptr<ResultInfo> indRes(new ResultInfo());
    indRes->resultType = COIL_INDUCTANCE;
    indRes->dofNames = "";
    indRes->unit = "Vs/A";
    indRes->definedOn = ResultInfo::COIL;
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
    cldRes->definedOn = ResultInfo::ELEMENT;
    cldRes->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> coreLossDensCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
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
    jld->definedOn = ResultInfo::ELEMENT;
    jld->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> jldCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
    DefineFieldResult( jldCoef, jld );


    // optimization results are provided in DesignSpace::ExtractResults()
    // copied from MechPDE
    // === MECH_PSEUDO_DENISTY ===
    shared_ptr<ResultInfo> mpd(new ResultInfo);
    mpd->resultType = MECH_PSEUDO_DENSITY;
    mpd->entryType = ResultInfo::SCALAR;
    mpd->definedOn = ResultInfo::ELEMENT;
    mpd->dofNames = "";
    mpd->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mpd); // the fe-function is only a dummy

    // === PHYSICAL_PSEUDO_DENISTY ===
    shared_ptr<ResultInfo> ppd(new ResultInfo);
    ppd->resultType = PHYSICAL_PSEUDO_DENSITY;
    ppd->entryType = ResultInfo::SCALAR;
    ppd->definedOn = ResultInfo::ELEMENT;
    ppd->dofNames = "";
    ppd->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), ppd);

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
      PtrCoefFct coreLossFct = actMat->GetScalCoefFncNonLin( CORE_LOSS, Global::REAL, GetCoefFct(MAG_FLUX_DENSITY) );
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
  
  void MagneticPDE::FinalizeAfterTimeStep() {
    
    //	  //check for hysteresis
    //	  if ( isHysteresis_ && isHysteresisFixPoint_ == false ) {
    //		  //set current values to previous values for hysteresis operator
    //		  //needed for the next time step
    //		  std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
    //		  std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
    //		  for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
    //			  it->second->SetPreviousHystVals();
    //		  }
    //	  }
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
      crSpaces[COIL_CURRENT] = FeSpace::CreateInstance(myParam_, voltSpaceNode, FeSpace::CONSTANT, ptGrid_);
      crSpaces[COIL_CURRENT]->Init(solStrat_);
    }
    return crSpaces;
  }
  
} // end namespace CoupledField
