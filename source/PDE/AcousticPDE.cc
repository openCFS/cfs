#include "AcousticPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/DivOperator.hh"

#include "FeBasis/FeFunctions.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionCurvilinearPML.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionImpedanceModel.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>


#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

   DEFINE_LOG(acousticpde, "pde.acoustic")


  AcousticPDE::AcousticPDE( Grid* aGrid, PtrParamNode paramNode,
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "acoustic";
    pdematerialclass_  = ACOUSTIC;
    nonLin_            = false;
    isMechCoupled_     = false;
    formulation_ = ACOU_PRESSURE;

    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;
    isTimeDomPML_      = false;
    isAPML_            = false;
    complexFluidFormulation_ = false;

    std::string pdeFormulation = myParam_->Get("formulation")->As<std::string>();
    //check for pressure or potential formulation
    sosAtLaplace_ = (pdeFormulation == "acouPressureSOSatLaplace")? true : false;

    if(pdeFormulation != "default" && pdeFormulation != "acouPressureSOSatLaplace"){
      formulation_ = SolutionTypeEnum.Parse(pdeFormulation);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> > AcousticPDE::CreateFeSpaces( const std::string&  formulation,
                  PtrParamNode infoNode ){

    if(this->analysistype_ == STATIC)
      EXCEPTION("No STATIC analysis in AcousticPDE");

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      std::string form = SolutionTypeEnum.ToString(formulation_);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[formulation_] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[formulation_]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of acoustic PDE is not known!");
    }

    // ===================================
    // Check for transient PML
    // ===================================
    if(this->analysistype_ == TRANSIENT && isTimeDomPML_){
      //now define the additional unknowns

      if(dim_==3 && !this->isAPML_){
        PtrParamNode scalarpml = infoNode->Get("TransientPMLScalarAuxVar");
        crSpaces[ACOU_PMLAUXSCALAR] =
            FeSpace::CreateInstance(myParam_,scalarpml,FeSpace::H1, ptGrid_);
        crSpaces[ACOU_PMLAUXSCALAR]->Init(solStrat_);
      }

      PtrParamNode vectorPML = infoNode->Get("TransientPMLVectorAuxVar");
      crSpaces[ACOU_PMLAUXVEC] =
          FeSpace::CreateInstance(myParam_,vectorPML,FeSpace::H1, ptGrid_);
      crSpaces[ACOU_PMLAUXVEC]->Init(solStrat_);

    }
    return crSpaces;
  }
  
  void AcousticPDE::ReadDampingInformation() {
    std::map<std::string, DampingType> idDampType;
    std::map<std::string, shared_ptr<RaylDampingData> > idRaylData;

    // try to get dampingList
    PtrParamNode dampListNode = myParam_->Get( "dampingList", ParamNode::PASS );
    if( dampListNode ) {

      // get specific damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();

      for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->As<std::string>();

        // determine type of damping
        DampingType actType;
        String2Enum( dampString, actType );

        if( actType == RAYLEIGH ) {
          // set data for Rayleigh damping
          shared_ptr<RaylDampingData> actRaylDamp(new RaylDampingData());
          actRaylDamp->alpha = "0.0";
          actRaylDamp->beta = "0.0";
          actRaylDamp->adjustDamping = true;
          actRaylDamp->ratioDeltaF = 0.01;
          actRaylDamp->freq = 0.0;
          dampNodes[i]->GetValue( "freq", actRaylDamp->freq, ParamNode::PASS);
          dampNodes[i]->GetValue( "ratioDeltaF", actRaylDamp->ratioDeltaF,
                                  ParamNode::PASS );
          dampNodes[i]->GetValue( "adjustDamping",actRaylDamp->adjustDamping,  
                                  ParamNode::PASS);
          idRaylData[actId] = actRaylDamp;        
        }

        // store damping type string
        idDampType[actId] = actType;

      }
    }

    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes =
        myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptGrid_->GetRegion().Parse( actRegionName );

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];
      if ( dampingList_[actRegionId] == RAYLEIGH ){
        RaylDampingData actRayl = *(idRaylData[actDampingId]);
        Double dampFreq;

        if( actRayl.freq == 0.0 ) {
          materials_[actRegionId]->GetScalar(dampFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        } else { 
          dampFreq = actRayl.freq;
        }
        // Compute Rayleigh damping parameters
        materials_[actRegionId]->
        ComputeRayleighDamping( actRayl.alpha, actRayl.beta,
                                dampFreq, actRayl.ratioDeltaF, 
                                actRayl.adjustDamping, isComplex_ );
        regionRaylDamping_[actRegionId] = actRayl;
      } else if(dampingList_[actRegionId] == PML &&
          analysistype_ == BasePDE::TRANSIENT ) {
        isTimeDomPML_ = true;
      }
      
    }
  }

  void AcousticPDE::DefineIntegrators(){

    RegionIdType actRegion;
    //type of geometry
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    LOG_DBG(acousticpde) << "DefineIntegrators BEGIN" <<  "\n" ;

    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++) {
      actRegion = regions_[iRegion];

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);

      //=======================================================================
      // Generate coefficient functions 4 (Speed of Sound)
      //=======================================================================
      PtrCoefFct c0;
      PtrCoefFct dens;
      PtrCoefFct blk;

      if ( complexFluidFormulation_ ) {
    	  dens = materials_[actRegion]->GetScalCoefFnc( DENSITY, Global::COMPLEX );
    	  blk = materials_[actRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::COMPLEX );
      }
      else {
    	  dens = materials_[actRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
    	  blk = materials_[actRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
      }

      LOG_DBG(acousticpde) << "DefineIntegrators : dens   = " << dens->ToString() << "\n";
      LOG_DBG(acousticpde) << "DefineIntegrators : blk    = " << blk->ToString() << "\n";


      // ====================================================================
      // check for temperature field, which effects the speed of sound
      // ====================================================================
      std::string tempId = curRegNode->Get("temperatureId")->As<std::string>();
      PtrCoefFct regionTemp;

      if ( tempId != "" ) {
    	  std::set<UInt> definedDofs;
    	  //just real valued temperature allowed
    	  shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);
    	  //Add the region information
    	  PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature","name",tempId.c_str());

    	  ReadUserFieldValues( actSDList, tempNode, tempInfo->dofNames, tempInfo->entryType,
                               false, regionTemp, definedDofs, updatedGeo_ );
    	  meanTemperatureCoef_->AddRegion( actRegion, regionTemp );
      }

      //compute speed of sound
      ComputeSOS(c0, dens, blk, regionTemp, tempId );
      // store coefficient functions
      matCoefs_[ELEM_DENSITY]->AddRegion(actRegion, dens);
      matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->AddRegion( actRegion, c0);

      // if pde couples with mechanic, we have to multiply the density by -1
      PtrCoefFct factor;
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      if ( isMechCoupled_ == true && formulation_ != ACOU_PRESSURE ) {
        if ( complexFluidFormulation_ )
          EXCEPTION("Complex fluid and coupled mechanical-acoustic simulation not allowed");

        // Important: In case of a general / quadratic EV problem, we must
        // ensure to have a "positive definite" matrix, i.e. we are not allowed
        // to multiply all matrices by -1!
        std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";

        factor = CoefFunction::Generate( mp_, Global::REAL,
                                        CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT ) );
      } else {
        factor = constOne;
      }

      //basic coeff-functions for mass and stiffness matrix
      PtrCoefFct coeffM, coeffK;

      if ( formulation_ == ACOU_PRESSURE && sosAtLaplace_ == true ) {
        if ( complexFluidFormulation_ )
          EXCEPTION("A complex fluid and sosAtLaplace-formulation not allowed!!");

        //pressure formulation with temperature depend speed of sound
        coeffM = constOne;
        coeffK = CoefFunction::Generate( mp_,  Global::REAL,
                               CoefXprBinOp(mp_, c0, c0, CoefXpr::OP_MULT ) );
      }
      else {
    	  if ( complexFluidFormulation_ ) {
    		  //in this case c0 is actually 1/c0^2!!
    		  //! coeffM: 1/compressionModulus
    		  //! coeffK = 1/density
    		  coeffM = CoefFunction::Generate( mp_, Global::COMPLEX,
    		      				  CoefXprBinOp(mp_, factor, blk, CoefXpr::OP_DIV ) );
    		  ///coeffM: 1/density
    		  coeffK = CoefFunction::Generate( mp_, Global::COMPLEX,
    				  CoefXprBinOp(mp_, factor, dens, CoefXpr::OP_DIV ) );
    	  }
    	  else {
    		  coeffK = factor;
    		  // build coefficient for mass matrix as (factor / (c0*c0))
    		  PtrCoefFct constValOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
    		  PtrCoefFct coeffa, coeffb;
    		  coeffa = CoefFunction::Generate( mp_,  Global::REAL,
    				                           CoefXprBinOp(mp_, constValOne, c0, CoefXpr::OP_DIV ) );
    		  coeffb = CoefFunction::Generate( mp_,  Global::REAL,
          		                               CoefXprBinOp(mp_, coeffa, coeffa, CoefXpr::OP_MULT ) );
    		  coeffM = CoefFunction::Generate( mp_, Global::REAL,
                		      CoefXprBinOp(mp_, factor, coeffb, CoefXpr::OP_MULT) );
    	  }
      }

      LOG_DBG(acousticpde) << "DefineIntegrators Fluid: coeffK = " << coeffK->ToString() << "\n";
      LOG_DBG(acousticpde) << "DefineIntegrators Fluid: coeffM = " << coeffM->ToString() << "\n";

      BaseBDBInt *stiffInt = NULL;
      BaseBDBInt *massInt = NULL;
        // ====================================================================
      // PML integrators
      // ====================================================================
      if (dampingList_[actRegion] == PML)
        DefinePMLIntegrators(actRegion, actSDList, curRegNode, c0, coeffK, coeffM, tempId, massInt, stiffInt);
      else {
        // ====================================================================
        // standard stiffness integrator
        // ====================================================================
        if( dim_ == 2 ) {
          if ( complexFluidFormulation_ ) {
            stiffInt = new BBInt<Complex>(new GradientOperator<FeH1,2>(), coeffK, 1.0, updatedGeo_ );
          } else {
            stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), coeffK, 1.0, updatedGeo_ );
         }
        } else { //dim_ == 3
          if ( complexFluidFormulation_ ) {
            stiffInt = new BBInt<Complex>(new GradientOperator<FeH1,3>(), coeffK, 1.0, updatedGeo_ );
          } else {
            stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), coeffK, 1.0, updatedGeo_ );
          }
      }

      // ====================================================================
      // standard mass integrator
      // ====================================================================
        if(dim_==2) {
          if ( complexFluidFormulation_ )
            massInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1,Double>,coeffM, 1.0, updatedGeo_ );
          else
            massInt = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>,coeffM, 1.0, updatedGeo_ );
        } else { // dim_==3
          if  ( complexFluidFormulation_ )
            massInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1,Double>, coeffM, 1.0, updatedGeo_ );
          else
            massInt = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>, coeffM, 1.0, updatedGeo_ );
        }
      }

      // ====================================================================
      // stiffness integrator context
      // ====================================================================
      stiffInt->SetName("LaplaceIntegrator");
      BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );

      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH ) {
        if ( complexFluidFormulation_ )
          EXCEPTION("Complex fluid region and Rayleigh damping not allowed!!");
        RaylDampingData & actDamp = (regionRaylDamping_[actRegion]);
        stiffIntDescr->SetSecDestMat(DAMPING, actDamp.beta );
      }

      feFunctions_[formulation_]->AddEntityList( actSDList );
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
      stiffInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

      // ====================================================================
      // mass integrator context
      // ====================================================================
      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace() );

      BiLinFormContext *massContext =  new BiLinFormContext(massInt, MASS );
      // Check for damping (mass part)
      if ( dampingList_[actRegion] == RAYLEIGH ) {
        if ( complexFluidFormulation_ )
          EXCEPTION("Complex fluid region and Rayleigh damping not allowed!!");
        RaylDampingData & actDamp = regionRaylDamping_[actRegion];
        massContext->SetSecDestMat( DAMPING, actDamp.alpha );
      }

      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
      assemble_->AddBiLinearForm( massContext );
      // Important: Add mass-integrator to global list, as we need them later
      // for calculation of postprocessing results
      massInts_[actRegion] = massInt;

      // ====================================================================
      // flow integrators
      // ====================================================================
      if (dim_ == 2) 
      {
        if (isComplex_) 
        {
          DefineConvectiveIntegrators<2, true>(actRegion, curRegNode, actSDList, coeffM);
        }

        else 
        {
          DefineConvectiveIntegrators<2, false>(actRegion, curRegNode, actSDList, coeffM);
        }
      } 
      
      else
      { /* if (dim_ == 3) */
        if (isComplex_)
        {
          DefineConvectiveIntegrators<3, true>(actRegion, curRegNode, actSDList, coeffM);
        } 

        else 
        {
          DefineConvectiveIntegrators<3, false>(actRegion, curRegNode, actSDList, coeffM);
        }
      }
    }
  }


  void AcousticPDE::DefinePMLIntegrators(RegionIdType actRegion, shared_ptr<ElemList>& actSDList, PtrParamNode& curRegNode,
                                         PtrCoefFct& c0, PtrCoefFct& coeffK, PtrCoefFct& coeffM, std::string& tempId,
                                         BaseBDBInt*& stiffInt,  BaseBDBInt*& massInt) {
    //variables used for PML .. mit nach PML
    std::string pmlDampId = "";  // dampId of a PML region
    std::string pmlFormul = "";  // formulation of the PML region ("classic", "shifted" or "curvilinear")
    PtrCoefFct coeffPMLDeterminant, coeffPMLTensor, coeffPMLVector, coeffPMLStiff, coeffPMLMass;

    // retrieve the dampingId of the current PML region
    curRegNode->GetValue("dampingId", pmlDampId);

    // check for PML formulation of current region: classic(=Cartesian) or curvilinear
    PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", pmlDampId.c_str());
    pmlFormul = pmlNode->Get("formulation")->As<std::string>();

    // Check if the analysis type is harmonic or inverse-source
    if (analysistype_ == HARMONIC || analysistype_ == BasePDE::INVERSESOURCE) {
      // For complex fluids, compute a real-valued speed of sound for PML definition
      shared_ptr<CoefFunction> densR, blkR, c0R;
      if (complexFluidFormulation_) {
        densR = materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
        blkR = materials_[actRegion]->GetScalCoefFnc(ACOU_BULK_MODULUS, Global::REAL);
        c0R = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blkR, densR, CoefXpr::OP_DIV), CoefXpr::OP_SQRT));
      } else {
        c0R = c0;
      }
      if (pmlFormul == "classic") {
        // here we only have a diagonal Jacobi matrix, so we store values as vector
        coeffPMLVector.reset(new CoefFunctionPML<Complex>(pmlNode, c0R, actSDList, regions_, true));
        coeffPMLDeterminant.reset(new CoefFunctionPML<Complex>(pmlNode, c0R, actSDList, regions_, false));
        // store pml factor for the postprocessing result
        matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLVector);
        // compute factors for the integrators
        coeffPMLStiff = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, coeffPMLDeterminant, coeffK, CoefXpr::OP_MULT));
        coeffPMLMass = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, coeffPMLDeterminant, coeffM, CoefXpr::OP_MULT));
        // + classic PML ints (stiff + mass)
        if (dim_ == 2) {
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,2,Complex>(),
                                        coeffPMLStiff, 1.0, updatedGeo_ );
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1,Complex>(), 
                                       coeffPMLMass, 1.0, updatedGeo_ );
        } else { // 3D
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,3,Complex>(),
                                        coeffPMLStiff, 1.0, updatedGeo_ );
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1,Complex>(), 
                                       coeffPMLMass, 1.0, updatedGeo_ );
        }
        // set coordinate stretching coefFunction
        stiffInt->SetBCoefFunctionOpB(coeffPMLVector);
      }
      else if (pmlFormul == "curvilinear") {
        if (complexFluidFormulation_)
          WARN("Curvilinear PML with complex fluid is not tested. Check results and consider adding a testcase.");
        if (dim_ == 2)
          EXCEPTION("Curvilinear PML currently only implemented in 3D!");

        // pointer to object that handles the computation of the curvilinear PML damping tensor
        // here we need the full Jacobi matrix
        coeffPMLTensor.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::TENSOR));
        coeffPMLDeterminant.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::DETERMINANT));

        // create some more CoefFunctionCurvilinearPMLs to store info about the PML parameters
        PtrCoefFct coeffPMLDampFactor, coeffPMLDistance;
        coeffPMLDampFactor.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::DAMP_FACTOR));
        coeffPMLDistance.reset(new CoefFunctionCurvilinearPML<Complex>(pmlNode, c0R, actSDList, regions_, OutputType::DISTANCE));

        // assign the coefFunctions to the matCoefs_ to make them available in the DefinePostProcResults()
        matCoefs_[PML_TENSOR]->AddRegion(actRegion, coeffPMLTensor);
        matCoefs_[PML_DETERMINANT]->AddRegion(actRegion, coeffPMLDeterminant);
        matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLDampFactor);
        matCoefs_[PML_DISTANCE]->AddRegion(actRegion, coeffPMLDistance);

        // the Jakobi determinant is used as scaling coefficient for the BBIntegrators
        coeffPMLStiff = CoefFunction::Generate(mp_, Global::COMPLEX,
          CoefXprBinOp(mp_, coeffPMLDeterminant, coeffK, CoefXpr::OP_MULT));
        coeffPMLMass = CoefFunction::Generate(mp_, Global::COMPLEX,
          CoefXprBinOp(mp_, coeffPMLDeterminant, coeffM, CoefXpr::OP_MULT));
        // + cuvilinear PML ints (stiff + mass)
        if (dim_ == 2) {
          // add 2D implementation here =)
        } else { // 3D
          // define integrators for curvilinear PML in 3D
          stiffInt = new BBInt<Complex>(new TensorScaledGradientOperator<FeH1,3,Complex>(),
                                        coeffPMLStiff, 1.0, updatedGeo_ );
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1,Complex>(), 
                                       coeffPMLMass, 1.0, updatedGeo_ );
        }
        // set coordinate stretching coefFunction
        stiffInt->SetBCoefFunctionOpB(coeffPMLTensor);
      } else { // when pmlFormul is invalid...
        EXCEPTION("Unknown PML formulation '" << pmlFormul << "' for AcousticPDE. Possible formulations: 'classic', 'curvilinear'.")
      }
    } else { // if(analysistype_ == HARMONIC || analysistype_ == BasePDE::INVERSESOURCE)
      if (pmlFormul == "classic") {
        if (dim_ == 2) {
          DefineTransientPMLInts<2>(actSDList, pmlDampId, actRegion, tempId); // add here
          // + standard stiff + mass ints
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), coeffK,
                                        1.0, updatedGeo_ );
          massInt = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>,coeffM, 
                                        1.0, updatedGeo_ );
        } else {
          DefineTransientPMLInts<3>(actSDList, pmlDampId, actRegion, tempId);
          // + standard stiff + mass ints
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), coeffK,
                                        1.0, updatedGeo_ );
          massInt = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>,coeffM, 
                                        1.0, updatedGeo_ );
        }
      } else {
        EXCEPTION("Transient PML is currently only implemented in 'classic' formulation.")
      }
    }
  }


  template<UInt DIM>
  void AcousticPDE::DefineTransientPMLInts(shared_ptr<ElemList> eList, std::string id,
		                                        RegionIdType actRegion, std::string tempId){

    //define some material coeffunction as above...
    PtrCoefFct dens = materials_[eList->GetRegion()]->GetScalCoefFnc( DENSITY, Global::REAL );
    PtrCoefFct blk = materials_[eList->GetRegion()]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
    PtrCoefFct one = CoefFunction::Generate( mp_, Global::REAL, "1.0");

    PtrCoefFct mechAcouFactor;
    if ( isMechCoupled_ == true && formulation_ != ACOU_PRESSURE ) {
      // Important: In case of a general / quadratic EV problem, we must
      // ensure to have a "positive definite" matrix, i.e. we are not allowed
      // to multiply all matrices by -1!
      std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";

      mechAcouFactor = CoefFunction::Generate( mp_, Global::REAL,
                                      CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT ) );
    } else {
      mechAcouFactor = CoefFunction::Generate( mp_, Global::REAL, "1.0");
    }

    PtrCoefFct regionTemp;
    std::set<UInt> definedDofs;
    if ( tempId != "" ) {
    	//just real valued temperature allowed
    	shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);
    	//Add the region information
    	PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature","name",tempId.c_str());

    	ReadUserFieldValues( eList, tempNode, tempInfo->dofNames, tempInfo->entryType,
                             false, regionTemp, definedDofs, updatedGeo_ );
    }

    //compute speed of sound
    PtrCoefFct c0;
    ComputeSOS(c0, dens, blk, regionTemp, tempId );

    //coeffc needed in PML
    shared_ptr<CoefFunction> coeffc;
    if ( sosAtLaplace_ )
    	// c0^2
    	ComputeSOS_SQR(coeffc, dens, blk, regionTemp, tempId);
    else {
    	// 1/c0^2
    	coeffc = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV) );
    }

    // check the PML node to retrieve information
    PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",id.c_str());
    // get the diagonal entries of the Jacobi matrix as vector
    shared_ptr<CoefFunction> coeffPMLVector;
    coeffPMLVector.reset(new CoefFunctionPML<Double>(pmlNode,c0,eList,regions_,true));

    // store pml factor
    matCoefs_[PML_DAMP_FACTOR]->AddRegion(eList->GetRegion(), coeffPMLVector);

    shared_ptr<CoefFunctionCompound<Double> > coefA(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > coefB(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > coefAlpha(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > coefBeta(new CoefFunctionCompound<Double>(mp_));


    // --- Set the FE ansatz for the current region ---
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",eList->GetName().c_str());
    std::string polyId = curRegNode->Get("polyId")->As<std::string>();
    std::string integId = curRegNode->Get("integId")->As<std::string>();

    shared_ptr<FeSpace> vecSpace = feFunctions_[ACOU_PMLAUXVEC]->GetFeSpace();
    vecSpace->SetRegionApproximation(eList->GetRegion(), polyId,integId);


    // ===> DEFINE PML DAMPINGFUNCTIONS
    std::map<std::string, PtrCoefFct> vars;
    std::map<std::string, PtrCoefFct> var;
    vars["a"] = coeffPMLVector;
    vars["b"] = coeffc;
    var["a"]  = coeffPMLVector;

    StdVector<std::string> matAReal;
    StdVector<std::string> matBReal;
    std::string alpha;
    std::string beta;

    if(DIM == 3 ) {
      if ( sosAtLaplace_ ) {
        const std::string Amat[] =  { "a_0_R", "0.0" , "0.0" , "0.0" , "a_1_R" , "0.0" , "0.0" , "0.0" , "a_2_R" };
        const std::string Bmat[] = {"( (a_0_R - a_1_R - a_2_R) * b_R )",
        		                    "0.0","0.0","0.0",
        		                    "( (a_1_R - a_0_R - a_2_R) * b_R )",
									"0.0","0.0","0.0",
									"( (a_2_R - a_1_R - a_0_R) * b_R )"};
        alpha = "(a_0_R + a_1_R + a_2_R )";
        beta  = "( (a_0_R * a_1_R) + (a_0_R * a_2_R) + (a_1_R * a_2_R) )";
        matAReal.Import(Amat,9);
        matBReal.Import(Bmat,9);
        coefA->SetTensor(matAReal,DIM,DIM,var);
        coefB->SetTensor(matBReal,DIM,DIM,vars);
        coefAlpha->SetScalar(alpha,var);
        coefBeta->SetScalar(beta,var);
      }
      else {
        const std::string Amat[] =  { "a_0_R", "0.0" , "0.0" , "0.0" , "a_1_R" , "0.0" , "0.0" , "0.0" , "a_2_R" };
        const std::string Bmat[] = {"( a_0_R - a_1_R - a_2_R ) ",
      		                        "0.0","0.0","0.0",
							        "( a_1_R - a_0_R - a_2_R )",
							        "0.0","0.0","0.0",
							        "( a_2_R - a_1_R - a_0_R) "};
        alpha = "( ( a_0_R + a_1_R + a_2_R) * b_R )";
        beta = " ( ( (a_0_R * a_1_R) + (a_0_R * a_2_R) + (a_1_R * a_2_R) ) * b_R )";
        matAReal.Import(Amat,9);
        matBReal.Import(Bmat,9);
        coefA->SetTensor(matAReal,DIM,DIM,var);
        coefB->SetTensor(matBReal,DIM,DIM,var);
        coefAlpha->SetScalar(alpha,vars);
        coefBeta->SetScalar(beta,vars);
      }
    }
    else {
      //dim = 2
      if ( sosAtLaplace_ ) {
        std::string Amat[] = { "a_0_R", "0.0" , "0.0" , "a_1_R" };
        std::string Bmat[] = {"( (a_0_R - a_1_R) * b_R )", "0.0" , "0.0",  "( (a_1_R - a_0_R) * b_R )"};
        alpha = "( a_0_R + a_1_R ) ";
        beta = " ( a_0_R * a_1_R ) ";
        matAReal.Import(Amat,4);
        matBReal.Import(Bmat,4);
        coefA->SetTensor(matAReal,DIM,DIM,var);
        coefB->SetTensor(matBReal,DIM,DIM,vars);
        coefAlpha->SetScalar(alpha,var);
        coefBeta->SetScalar(beta,var);
      }
      else{
        const std::string Amat[] =  { "a_0_R", "0.0" , "0.0" , "a_1_R" };
        const std::string Bmat[] = {"( a_0_R - a_1_R)", "0.0" , "0.0",  "( a_1_R - a_0_R )"};
        alpha = "( (a_0_R + a_1_R ) * b_R )";
        beta = " ( (a_0_R * a_1_R ) * b_R )";
        matAReal.Import(Amat,4);
        matBReal.Import(Bmat,4);
        coefA->SetTensor(matAReal,DIM,DIM,var);
        coefB->SetTensor(matBReal,DIM,DIM,var);
        coefAlpha->SetScalar(alpha,vars);
        coefBeta->SetScalar(beta,vars);
      }
    }

    ///now lets define some integrators
    //the vectorial auxiliary variable is called U...
    //in case of mechacoucoupling, we need to multiply all integrators which couple into acouPDE with -density or density
    //if eigenfrequency

    PtrCoefFct mAcouCorrect_CoefAlpha =
        CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp( mp_, mechAcouFactor, coefAlpha, CoefXpr::OP_MULT ) );
    PtrCoefFct mAcouCorrect_CoefBeta =
        CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp( mp_, mechAcouFactor, coefBeta, CoefXpr::OP_MULT ) );

    BaseBDBInt *     dampdPdt = new BBInt<>(new IdentityOperator<FeH1,DIM>(), mAcouCorrect_CoefAlpha, 1.0, updatedGeo_ );
    BaseBDBInt *     dampP    = new BBInt<>(new IdentityOperator<FeH1,DIM>(), mAcouCorrect_CoefBeta, 1.0, updatedGeo_ );
    //this is already integrated by parts...
    BaseBDBInt *     divU     = new ABInt<>(new GradientOperator<FeH1,DIM>(), new IdentityOperator<FeH1,DIM,DIM>(),
    		                                mechAcouFactor,1.0,updatedGeo_);

    BaseBDBInt *     dUdt     = new BBInt<>(new IdentityOperator<FeH1,DIM,DIM>(), one, 1.0, updatedGeo_ );
    BaseBDBInt *     AU       = new BDBInt<>(new IdentityOperator<FeH1,DIM,DIM>(), coefA,1.0, updatedGeo_ );
    BaseBDBInt *     BgradP   = new ADBInt<>(new IdentityOperator<FeH1,DIM,DIM>(),
    		                                 new GradientOperator<FeH1,DIM>(),coefB,1.0,updatedGeo_);

    dampdPdt->SetName("acouPML_dampdPdt");
    dampP->SetName("acouPML_dampP");
    divU->SetName("acouPML_divU");
    dUdt->SetName("acouPML_dUdt");
    AU->SetName("acouPML_AU");
    BgradP->SetName("acouPML_BgradP");

    BiLinFormContext * Context_dampdPdt   = new BiLinFormContext(dampdPdt, DAMPING );
    BiLinFormContext * Context_dampP      = new BiLinFormContext(dampP, STIFFNESS );
    BiLinFormContext * Context_divU       = new BiLinFormContext(divU, STIFFNESS );
    BiLinFormContext * Context_dUdt       = new BiLinFormContext(dUdt, DAMPING );
    BiLinFormContext * Context_AU         = new BiLinFormContext(AU, STIFFNESS );
    BiLinFormContext * Context_BgradP     = new BiLinFormContext(BgradP , STIFFNESS );

    Context_dampdPdt->SetEntities( eList, eList );
    Context_dampdPdt->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
    Context_dampP->SetEntities( eList, eList );
    Context_dampP->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
    Context_divU->SetEntities( eList, eList );
    Context_divU->SetFeFunctions(feFunctions_[formulation_],feFunctions_[ACOU_PMLAUXVEC]);
    Context_dUdt->SetEntities( eList, eList );
    Context_dUdt->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC],feFunctions_[ACOU_PMLAUXVEC]);
    Context_AU->SetEntities( eList, eList );
    Context_AU->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC],feFunctions_[ACOU_PMLAUXVEC]);
    Context_BgradP->SetEntities( eList, eList );
    Context_BgradP->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC],feFunctions_[formulation_]);


    assemble_->AddBiLinearForm( Context_dampdPdt );
    assemble_->AddBiLinearForm( Context_dampP );
    assemble_->AddBiLinearForm( Context_divU );
    assemble_->AddBiLinearForm( Context_dUdt );
    assemble_->AddBiLinearForm( Context_AU );
    assemble_->AddBiLinearForm( Context_BgradP );

    feFunctions_[ACOU_PMLAUXVEC]->AddEntityList( eList );

    if(!this->isAPML_ && DIM == 3){
      shared_ptr<FeSpace> scalSpace = feFunctions_[ACOU_PMLAUXSCALAR]->GetFeSpace();
      scalSpace->SetRegionApproximation(eList->GetRegion(), polyId,integId);
      feFunctions_[ACOU_PMLAUXSCALAR]->AddEntityList( eList);

      shared_ptr<CoefFunctionCompound<Double> > coefGamma(new CoefFunctionCompound<Double>(mp_));
      shared_ptr<CoefFunctionCompound<Double> > coefC(new CoefFunctionCompound<Double>(mp_));

      StdVector<std::string> matCReal;

      if ( sosAtLaplace_ ) {
        std::string gamma = "( a_0_R * a_1_R * a_2_R ) ";
        const std::string Cmat[]=  { " (a_1_R * a_2_R * b_R) " ,
        		                     "0.0" , "0.0" , "0.0" ,
									 " (a_0_R * a_2_R * b_R) " ,
									 "0.0", "0.0", "0.0",
									 " (a_0_R * a_1_R * b_R) "};
        matCReal.Import(Cmat,9);
        coefC->SetTensor(matCReal,3,3,vars);
        coefGamma->SetScalar(gamma,var);
      }
      else {
        std::string gamma = "( (a_0_R * a_1_R * a_2_R) * b_R )";
        const std::string Cmat[]=  { " (a_1_R * a_2_R) " , "0.0" , "0.0" , "0.0" , " (a_0_R * a_2_R) " , "0.0", "0.0", "0.0", " (a_0_R * a_1_R) "};
        matCReal.Import(Cmat,9);
        coefC->SetTensor(matCReal,3,3,var);
        coefGamma->SetScalar(gamma,vars);
      }

      PtrCoefFct mAcouCorrect_CoefGamma =
          CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp( mp_, mechAcouFactor, coefGamma, CoefXpr::OP_MULT ) );

      //the scalar auxiliary variable is called Nu...
      BaseBDBInt *     dampNu   = new BBInt<>(new IdentityOperator<FeH1,DIM>(), mAcouCorrect_CoefGamma, 1.0, updatedGeo_ );
      BaseBDBInt *     CgradNu  = new ADBInt<>(new IdentityOperator<FeH1,DIM,DIM>(),
    		                                   new GradientOperator<FeH1,DIM>(),coefC,-1.0,updatedGeo_);
      BaseBDBInt *     dNudt    = new BBInt<>(new IdentityOperator<FeH1,DIM>(), one, 1.0, updatedGeo_ );
      BaseBDBInt *     P        = new BBInt<>(new IdentityOperator<FeH1,DIM>(), one, -1.0, updatedGeo_ );

      dampNu->SetName("acouPML_dampNu");
      CgradNu->SetName("acouPML_CgIdentityOperaradNu");
      dNudt->SetName("acouPML_dNudt");
      P->SetName("acouPML_P");

      BiLinFormContext * Context_dampNu     = new BiLinFormContext(dampNu, STIFFNESS );
      BiLinFormContext * Context_CgradNu    = new BiLinFormContext(CgradNu, STIFFNESS );
      BiLinFormContext * Context_dNudt      = new BiLinFormContext(dNudt, DAMPING );
      BiLinFormContext * Context_P          = new BiLinFormContext(P, STIFFNESS );

      Context_dampNu->SetEntities( eList, eList );
      Context_dampNu->SetFeFunctions(feFunctions_[formulation_],feFunctions_[ACOU_PMLAUXSCALAR]);
      Context_CgradNu->SetEntities( eList, eList );
      Context_CgradNu->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC],feFunctions_[ACOU_PMLAUXSCALAR]);
      Context_dNudt->SetEntities( eList, eList );
      Context_dNudt->SetFeFunctions(feFunctions_[ACOU_PMLAUXSCALAR],feFunctions_[ACOU_PMLAUXSCALAR]);
      Context_P->SetEntities( eList, eList );
      Context_P->SetFeFunctions(feFunctions_[ACOU_PMLAUXSCALAR],feFunctions_[formulation_]);

      assemble_->AddBiLinearForm( Context_dampNu );
      assemble_->AddBiLinearForm( Context_CgradNu );
      assemble_->AddBiLinearForm( Context_dNudt );
      assemble_->AddBiLinearForm( Context_P );
    }
  }

  void AcousticPDE::DefineNcIntegrators() {
//	if ( complexFluidFormulation_ )
//		EXCEPTION("Complex fluid and NC-interfaces currently not allowed");

    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        if (dim_ == 2)
          DefineMortarCoupling<2,1>(formulation_, *ncIt);
        else
          DefineMortarCoupling<3,1>(formulation_, *ncIt);
        break;
      case NC_NITSCHE:
        if (dim_ == 2)
          DefineNitscheCoupling<2,1>(formulation_, *ncIt );
        else
          DefineNitscheCoupling<3,1>(formulation_, *ncIt );
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }
  
  void AcousticPDE::DefineSurfaceIntegrators( ){
    //========================================================================================
    // ABC boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    LOG_DBG(acousticpde) << "Define Surface Integrator BEGIN"<< "\n";

    if( bcNode ) {
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );
      LOG_DBG(acousticpde) << "ABCs count :  " << abcNodes.GetSize() <<  "\n" ;

      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        LOG_DBG(acousticpde) << "ABCs volRegName :  " << volRegName <<  "\n" ;

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

        //check, if region has complex fluid
        PtrParamNode curRegNode =
            myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());

        // c0 = sqrt(bulk_modulus / density)
        PtrCoefFct dens;
        PtrCoefFct blk;
        if( complexFluidFormulation_) {
          dens = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::COMPLEX );
          blk = materials_[aRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::COMPLEX );
        }
        else {
          dens = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
          blk = materials_[aRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
        }

        LOG_DBG(acousticpde) << "ABC: dens = " << dens->ToString() << "\n";
        LOG_DBG(acousticpde) << "ABC: blk  = " << blk->ToString() << "\n";

        PtrCoefFct c0;
        
        //check for temperature dependency
        curRegNode = myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());
        std::string tempId = curRegNode->Get("temperatureId")->As<std::string>();
        if (tempId != "") {
          // Get result info object for flow
          shared_ptr<ResultInfo> tempInfo = GetResultInfo(HEAT_MEAN_TEMPERATURE);

          //Add the region information
          PtrParamNode tempNode = myParam_->Get("temperatureList")->GetByVal("temperature","name",tempId.c_str());

          // Read flow coefficient function for the region and add it to flow functor
          PtrCoefFct regionTemp;
          std::set<UInt> definedDofs;
          //we assume that the temperature values are real (not complex!)
          ReadUserFieldValues( actSDList, tempNode, tempInfo->dofNames, tempInfo->entryType,
                             false, regionTemp, definedDofs, updatedGeo_ );
          // gasR=287.058 J/kg K   ... universal gas constant
          // kappa=1.402, adabatic exponent for air
          PtrCoefFct constVal = CoefFunction::Generate( mp_, Global::REAL, "402.4553160");
          c0 = CoefFunction::Generate( mp_,  Global::REAL,
                           CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT),
                           CoefXpr::OP_SQRT) );
          LOG_DBG(acousticpde) << "Define Surface Integrators standard c0 =" << c0->ToString() << "\n";
        }
        else {
        	if (complexFluidFormulation_ ) {
        		// correct factor is 1/(density*c0) = 1/sqrt(dens*compressModulus)
        		// so we modify c0 to : sqrt(density*compressModulus)
        		// division is done later
        		// note that c0 in this case does not represent the speed of sound
        		c0 = CoefFunction::Generate( mp_,  Global::COMPLEX, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens,
                                            CoefXpr::OP_MULT), CoefXpr::OP_SQRT) );
        	}
        	else
        		c0 =  CoefFunction::Generate( mp_,  Global::REAL,
                                   CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens,
                                                  CoefXpr::OP_DIV), CoefXpr::OP_SQRT) );
        }
        LOG_DBG(acousticpde) << "Def Surface Integrator:  c0 =" << c0->ToString() << "\n";


        // the following part was missing which is why abc did not function for acouPotential + mechanic
	  // if pde couples with mechanic, we have to multiply the density by -1
	  PtrCoefFct factor;
	  if ( isMechCoupled_ == true && formulation_ != ACOU_PRESSURE ) {
	    // Important: In case of a general / quadratic EV problem, we must
	    // ensure to have a "positive definite" matrix, i.e. we are not allowed
	    // to multiply all matrices by -1!
	    std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";

	    factor = CoefFunction::Generate( mp_, Global::REAL,
	 					    CoefXprBinOp(mp_, dens, stringFac, CoefXpr::OP_MULT ) );
	  } else {
	    factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");
	  }
	  LOG_DBG(acousticpde) << "Def Surface Integrator: factor =" << factor->ToString() << "\n";

        PtrCoefFct coeffDamp;
        if ( sosAtLaplace_ ) {
          // factor for damping matrix: factor * c0
          coeffDamp = CoefFunction::Generate( mp_, Global::REAL,
                         			CoefXprBinOp(mp_, factor, c0, CoefXpr::OP_MULT ) );
        }
        else {
          // factor for damping matrix: factor / c0
         if (complexFluidFormulation_ )
           coeffDamp = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, factor, c0, CoefXpr::OP_DIV ) );
         else
           coeffDamp = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, factor, c0, CoefXpr::OP_DIV ) );

        }
        LOG_DBG(acousticpde) << "Define Surface Integrator: coeffDamp =" << coeffDamp->ToString() << "\n";
        
        BiLinearForm * abcInt = NULL;
        if (complexFluidFormulation_) {
          if( dim_ == 2 )
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1, Complex>(), coeffDamp, 1.0, updatedGeo_ );
          else
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1, Complex>(), coeffDamp, 1.0, updatedGeo_ );
        }

        else {
          if( dim_ == 2 )
            abcInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), coeffDamp, 1.0, updatedGeo_ );
          else
        	abcInt = new BBInt<>(new IdentityOperator<FeH1,3,1>(), coeffDamp, 1.0, updatedGeo_ );

        }

        FEMatrixType targetMatrix = DAMPING;
        if(updatedGeo_){
          targetMatrix = DAMPING_UPDATE;
        }

        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, targetMatrix );

        abcContext->SetEntities( actSDList, actSDList );
        abcContext->SetFeFunctions( feFunctions_[formulation_] , feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }

      //========================================================================================
      // Impedance boundaries
      // TODO: implement impedance BC
      //========================================================================================
      ParamNodeList impedNodes = bcNode->GetList( "impedance" );

      for( UInt i = 0; i < impedNodes.GetSize(); ++i ) {
    	if ( complexFluidFormulation_ )
    		EXCEPTION("Complex fluid and impedance currently not working");

        BiLinearForm * impedInt = NULL;
        std::string regionName = impedNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, regionName );
        std::string volRegName = impedNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        bool isNormalised = impedNodes[i]->Get("isNormalised")->As<bool>();
        std::string impId = impedNodes[i]->Get("impedanceId")->As<std::string>();

        shared_ptr<CoefFunctionImpedanceModel<Complex> >
                   Z_impMod(new CoefFunctionImpedanceModel<Complex>(mp_, materials_[aRegion], isNormalised));

        PtrParamNode impListNodes = myParam_->Get( "impedanceList", ParamNode::PASS );
        if ( !impListNodes ) {
          EXCEPTION("No impedance list in *.xml");
        }
        ParamNodeList impChildren = impListNodes->GetChildren();
        UInt i_imp = 0;
        for ( ; i_imp < impChildren.GetSize(); ++i_imp) {
          std::string impType = impChildren[i_imp]->GetName();
          std::string actId = impChildren[i_imp]->Get("id")->As<std::string>();

          if (actId != impId) {
            continue;
          }
          if (impType == "mpp") {
            std::string impSubType = impChildren[i_imp]->Get("subtype")->As<std::string>();
            if (impSubType == "slit") {
              Z_impMod->GenerateSlitMpp(impChildren[i_imp]);
            } else if (impSubType == "circ") {
              Z_impMod->GenerateCircMpp(impChildren[i_imp]);
            } else {
              EXCEPTION("No such impedance type: " << impType << " with subtype: " << impSubType);
            }
          } else if (impType == "interpolImpedance"){
            Z_impMod->GenerateInterpolImpedance(impChildren[i_imp]);
          } else if (impType == "fctImpedance") {
            Z_impMod->GenerateImpedanceFct(impChildren[i_imp]);
          } else {
            EXCEPTION("No such impedance type: " << impType);
          }
          break;

        }
        if (i_imp == impChildren.GetSize()) {
          EXCEPTION("No such impedance id found: " << impId);
        }

        if( dim_ == 2 ) {
          impedInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1, Complex>(), Z_impMod, 1.0, updatedGeo_ );
        } else {
          impedInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1, Complex>(), Z_impMod, 1.0, updatedGeo_ );
        }

        FEMatrixType targetMatrix = DAMPING;
        if(updatedGeo_){
          targetMatrix = DAMPING_UPDATE;
        }

        impedInt->SetName("impedIntegrator");
        BiLinFormContext *impedContext = new BiLinFormContext(impedInt, targetMatrix );

        impedContext->SetEntities( actSDList, actSDList );
        impedContext->SetFeFunctions( feFunctions_[formulation_] , feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( impedContext );
      }

      //========================================================================================
      // boundary Layers
      //========================================================================================
      ParamNodeList blNodes = bcNode->GetList( "boundaryLayer" );

      for( UInt i = 0; i < blNodes.GetSize(); ++i ) {
        if ( !(formulation_ == ACOU_PRESSURE) || !(this->analysistype_ == HARMONIC) ) EXCEPTION("bounadryLayer only available for harmonic pressure formulation");

        std::string regionName = blNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST, regionName );
        std::string volRegName = blNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType volRegion = ptGrid_->GetRegion().Parse(volRegName);

        PtrCoefFct rho0 = materials_[volRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
        PtrCoefFct K = materials_[volRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
        PtrCoefFct cp = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("cp")->As<std::string>() );
        PtrCoefFct cv = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("cv")->As<std::string>() );
        PtrCoefFct nu = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("nu")->As<std::string>() );
        PtrCoefFct k = CoefFunction::Generate( mp_,  Global::REAL, blNodes[i]->Get("k")->As<std::string>() );

        LOG_DBG(acousticpde) << "Define Surface Integrator: regionName =" << regionName << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: volRegName =" << volRegName << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: rho0       =" << rho0->ToString() << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: K          =" << K->ToString() << "\n";
        LOG_DBG(acousticpde) << "Define Surface Integrator: k          =" << k->ToString() << "\n";

        PtrCoefFct omegaHalv = CoefFunction::Generate( mp_,  Global::REAL, "pi*f");//
        // deltaV = sqrt( 2*nu/omega )
        PtrCoefFct deltaV = CoefFunction::Generate(mp_,Global::REAL, CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, nu, omegaHalv, CoefXpr::OP_DIV ) , CoefXpr::OP_SQRT));
        // deltaT = sqrt( 2*k/(omega*rho0*cp) )
        PtrCoefFct deltaT = CoefFunction::Generate(mp_,Global::COMPLEX,
            CoefXprUnaryOp(mp_,
                CoefXprBinOp(mp_,
                    k, CoefXprBinOp(mp_,CoefXprBinOp(mp_, cp, rho0, CoefXpr::OP_MULT ),omegaHalv, CoefXpr::OP_MULT ),
                    CoefXpr::OP_DIV ),
                    CoefXpr::OP_SQRT));
        // 1-i ... common factor for both terms
        PtrCoefFct oneMinusI =  CoefFunction::Generate( mp_,  Global::COMPLEX, "1", "-1");
        // gamma = cp/cv -> gamma-1 = (cp-cv)/cv
        PtrCoefFct gammaMinusOne = CoefFunction::Generate(mp_,Global::COMPLEX, CoefXprBinOp(mp_, CoefXprBinOp(mp_,cp,cv,CoefXpr::OP_SUB), cv, CoefXpr::OP_DIV));

        // Mass matrix
        PtrCoefFct coefM =  CoefFunction::Generate(mp_,Global::COMPLEX,
            CoefXprBinOp(mp_,
                CoefXprBinOp(mp_,
                    CoefXprBinOp(mp_, deltaT, CoefXprBinOp(mp_, rho0, K, CoefXpr::OP_DIV), CoefXpr::OP_MULT),
                    oneMinusI,
                    CoefXpr::OP_MULT),
                    gammaMinusOne,
                    CoefXpr::OP_MULT)
        );
        BiLinearForm * blmInt = NULL;
        if( dim_ == 2 ) {
          blmInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1, Complex>(), coefM, 0.5, updatedGeo_ );
        } else {
          blmInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1, Complex>(), coefM, 0.5, updatedGeo_ );
        }
        blmInt->SetName("blMassIntegrator");
        BiLinFormContext *blmContext = new BiLinFormContext(blmInt, MASS);
        blmContext->SetEntities(actSDList, actSDList);
        blmContext->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(blmContext);

        LOG_DBG(acousticpde) << "Define Surface Integrator boundary layer: coefM =" << coefM->ToString() << "\n";

        // Stiffness matrix
        PtrCoefFct coefK = CoefFunction::Generate(mp_,Global::COMPLEX, CoefXprBinOp(mp_, deltaV, oneMinusI, CoefXpr::OP_MULT) );
        BiLinearForm * blkInt = NULL;
        std::set<RegionIdType> volRegions;
        volRegions.insert(volRegion);
        if( dim_ == 2 ) {
          blkInt = new BBInt<Complex>(new GradientOperator<FeH1,2,1, Complex>(), coefK, -0.5, updatedGeo_ );
        } else {
          blkInt = new BBInt<Complex>(new GradientOperator<FeH1,3,1, Complex>(), coefK, -0.5, updatedGeo_ );
        }
        blkInt->SetName("blStiffnessIntegrator");
        BiLinFormContext *blkContext = new BiLinFormContext(blkInt, STIFFNESS);
        blkContext->SetEntities(actSDList, actSDList);
        blkContext->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList(actSDList);
        assemble_->AddBiLinearForm(blkContext);

        LOG_DBG(acousticpde) << "Define Surface Integrator boundary layer: coefK =" << coefK->ToString() << "\n";
      } // boundary Layers
    } // end if ( bcNode )
    LOG_DBG(acousticpde) << "Define Surface Integrator END"<< "\n";
  } // DefineSurfaceIntegrators

  void AcousticPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(acousticpde) << "Defining rhs load integrators for acoustic PDE";
    
    // Get FESpace and FeFunction
    shared_ptr<BaseFeFunction> myFct = feFunctions_[formulation_];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> empty;
    LinearForm * lin = NULL;

    // obtain density
    shared_ptr<CoefFunctionMulti> densFct  = matCoefs_[ELEM_DENSITY];
    shared_ptr<CoefFunctionSurf> surfDens(new CoefFunctionSurf(false));
    surfDens->SetVolumeCoefs( densFct->GetRegionCoefs() );
    

    // In the case of acou-mech coupling we have to multiply the
    // integrators by -density
    Double scalFactor = 1.0;
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

    bool coefUpdateGeo;
    // ===========================
    //  NORMAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading normal velocity"; 
    
    ReadRhsExcitation( "normalVelocity", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("NormalVelocity-BC and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }
      PtrCoefFct exValue;
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
        scalFactor = -1.0;
        exValue = 
            CoefFunction::Generate( mp_, part,
                                   CoefXprBinOp(mp_, coef[i],surfDens, CoefXpr::OP_MULT) );
      } else {
        exValue = coef[i];
      }
      
      if( formulation_ == ACOU_POTENTIAL ) {
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,2,1>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
                                                  1.0, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,3,1>(),
                                                 1.0, exValue , volRegions, coefUpdateGeo);
          }
        }
      } else if( formulation_ == ACOU_PRESSURE && isComplex_) {

        // in this case the pressure can be related to the 
        // normal velocity as 
        // p_n = - j*Omega*v_n*rho
    	PtrCoefFct tmp = CoefFunction::Generate( mp_, part,"0.0",
    	                                                "-2*pi*f");
    	if ( complexFluidFormulation_ ) {
    		//do not need multiplication with density
    		exValue = CoefFunction::Generate( mp_, part,
                    CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
    	}
    	else {
    		PtrCoefFct tmp2 = CoefFunction::Generate( mp_, part,
                                 CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
    		exValue = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, tmp2,surfDens,
                                              CoefXpr::OP_MULT) );
    	}

        if( dim_ == 2) {
          lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
                                                1.0, exValue, volRegions, coefUpdateGeo);
        } else  {
          lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
                                                1.0, exValue, volRegions, coefUpdateGeo);
        }

      } else {
        EXCEPTION( "Normal velocity can only be prescribed for potential "
            << "formulation or for harmonic pressure formulation" )
      }
      lin->SetName("NormalVelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

    
    // ===========================
    //  TOTAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading total velocity";
    StdVector<std::string> vecDofNames;
    if(dim_ == 3)
      vecDofNames = "x", "y", "z";
    if(dim_ == 2 && !isaxi_)
      vecDofNames = "x", "y";
    if(dim_ == 2 && isaxi_)
      vecDofNames = "r", "z";
    ReadRhsExcitation( "velocity", vecDofNames, ResultInfo::VECTOR, isComplex_,
                       ent, coef,coefUpdateGeo );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }

      PtrCoefFct exValue;
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
        scalFactor = -1.0;
        exValue = 
            CoefFunction::Generate( mp_, part,
                                   CoefXprBinOp(mp_, coef[i],surfDens, CoefXpr::OP_MULT) );
      } else {
        exValue = coef[i];
      }
      if( formulation_ == ACOU_POTENTIAL ) {
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,2>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(),
                                                  scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,3>(),
                                                 scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        }
      } else if( formulation_ == ACOU_PRESSURE && isComplex_) {

        // in this case the pressure can be related to the 
        // normal velocity as 
        // p_n = - j*Omega*v_n*rho
        PtrCoefFct tmp = CoefFunction::Generate( mp_, part,"0.0",
                                                "-2*pi*f");
        if ( complexFluidFormulation_ ) {
        	//do not need multiplication with density
        	exValue = CoefFunction::Generate( mp_, part,
        			CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
        }
        else {
          PtrCoefFct tmp2 = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
          exValue = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, surfDens, tmp2, CoefXpr::OP_MULT) );
        }

        if( dim_ == 2) {
          lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(),
                                                scalFactor, exValue, volRegions, coefUpdateGeo);
        } else  {
          lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(),
                                                scalFactor, exValue, volRegions, coefUpdateGeo);
        }
      } else {
        EXCEPTION( "Normal velocity can only be prescribed for potential "
            << "formulation or for harmonic pressure formulation" )
      }
      lin->SetName("VelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }
    

    // ===========================
    //  NORMAL ACCELERATION (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading normal acceleration";

    ReadRhsExcitation( "normalAcceleration", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal acceleration can only be defined on surface elements");
      }
      scalFactor = 1.0;
      PtrCoefFct exValue;
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
    	  EXCEPTION( "Normal acceleration can only be prescribed for pressure"
    	              << "formulation" );
      } else {
        exValue = coef[i];
      }

      if( formulation_ == ACOU_PRESSURE ) {
    	  exValue =
    			  CoefFunction::Generate( mp_, part,
    					  CoefXprBinOp(mp_, coef[i],surfDens, CoefXpr::OP_MULT) );

        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,2,1>(),
            		scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,2,1>(),
            		scalFactor, exValue, volRegions, coefUpdateGeo);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex,true>( new IdentityOperator<FeH1,3,1>(),
            			scalFactor, exValue, volRegions, coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double,true>( new IdentityOperator<FeH1,3,1>(),
            		scalFactor, exValue , volRegions, coefUpdateGeo);
          }
        }
      } else if( formulation_ == ACOU_POTENTIAL) {
        EXCEPTION( "Normal acceleration can only be prescribed for pressure "
            << "formulation" )
      }
      lin->SetName("NormalAccelerationIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }


    // ===========================
    //  PRESSURE (surface) 
    // ===========================
    // Only possible for harmonic potential formulation
    ReadRhsExcitation( "pressure", empty, ResultInfo::SCALAR, isComplex_,
                       ent, coef, coefUpdateGeo );
    if( formulation_ == ACOU_POTENTIAL && isComplex_ ) { 
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
         if ( sosAtLaplace_)
          EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

        // ensure that list contains only surface elements
        EntityIterator it = ent[i]->GetIterator();
        UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
        if( elemDim != (dim_-1) ) {
          EXCEPTION("Pressure can only be defined on surface elements");
        }
        PtrCoefFct exValue;
        if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
          scalFactor = -1.0;
          exValue = 
              CoefFunction::Generate( mp_, part,
                                     CoefXprBinOp(mp_, coef[i],surfDens, CoefXpr::OP_MULT) );
        } else {
          exValue = coef[i];
        }
        // psi_n = j * 1 / (omega*rho) * p
        PtrCoefFct tmp = 
             CoefFunction::Generate( mp_, Global::COMPLEX, "0","1/(2*pi*f)");
         PtrCoefFct tmp2 = 
             CoefFunction::Generate( mp_, Global::COMPLEX,
                                    CoefXprBinOp(mp_, tmp, surfDens, CoefXpr::OP_DIV ) );
         exValue = CoefFunction::Generate( mp_, part, 
                                          CoefXprBinOp(mp_, exValue, tmp2, 
                                                       CoefXpr::OP_MULT) );
         shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
         actBc->entities = ent[i];
         actBc->result = myFct->GetResultInfo();
         actBc->value = exValue;
         actBc->dofs.insert(0);
         actBc->updatedGeo = coefUpdateGeo;
         myFct->AddInhomDirichletBc(actBc);
      } // loop entities
    } // if

    // =====================
    //  ACOUSTIC SOURCE DENSITY
    // =====================
    //LOG_DBG(heatcondpde) << "Reading acoustic source density";

    ReadRhsExcitation( "acouSourceDensity", empty,
                       ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Acoustic source density must be defined on elements")
      }
      EntityIterator it = ent[i]->GetIterator();
      it.Begin();

      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(),
                                          Complex(1.0), coef[i], coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1>(),
                                         1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("AcousticSourceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);

      RegionIdType regionId = ent[i]->GetRegion();
      acousticSourceDensityCoef_->AddRegion(  regionId, coef[i] );
    } // for

    // =====================================
    //  rhsValues for e.g. for aeroacoustics
    // =====================================
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if ( sosAtLaplace_)
        EXCEPTION("RHS and speed of sound at Laplace operator currently not possible");

      coef[i]->SetConservative(true);
      //std::cout << "ADD" << std::endl;
      coef[i]->SetFeFunction( feFunctions_[formulation_], formulation_);

      if ( !approxSourceWithDeltaFnc_ && coef[i]->GetInverseType() == CoefFunction::INVSOURCE )  {
    	  std::cout << "Add BUIntegrator SRC" << std::endl;
    	  // check type of entitylist
          if (ent[i]->GetType() == EntityList::NODE_LIST) {
            EXCEPTION("Acoustic source density must be defined on elements")
          }
          EntityIterator it = ent[i]->GetIterator();
          it.Begin();

          if(isComplex_) {
            lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(),
                                              Complex(1.0), coef[i], coefUpdateGeo);
          } else  {
            lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1>(),
                                             1.0, coef[i], coefUpdateGeo);
          }
          lin->SetName("AcousticSourceDensityInverseInt");
          LinearFormContext *ctx = new LinearFormContext( lin );
          ctx->SetEntities( ent[i] );
          ctx->SetFeFunction(myFct);
          assemble_->AddLinearForm(ctx);

          RegionIdType regionId = ent[i]->GetRegion();
          acousticSourceDensityCoef_->AddRegion(  regionId, coef[i] );
      }

      //if ( coef[i]->GetInverseType() == CoefFunction::INVSOURCE )
      this->rhsFeFunctions_[formulation_]->AddLoadCoefFunction(coef[i], ent[i]);
    }
    
       // ================
       //  RHS DENSITY Vector
       //	For reading vectoriel source terms of the lightill analogy (div(T)) 
       //      Diss Freidhager 2022
       // ================
       LOG_DBG(acousticpde) << "Reading rhsDensityVector";
       ReadRhsExcitation( "rhsDensityVector", vecDofNames, ResultInfo::VECTOR,
    		   isComplex_, ent, coef, updatedGeo_ );

       //perform checks



       for( UInt i = 0; i < ent.GetSize(); ++i ) {
    	   	   	   	   	   	   	   	   std::cout << "ausgabe     " << i << std::endl;
         // check type of entitylist
         //if (ent[i]->GetType() == EntityList::NODE_LIST) {
         //  EXCEPTION("Rhs density Vector must be defined on elements, because BUIntegrator requires a density")
         //}
         EntityIterator it = ent[i]->GetIterator();
         it.Begin();

         if( formulation_ == ACOU_POTENTIAL){
      	   EXCEPTION("Using rhsDensityVector based on Lighthill's equation -> AcouPressure Formulation!")
         } //closes if AcouPotential
         else{
        	 if(isComplex_) {
           		 std::cout << "\t \t \t Using complexe vector formulation of Lighthill div(T) in grad(phi)" << std::endl;
                 scalFactor = 1.0;
        		 if( dim_ == 2) {

        		 	 lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex()>(),Complex(scalFactor), coef[i], updatedGeo_);
        		 }else{

            		 lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,3,1,Complex()>(),Complex(scalFactor), coef[i], updatedGeo_);
        		 }

				} else  { //Transiente Simulation

				std::cout << "\t \t \t Using vector formulation of Lighthill div(T) in grad(phi)" << std::endl;
				scalFactor = 1.0;
				if( dim_ == 2) {
					lin = new BUIntegrator<Double,false>( new GradientOperator<FeH1,2,1>(),scalFactor, coef[i], updatedGeo_);
				}else{
					lin = new BUIntegrator<Double,false>( new GradientOperator<FeH1,3,1>(),scalFactor, coef[i], updatedGeo_);
				}
				}
         }

         lin->SetName("RhsDensityVectorInt");
         LinearFormContext *ctx = new LinearFormContext( lin );
         ctx->SetEntities( ent[i] );
         ctx->SetFeFunction(myFct);
         assemble_->AddLinearForm(ctx);
     }


       // ================
       //  BC Vector
       //  Requires a vector v for which v in n is computed on BC surfaces.
       //  Therefore the vector v has to be interpolated onto BC surfaces first
       //  v has to be on elements, since it is a numeric density
       //  For Lighthill thsi means computing div(T)+grad(c_0^2 rho) as BC Vector
       //  Diss Freidhager 2022
       // ================
       LOG_DBG(acousticpde) << "Reading BCVector";
             ReadRhsExcitation( "BCVector", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );

             for( UInt i = 0; i < ent.GetSize(); ++i ) {
               EntityIterator it = ent[i]->GetIterator();
               it.Begin();
               // check type of entitylist
		// if (ent[i]->GetType() == EntityList::NODE_LIST) {
		// EXCEPTION("BCVector must be defined on elements, because BUIntegrator requires a density")
		// }

               // ensure that list contains only surface elements
               UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
               if( elemDim != (dim_-1) ) {
                 EXCEPTION("BCVector can only be defined on surface elements");
               }


               if ( sosAtLaplace_) {
            	   EXCEPTION("BCVector not defined for sosAtLaplace");
             	}//closes sosAtLaplace
              else{
            	  if(isComplex_) {
            		  std::cout << "\t \t \t Surface Term of vector formulation of Lighthill's analogy" << std::endl;
            		  scalFactor = 1.0;
		          if( dim_ == 2) {
			      lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,2>(),
								Complex(scalFactor), coef[i],volRegions, coefUpdateGeo); //IdentityOperatorNormal

			   } else  {
			      lin = new BUIntegrator<Complex,true>( new IdentityOperatorNormal<FeH1,3>(),
								Complex(scalFactor), coef[i],volRegions, coefUpdateGeo);
			   }
            	  }
            	  else {

            	          std::cout << "\t \t \t Surface Term of vector formulation of Lighthill's analogy" << std::endl;
            	          scalFactor = 1.0;
                         if( dim_ == 2) {
            	         	lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,2>(),
            	         			scalFactor, coef[i],volRegions, coefUpdateGeo);
                         } else  {
                              lin = new BUIntegrator<Double,true>( new IdentityOperatorNormal<FeH1,3>(),
            	         			scalFactor, coef[i],volRegions, coefUpdateGeo);
            	         }
            	  }
               }

               lin->SetName("BCVector");
               LinearFormContext *ctx = new LinearFormContext( lin );
               ctx->SetEntities( ent[i] );
               ctx->SetFeFunction(myFct);
               assemble_->AddLinearForm(ctx);
               myFct->AddEntityList(ent[i]);

//               RegionIdType regionId = ent[i]->GetRegion();
//               BCSurfaceVectorCoef_->AddRegion(  regionId, coef[i] );

             }
                      

    // ================
    //  RHS DENSITY
    // ================
    ReadRhsExcitation( "rhsDensity", empty,
                       ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Rhs density must be defined on elements")
      }
      if ( sosAtLaplace_) {
    	//get region id
    	std::string volRegName = ent[i]->GetName();
    	RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

    	//get adiabatic exponent
    	PtrCoefFct density = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );

    	//get speed of sound and square it
    	std::map<RegionIdType,PtrCoefFct > c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
    	PtrCoefFct regionC0 = c0Fcts[aRegion];
    	PtrCoefFct regionC02 = CoefFunction::Generate( mp_, Global::REAL,
    			CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT ) );

    	PtrCoefFct tmpCoef = CoefFunction::Generate( mp_, Global::REAL,
              		  CoefXprBinOp(mp_, regionC02, density, CoefXpr::OP_MULT ) );

    	if(isComplex_) {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::COMPLEX,
    	      	    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
    				  Complex(scalFactor), coefRHS, updatedGeo_);
    	}
    	else {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::REAL,
    		         		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
    			                          scalFactor, coefRHS, updatedGeo_);
    	}
      }
      else {
        if(isComplex_) {
        	  lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
        				  Complex(scalFactor), coef[i], updatedGeo_);
        }
        else {
        	lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
        			                          scalFactor, coef[i], updatedGeo_);
        }
      }
      lin->SetName("RhsDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for


    // ================
    //  RHS Mass time derivative
    // ================
    ReadRhsExcitation( "rhsMassTimeDeriv", empty,
                       ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("rhsMassTimeDeriv must be defined on elements");
      }
      if ( !sosAtLaplace_ )
    	  EXCEPTION("rhsMassTimeDeriv currently just available for speed of sound at Laplace operator");

      //source term is partial time derivative of speed of sound squared
      //multiplied by mass, which has been evaluated by "ReadRhsExcitation"
      //and stored in coef[i]
//      if ( this->analysistype_ == HARMONIC )
//    	  EXCEPTION("rhsMassTimeDeriv in harmonic case for speed of sound at Laplace operator\
//    		         currently not available");

      //get region id
      std::string volRegName = ent[i]->GetName();
      RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

      //get adiabatic exponent
      PtrCoefFct aExp = materials_[aRegion]->GetScalCoefFnc( FLUID_ADIABATIC_EXPONENT, Global::REAL );

      //get speed of sound and square it
      std::map<RegionIdType,PtrCoefFct > c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
      PtrCoefFct regionC0 = c0Fcts[aRegion];
      PtrCoefFct regionC02 = CoefFunction::Generate( mp_, Global::REAL,
    		  CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT ) );

      PtrCoefFct tmpCoef = CoefFunction::Generate( mp_, Global::REAL,
    		  CoefXprBinOp(mp_, regionC02, aExp, CoefXpr::OP_DIV ) );

      if (isComplex_) {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::COMPLEX,
    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );

    	  lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
        		                          scalFactor, coefRHS, updatedGeo_);
      }
      else {
          PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::REAL,
        		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );

          lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
            		                          scalFactor, coefRHS, updatedGeo_);
      }
      lin->SetName("RhsMassTimeDerivInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for


    // ================
    //  RHS Mass convective term
    // ================
    ReadRhsExcitation( "rhsMassConvective", empty,
                       ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("rhsMassConvective must be defined on elements");
      }

      //source term is convective derivative of speed of sound squared
      //multiplied by mass, which has been evaluated by "ReadRhsExcitation"
      //and stored in coef[i]

      //get region id
      std::string volRegName = ent[i]->GetName();
      RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);

      //get adiabatic exponent
      PtrCoefFct aExp = materials_[aRegion]->GetScalCoefFnc( FLUID_ADIABATIC_EXPONENT, Global::REAL );

      //get speed of sound and square it
      std::map<RegionIdType,PtrCoefFct > c0Fcts = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->GetRegionCoefs();
      PtrCoefFct regionC0 = c0Fcts[aRegion];
      PtrCoefFct regionC02 = CoefFunction::Generate( mp_, Global::REAL,
    		  CoefXprBinOp(mp_, regionC0, regionC0, CoefXpr::OP_MULT ) );

      PtrCoefFct tmpCoef = CoefFunction::Generate( mp_, Global::REAL,
          		  CoefXprBinOp(mp_, regionC02, aExp, CoefXpr::OP_DIV ) );

      BaseBOperator* bOp;
      if ( dim_ == 2)
    	  bOp = new ConvectiveOperator<FeH1,2,1>();
      else
    	  bOp = new ConvectiveOperator<FeH1,3,1>();
      bOp->SetCoefFunction( meanFlowCoef_ ); //meanFlow )

      if (isComplex_) {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::COMPLEX,
    	    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Complex>( bOp, scalFactor, coefRHS, updatedGeo_);
      }
      else {
    	  PtrCoefFct coefRHS = CoefFunction::Generate( mp_, Global::REAL,
    	    		  CoefXprBinOp(mp_, tmpCoef, coef[i], CoefXpr::OP_MULT ) );
    	  lin = new BUIntegrator<Double>( bOp, scalFactor, coefRHS, updatedGeo_);
      }

      lin->SetName("RhsMassConvectiveInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for

  }

  void AcousticPDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticPDE::DefinePrimaryResults(){
	//check for complex fluid formulation
	RegionIdType actRegion;
	std::map<RegionIdType, BaseMaterial*>::iterator it;
	for ( it = materials_.begin(); it != materials_.end(); it++ ) {
		actRegion = it->first;
		std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
		PtrParamNode curRegNode =
		   			myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
		if ( curRegNode->Get("complexFluid")->As<std::string>() == "yes" ) {
			complexFluidFormulation_ = true;
			isMaterialComplex_ = true;
			if ( this->analysistype_ != HARMONIC )
				EXCEPTION("Complex fluid region just allowed in harmonic analysis");
	   		//need an acoustic pressure formulation
			if ( formulation_ != ACOU_PRESSURE )
				EXCEPTION("Complex fluid needs acoustic pressure formulation");
		}
	}

    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    if ( formulation_ ==  ACOU_PRESSURE) {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "";
      res1->unit = MapSolTypeToUnit(ACOU_PRESSURE);
    } else {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "";
      res1->unit = MapSolTypeToUnit(ACOU_POTENTIAL);
    }
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[formulation_]->SetResultInfo(res1);
    results_.Push_back( res1 );
    res1->SetFeFunction(feFunctions_[formulation_]);
    DefineFieldResult( feFunctions_[formulation_], res1 );
    

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    if ( formulation_ ==  ACOU_PRESSURE) {
      hdbcSolNameMap_[ACOU_PRESSURE] = "soundSoft";
      idbcSolNameMap_[ACOU_PRESSURE] = "pressure";
    } else {
      hdbcSolNameMap_[ACOU_POTENTIAL] = "soundSoft";
      idbcSolNameMap_[ACOU_POTENTIAL] = "potential";
    }
    
    // === ACOUSTIC RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ACOU_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    this->rhsFeFunctions_[formulation_]->SetResultInfo(rhs);
    DefineFieldResult( this->rhsFeFunctions_[formulation_], rhs );
    results_.Push_back( rhs );
    availResults_.insert( rhs );

    //creates the mean flow
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

    //create result info for mean flow velocity
    CreateMeanFlowFunction(vecDofNames);
    
    //create result info for mean temperature field
    //(mean in the sense of time averaged)
    shared_ptr<ResultInfo> temperature( new ResultInfo);
    temperature->resultType = HEAT_MEAN_TEMPERATURE;
    temperature->dofNames = "";
    temperature->unit = MapSolTypeToUnit(HEAT_MEAN_TEMPERATURE);

    temperature->definedOn = ResultInfo::NODE;
    temperature->entryType = ResultInfo::SCALAR;

   	meanTemperatureCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR,1,1,false));
    DefineFieldResult( meanTemperatureCoef_, temperature );

    // ==============================================
    //  Define CoefFunctions for material parameters
    // ==============================================
    
    // === DENSITY ===
    // Density \rho_{\mathrm{a}} = {\frac{p_{\mathrm{a}}} {c^{2}}}
    shared_ptr<ResultInfo> density ( new ResultInfo );
    density->resultType = ELEM_DENSITY;
    density->dofNames = "";
    density->unit = MapSolTypeToUnit(ELEM_DENSITY);
    density->definedOn = ResultInfo::ELEMENT;
    density->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> densFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,
    		                              complexFluidFormulation_));
    matCoefs_[ELEM_DENSITY] = densFct;
    DefineFieldResult(densFct, density);
    
    // === SPEED OF SOUND ===
    // Speed of sound c = \sqrt{\kappa \frac{p_0} {\rho_0}} = \sqrt{\kappa R T_{0}}
    shared_ptr<ResultInfo> sos ( new ResultInfo );
    sos->resultType = ACOU_ELEM_SPEED_OF_SOUND;
    sos->dofNames = "";
    sos->unit = MapSolTypeToUnit(ACOU_ELEM_SPEED_OF_SOUND);
    sos->definedOn = ResultInfo::ELEMENT;
    sos->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> sosFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,
    		                             complexFluidFormulation_));
    matCoefs_[ACOU_ELEM_SPEED_OF_SOUND] = sosFct;
    DefineFieldResult(sosFct, sos);
  }
  
  void AcousticPDE::FinalizePostProcResults(){
    //first call base class method
    SinglePDE::FinalizePostProcResults();
  }

  void AcousticPDE::DefinePostProcResults(){
    
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[formulation_];
    shared_ptr<ResultInfo> res1 = feFct->GetResultInfo();
    
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

    StdVector<std::string> tensorDofNames;
    if( ptGrid_->GetDim() == 3 ) {
      tensorDofNames.Resize(9);
      tensorDofNames[0] = "xx";
      tensorDofNames[1] = "xy";
      tensorDofNames[2] = "xz";
      tensorDofNames[3] = "yx";
      tensorDofNames[4] = "yy";
      tensorDofNames[5] = "yz";
      tensorDofNames[6] = "zx";
      tensorDofNames[7] = "zy";
      tensorDofNames[8] = "zz";
    } else {
      if( ptGrid_->IsAxi() ) {
        // Tensor DOF names in are not yet used in axisymmetric problems.
        tensorDofNames.Resize(4);
        tensorDofNames[0] = "rr";
        tensorDofNames[1] = "rz";
        tensorDofNames[2] = "zr";
        tensorDofNames[3] = "zz";
      } else { // 2D case
        tensorDofNames.Resize(4);
        tensorDofNames[0] = "xx";
        tensorDofNames[1] = "xy";
        tensorDofNames[2] = "yx";
        tensorDofNames[3] = "yy";
      }
    }

    // === PRESSURE / POTENTIAL - 1.DERIVATIVE ===
    // first time derivative of potential 
    // \frac{\partial \psi_mathrm{a}} {\partial t} or j \omega \psi_mathrm{a}
    // or of pressure 
    // \frac{\partial p_{\mathrm{a}}} {\partial t} or j \omega p_{\mathrm{a}}
    shared_ptr<ResultInfo> deriv1(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv1->resultType = ACOU_POTENTIAL_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = MapSolTypeToUnit(ACOU_POTENTIAL_DERIV_1);
    } else {
      deriv1->resultType = ACOU_PRESSURE_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = MapSolTypeToUnit(ACOU_PRESSURE_DERIV_1);
    }
    deriv1->entryType = res1->entryType;
    deriv1->definedOn = res1->definedOn;
    availResults_.insert( deriv1 );
    DefineTimeDerivResult( deriv1->resultType, 1, formulation_ );

    // === PRESSURE / POTENTIAL - 2.DERIVATIVE ===
    // second time derivative of potential 
    // \frac{\partial^{2} \psi_mathrm{a}} {\partial t^{2}} \mathrm or -\omega^{2} \psi_mathrm{a}
    // or of pressure 
    // \frac{\partial^{2} p_{\mathrm{a}}} {\partial t^{2}} or -\omega^{2} p_{\mathrm{a}}
    shared_ptr<ResultInfo> deriv2(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv2->resultType = ACOU_POTENTIAL_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = MapSolTypeToUnit(ACOU_POTENTIAL_DERIV_2);
    } else {
      deriv2->resultType = ACOU_PRESSURE_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = MapSolTypeToUnit(ACOU_PRESSURE_DERIV_2);
    }
    deriv2->entryType = res1->entryType;
    deriv2->definedOn = res1->definedOn;
    availResults_.insert( deriv2 );
    DefineTimeDerivResult( deriv2->resultType, 2, formulation_ );

    
    // define acoustic pressure, depending on formulation
    PtrCoefFct presFct;
    if( formulation_ == ACOU_POTENTIAL ) {
      shared_ptr<ResultInfo> pres(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = MapSolTypeToUnit(ACOU_PRESSURE);
      pres->entryType = ResultInfo::SCALAR;
      pres->definedOn = ResultInfo::ELEMENT;
      // Define pressure as p = rho * dPsi/dt
      PtrCoefFct potD1Fct= this->GetCoefFct( ACOU_POTENTIAL_DERIV_1 );
      PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
      presFct  =
               CoefFunction::Generate( mp_, part,
                                      CoefXprBinOp(mp_, potD1Fct, densFct, CoefXpr::OP_MULT ) );

     //in case of flow, we need the substantial derivative
      // THIS WOULD WORK IF ONLY CefFunctionBOp WOULD BE DEFINED IN A CONSISTENT MANNER!
      // USAGE OF THIS COEFFUNCTION in MechPDE and FluidPDE PROHIBIT THIS AT THE MOMENT!
     /*if( isComplex_ ) {
       convectiveCoef_.reset(new CoefFunctionBOp<Complex>(feFct, res1));
     } else {
       convectiveCoef_.reset(new CoefFunctionBOp<Double>(feFct, res1));
     }


     PtrCoefFct p2Fct  =
              CoefFunction::Generate( mp_, part,
                                     CoefXprBinOp(mp_, densFct,convectiveCoef_, CoefXpr::OP_MULT ) );

      presFct  =
          CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, p1Fct, p2Fct, CoefXpr::OP_ADD ) );
*/
      DefineFieldResult(presFct, pres);
    } else{
      presFct = this->GetCoefFct(ACOU_PRESSURE);
    }
    
    // === ACOU_SURFACE_PRESSURE ===
    shared_ptr<ResultInfo> surfPres(new ResultInfo);
    surfPres->resultType = ACOU_SURFPRESSURE;
    surfPres->dofNames = "";
    surfPres->unit = MapSolTypeToUnit(ACOU_SURFPRESSURE);
    surfPres->entryType = ResultInfo::SCALAR;
    surfPres->definedOn = ResultInfo::SURF_ELEM;
    DefineFieldResult(presFct, surfPres);
  
    shared_ptr<CoefFunctionFormBased> presGradFct;
    
    // some results are only available in potential formulation (both transient and harmonic) or
    // pressure formulation (harmonic only)
    if( formulation_ == ACOU_POTENTIAL || isComplex_ ){
      // === ACOU_VELOCITY ===
      // for potential formulation
      // Velocity \bm{v}_{\mathrm{a}} = -\nabla \psi_\mathrm{a}
      // or for pressure formulation
      // Velocity \bm{v}_\mathrm{a} = j \frac{1} {\omega* \rho} \nabla p_\mathrm{a}
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = ACOU_VELOCITY;
      vel->dofNames = vecDofNames;
      vel->unit = MapSolTypeToUnit(ACOU_VELOCITY);
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::ELEMENT;
      PtrCoefFct velFct;
      shared_ptr<CoefFunctionFormBased> velFctPot;
      if( formulation_ == ACOU_POTENTIAL) {
        // Velocity v = - grad Psi
        if( isComplex_ ) {
          velFctPot.reset(new CoefFunctionBOp<Complex>(feFct, vel, -1.0));
        } else {
          velFctPot.reset(new CoefFunctionBOp<Double>(feFct, vel, -1.0));
        }
        velFct = velFctPot;
        stiffFormCoefs_.insert(velFctPot);
        DefineFieldResult( velFct, vel );
      }
      else if( formulation_ == ACOU_PRESSURE && isComplex_ ) {
        // Velocity v = j* 1/(omega*rho) * grad(p)
        // a) define gradient of pressure
        presGradFct.reset(new CoefFunctionBOp<Complex>(feFct, vel, 1.0));  
        stiffFormCoefs_.insert(presGradFct);
        
        // b) multiply by factor
        PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
        PtrCoefFct factor = 
            CoefFunction::Generate( mp_, Global::COMPLEX, "0","1/(2*pi*f)");
        PtrCoefFct factor2 = 
            CoefFunction::Generate( mp_, Global::COMPLEX,
                                  CoefXprBinOp(mp_,factor, densFct, CoefXpr::OP_DIV ) );
        velFct = 
            CoefFunction::Generate( mp_,  part,
                                    CoefXprBinOp( mp_, factor2, presGradFct, CoefXpr::OP_MULT ) );
        DefineFieldResult( velFct, vel );
      }

      // === ACOU_NORMAL_VELOCITY ===
      // Normal velocity v_{\mathrm{a,n}} = -\nabla \psi_mathrm{a} \cdot \bm{n}
      shared_ptr<ResultInfo> velNormal(new ResultInfo);
      velNormal->resultType = ACOU_NORMAL_VELOCITY;
      velNormal->dofNames = "";
      velNormal->unit = MapSolTypeToUnit(ACOU_NORMAL_VELOCITY);
      velNormal->entryType = ResultInfo::SCALAR;
      velNormal->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> velFctNormal;
      velFctNormal.reset(new CoefFunctionSurf(true, 1.0, velNormal));
      DefineFieldResult(velFctNormal, velNormal);
      if( formulation_ == ACOU_POTENTIAL ) {
        surfCoefFcts_[velFctNormal] = velFctPot;
      }
      else if ( formulation_ == ACOU_PRESSURE && isComplex_ ) {
        surfCoefFcts_[velFctNormal] = velFct;
      }
      
      // === ACOU_INTENSITY ===
      // Intensity \bm{\mathrm{I}}_{\mathrm{a}} = p_{\mathrm{a}} \overline{\bm{v}_{\mathrm{a}}}
      shared_ptr<ResultInfo> intensity(new ResultInfo);
      intensity->resultType = ACOU_INTENSITY;
      intensity->dofNames = vecDofNames;
      intensity->unit = MapSolTypeToUnit(ACOU_INTENSITY);
      intensity->entryType = ResultInfo::VECTOR;
      intensity->definedOn = ResultInfo::ELEMENT;
      PtrCoefFct intensFct;
      // Intensity I = p * conj(v)
      intensFct = 
          CoefFunction::Generate( mp_, part,
                                CoefXprBinOp(mp_, presFct, velFct, CoefXpr::OP_MULT_CONJ ) );
      DefineFieldResult(intensFct, intensity);
      
      // === ACOU_SURFINTENSITY ===
      // Surface intensity \bm{I}_{\mathrm{a,surf}} =  p_{\mathrm{a}} \bm{v_{\mathrm{a}}} \cdot \bm{n}
      shared_ptr<ResultInfo> surfIntensity(new ResultInfo);
      surfIntensity->resultType = ACOU_SURFINTENSITY;
      surfIntensity->dofNames = vecDofNames;
      surfIntensity->unit = MapSolTypeToUnit(ACOU_SURFINTENSITY);
      surfIntensity->entryType = ResultInfo::VECTOR;
      surfIntensity->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> sIntens;
      sIntens.reset(new CoefFunctionSurf(false, 1.0, surfIntensity));
      DefineFieldResult(sIntens, surfIntensity);
      surfCoefFcts_[sIntens] = intensFct;

      // === ACOU_NORMAL_INTENSITY ===
      // Normal intensity I_{\mathrm{a}} = p_{\mathrm{a}} \bm{v_{\mathrm{a}}} \cdot \bm{n}
      shared_ptr<ResultInfo> intensNormal(new ResultInfo);
      intensNormal->resultType = ACOU_NORMAL_INTENSITY;
      intensNormal->dofNames = "";
      intensNormal->unit = MapSolTypeToUnit(ACOU_NORMAL_INTENSITY);
      intensNormal->entryType = ResultInfo::SCALAR;
      intensNormal->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> sNormIntens;
      sNormIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult( sNormIntens, intensNormal );
      surfCoefFcts_[sNormIntens] = intensFct;

      // === ACOU_POWER ===
      // Power P_{\mathrm{a}} = \int_{\Gamma} I_{\mathrm{a}} \ \mathrm{d} \Gamma
      shared_ptr<ResultInfo> power(new ResultInfo);
      power->resultType = ACOU_POWER;
      power->dofNames = "";
      power->unit = MapSolTypeToUnit(ACOU_POWER);
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
      // Power p = \int_Gamma I *n dGamma
      // Integrate normal intensity
      shared_ptr<ResultFunctor> powerFct;
      if( isComplex_ ) {
        powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens, 
                                                            feFct, power ) );
      } else {
        powerFct.reset(new ResultFunctorIntegrate<Double>(sNormIntens, 
                                                          feFct, power ) );
      }
      resultFunctors_[ACOU_POWER] = powerFct;
      availResults_.insert(power);
    }

    // some results are only available in pressure formulation with harmonic simulation
    if ( formulation_ == ACOU_PRESSURE && isComplex_ ) {
      // --------------------------------
      //  COMPLEX & PRESSURE FORMULATION
      // --------------------------------
      
      // === ACOU_POSITION ===
      // Particle Position u_\mathrm{a} = \frac{1} {\rho* \omega^2} \nabla p_\mathrm{a}
      shared_ptr<ResultInfo> pos(new ResultInfo);
      pos->resultType = ACOU_POSITION;
      pos->dofNames = vecDofNames;
      pos->unit = MapSolTypeToUnit(ACOU_POSITION);
      pos->entryType = ResultInfo::VECTOR;
      pos->definedOn = ResultInfo::ELEMENT;
      // u = 1/(rho*omega^2) * grad(p)
      PtrCoefFct one = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
      PtrCoefFct oneOverOmega2rho = CoefFunction::Generate( mp_, part,
          CoefXprBinOp( mp_, one,
            CoefXprBinOp(mp_,CoefFunction::Generate( mp_, Global::REAL, "4*pi*pi*f*f"), densFct, CoefXpr::OP_MULT ),
          CoefXpr::OP_DIV ));
      PtrCoefFct posFct = CoefFunction::Generate( mp_,  part, CoefXprBinOp( mp_, oneOverOmega2rho, presGradFct, CoefXpr::OP_MULT ) );
      DefineFieldResult( posFct, pos );

      // === ACOU_NORMAL_INTENSITY ===
      // Normal intensity I_{\mathrm{a}} = p_{\mathrm{a}} \bm{v_{\mathrm{a}}} \cdot \bm{n}
      shared_ptr<ResultInfo> intensNormalPW;
      intensNormalPW.reset(new ResultInfo);
      intensNormalPW->resultType = ACOU_NORMAL_INTENSITY_PLANEWAVE;
      intensNormalPW->dofNames = "";
      intensNormalPW->unit = MapSolTypeToUnit(ACOU_NORMAL_INTENSITY_PLANEWAVE);
      intensNormalPW->entryType = ResultInfo::SCALAR;
      intensNormalPW->definedOn = ResultInfo::SURF_ELEM;

      //compute 1/(2*rho*c0)
      //shared_ptr<CoefFunctionMulti> c0Fct  = matCoefs_[ACOU_ELEM_SPEED_OF_SOUND];
      //shared_ptr<CoefFunctionSurf> surfC0(new CoefFunctionSurf(false));
      //surfC0->SetVolumeCoefs( c0Fct->GetRegionCoefs() );

      PtrCoefFct c0Fct = this->GetCoefFct( ACOU_ELEM_SPEED_OF_SOUND);
      PtrCoefFct constVal = CoefFunction::Generate( mp_, Global::REAL, "0.5");
      PtrCoefFct val1 =
                CoefFunction::Generate( mp_, Global::COMPLEX,
                                      CoefXprBinOp(mp_,c0Fct, densFct, CoefXpr::OP_MULT ) );
      PtrCoefFct val2 =
                      CoefFunction::Generate( mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_,constVal, val1, CoefXpr::OP_DIV ) );

      PtrCoefFct velFctPW =
          CoefFunction::Generate( mp_,  part,
                                  CoefXprBinOp( mp_, val2, presFct, CoefXpr::OP_MULT ) );

      PtrCoefFct intensPWfnc =
                      CoefFunction::Generate( mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, velFctPW, presFct, CoefXpr::OP_MULT_CONJ ) );

      // normal acoustic intensity for plane waves
      shared_ptr<CoefFunctionSurf> sNormIntensPW;
      sNormIntensPW.reset(new CoefFunctionSurf(false, 1.0, intensNormalPW));
      DefineFieldResult( sNormIntensPW, intensNormalPW );
      surfCoefFcts_[sNormIntensPW] = intensPWfnc;

      // === ACOU_POWER WITH PLANE WAVE ASSUMPTION: p/vn = \rho c ===
      // Power P = \int_{\Gamma} \frac{p_{\mathrm{a}}^2} {2 \rho c} \ \mathrm{d} \Gamma
      shared_ptr<ResultInfo> powerPW;
      powerPW.reset(new ResultInfo);
      powerPW->resultType = ACOU_POWER_PLANEWAVE;
      powerPW->dofNames = "";
      powerPW->unit = MapSolTypeToUnit(ACOU_POWER_PLANEWAVE);
      powerPW->entryType = ResultInfo::SCALAR;
      powerPW->definedOn = ResultInfo::SURF_REGION;
      shared_ptr<ResultFunctor> powerFctPW;
      powerFctPW.reset(new ResultFunctorIntegrate<Complex>(sNormIntensPW,
                                                        feFct, powerPW ) );
      resultFunctors_[ACOU_POWER_PLANEWAVE] = powerFctPW;
      availResults_.insert(powerPW);
    }
    
    // === ACOUSTIC KINETIC ENERGY ===
    // kinetic Energy w_{\mathrm{a}}^{\mathrm{kin}} = \frac{1} {2}\rho_{0} \bm{v}_{\mathrm{a}} \cdot \bm{v}_{\mathrm{a}}
    shared_ptr<ResultFunctor> keFunc;
    shared_ptr<ResultInfo> kinEnergy(new ResultInfo);
    kinEnergy->resultType = ACOU_KIN_ENERGY;
    kinEnergy->dofNames = "";
    kinEnergy->unit = MapSolTypeToUnit(ACOU_KIN_ENERGY);
    kinEnergy->entryType = ResultInfo::SCALAR;
    kinEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert ( kinEnergy );

    shared_ptr<BaseFeFunction> deriv1vFct;
    if ( formulation_ == ACOU_PRESSURE )
      deriv1vFct = timeDerivFeFunctions_[ACOU_PRESSURE_DERIV_1];
    else
      deriv1vFct = timeDerivFeFunctions_[ACOU_POTENTIAL_DERIV_1];

    if( isComplex_ ) {
      keFunc.reset(new EnergyResultFunctor<Complex>(deriv1vFct, kinEnergy, 0.5));
    } else {
      keFunc.reset(new EnergyResultFunctor<Double>(deriv1vFct, kinEnergy, 0.5));
    }
    resultFunctors_[ACOU_KIN_ENERGY] = keFunc;
    massFormFunctors_.insert(keFunc);

    // === ACOUSTIC POTENTIAL ENERGY ===
    // potential Energy w_{\mathrm{a}}^{\mathrm{pot}} = \frac{p_{\mathrm{a}}^{2}} {2 \rho_{0} c^{2}}
    shared_ptr<ResultFunctor> keFuncPot;
    shared_ptr<ResultInfo> potEnergy(new ResultInfo);
    potEnergy->resultType = ACOU_POT_ENERGY;
    potEnergy->dofNames = "";
    potEnergy->unit = MapSolTypeToUnit(ACOU_POT_ENERGY);
    potEnergy->entryType = ResultInfo::SCALAR;
    potEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert ( potEnergy );
    if( isComplex_ ) {
      keFuncPot.reset(new EnergyResultFunctor<Complex>(feFct, potEnergy, 0.5));
    } else {
      keFuncPot.reset(new EnergyResultFunctor<Double>(feFct, potEnergy, 0.5));
    }
    resultFunctors_[ACOU_POT_ENERGY] = keFuncPot;
    stiffFormFunctors_.insert(keFuncPot);

    //== ACOUSTIC LOAD DDENSITY  ====
    shared_ptr<ResultInfo> loadDensity( new ResultInfo);
    loadDensity->resultType = ACOU_RHS_LOAD_DENSITY;
    loadDensity->dofNames = "";
    loadDensity->unit = "";

    loadDensity->definedOn = ResultInfo::NODE;
    loadDensity->entryType = ResultInfo::SCALAR;

    acousticSourceDensityCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,isComplex_));
    DefineFieldResult( acousticSourceDensityCoef_,loadDensity );

    // === PML Output Parameters ===
    // Vector holding the eigenvalues of the PML matrix
    // this are the diagonal elements for CLASSIC formulation, 
    // or the diagonal elements of the inner (stretching) matrix 
    // for the CURVILINEAR PML formulation 
    shared_ptr<ResultInfo> pmlDampFactor ( new ResultInfo );
    pmlDampFactor->resultType = PML_DAMP_FACTOR;
    pmlDampFactor->dofNames = vecDofNames;
    pmlDampFactor->unit = "";
    pmlDampFactor->definedOn = ResultInfo::ELEMENT;
    pmlDampFactor->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> pmlDampFactorCoefFct(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    matCoefs_[PML_DAMP_FACTOR] = pmlDampFactorCoefFct;
    DefineFieldResult(pmlDampFactorCoefFct, pmlDampFactor);

    // Curvilinear PML: currently only available for harmonic simulation
    if (this->isComplex_) {
      // Tensor holding the complete the PML matrix. 
      // only set for CURVILINEAR PML formulation
      shared_ptr<ResultInfo> pmlTensor ( new ResultInfo );
      pmlTensor->resultType = PML_TENSOR;
      pmlTensor->dofNames = tensorDofNames;
      pmlTensor->unit = "";
      pmlTensor->definedOn = ResultInfo::ELEMENT;
      pmlTensor->entryType = ResultInfo::TENSOR;
      shared_ptr<CoefFunctionMulti> pmlTensorCoefFct(new CoefFunctionMulti(CoefFunction::TENSOR, dim_, dim_, isComplex_));
      matCoefs_[PML_TENSOR] = pmlTensorCoefFct;
      DefineFieldResult(pmlTensorCoefFct, pmlTensor);

      // Scalar holding the determinant of the PML matrix. 
      // only set for CURVILINEAR PML formulation
      shared_ptr<ResultInfo> pmlDeterminant ( new ResultInfo );
      pmlDeterminant->resultType = PML_DETERMINANT;
      pmlDeterminant->dofNames = "";
      pmlDeterminant->unit = "";
      pmlDeterminant->definedOn = ResultInfo::ELEMENT;
      pmlDeterminant->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionMulti> pmlDetCoefFct(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
      matCoefs_[PML_DETERMINANT] = pmlDetCoefFct;
      DefineFieldResult(pmlDetCoefFct, pmlDeterminant);

      // Scalar holding the distance between points and the PML interface. 
      // only set for CURVILINEAR PML formulation
      shared_ptr<ResultInfo> pmlDistance ( new ResultInfo );
      pmlDistance->resultType = PML_DISTANCE;
      pmlDistance->dofNames = "";
      pmlDistance->unit = "";
      pmlDistance->definedOn = ResultInfo::ELEMENT;
      pmlDistance->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionMulti> pmlDistanceCoefFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      matCoefs_[PML_DISTANCE] = pmlDistanceCoefFct;
      DefineFieldResult(pmlDistanceCoefFct, pmlDistance);
    }

    // === AUX Variables for transient PML===
    if(this->isTimeDomPML_){
      if(!this->isAPML_ && dim_ == 3){
        shared_ptr<ResultInfo> pmlScal ( new ResultInfo );
        pmlScal->resultType = ACOU_PMLAUXSCALAR;
        pmlScal->dofNames = "";
        pmlScal->unit = "-";
        pmlScal->definedOn = ResultInfo::NODE;
        pmlScal->entryType = ResultInfo::SCALAR;
        feFunctions_[ACOU_PMLAUXSCALAR]->SetResultInfo(pmlScal);
        results_.Push_back( pmlScal );
        pmlScal->SetFeFunction(feFunctions_[ACOU_PMLAUXSCALAR]);
        DefineFieldResult( feFunctions_[ACOU_PMLAUXSCALAR], pmlScal );
      }

      shared_ptr<ResultInfo> pmlVec ( new ResultInfo );
      pmlVec->resultType = ACOU_PMLAUXVEC;
      pmlVec->dofNames = vecDofNames;
      pmlVec->unit = "-";
      pmlVec->definedOn = ResultInfo::NODE;
      pmlVec->entryType = ResultInfo::VECTOR;
      feFunctions_[ACOU_PMLAUXVEC]->SetResultInfo(pmlVec);
      results_.Push_back( pmlVec );
      pmlVec->SetFeFunction(feFunctions_[ACOU_PMLAUXVEC]);
      DefineFieldResult( feFunctions_[ACOU_PMLAUXVEC], pmlVec );
    }
  }


   void AcousticPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
     //// === MEAN FLUIDMECH VELOCITY ===
     shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
     flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
     flowvelocity->dofNames = dofNames;
     flowvelocity->unit = MapSolTypeToUnit(MEAN_FLUIDMECH_VELOCITY);

     flowvelocity->definedOn = ResultInfo::NODE;
     flowvelocity->entryType = ResultInfo::VECTOR;
     
     meanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
     DefineFieldResult( meanFlowCoef_, flowvelocity );

     results_.Push_back( flowvelocity );
     availResults_.insert( flowvelocity );
     if(myParam_->Get("flowFormulation")->As<std::string>() == "withDivergence"){
       //// === DIVERGENCE OF MEAN FLOW ===
       shared_ptr<ResultInfo> divflowvelocity( new ResultInfo);
       divflowvelocity->resultType = DIV_MEAN_FLUIDMECH_VELOCITY;
       divflowvelocity->dofNames = "";
       divflowvelocity->unit = MapSolTypeToUnit(DIV_MEAN_FLUIDMECH_VELOCITY);

       divflowvelocity->definedOn = ResultInfo::ELEMENT;
       divflowvelocity->entryType = ResultInfo::SCALAR;

       divMeanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1,isComplex_));
       DefineFieldResult( divMeanFlowCoef_, divflowvelocity );

       results_.Push_back( divflowvelocity );
       availResults_.insert( divflowvelocity );
     }
   }

  //! Init the time stepping
  void AcousticPDE::InitTimeStepping(){

    Double alpha = this->myParam_->Get("timeStepAlpha")->As<Double>();
    GLMScheme * scheme1 = new Newmark(0.5,0.25,alpha);
    shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(scheme1,0));
    feFunctions_[formulation_]->SetTimeScheme(acouScheme);

    if(this->isTimeDomPML_){
      //the choice of alpha schemes depends on the damping matrix (e.g. PML)
      GLMScheme * scheme2 = new Newmark(0.5,0.25,alpha);
      shared_ptr<BaseTimeScheme> vecScheme(new TimeSchemeGLM(scheme2,0));
      feFunctions_[ACOU_PMLAUXVEC]->SetTimeScheme(vecScheme);

      if(!this->isAPML_ && dim_ == 3){
        GLMScheme * scheme3 = new Newmark(0.5,0.25,alpha);
        shared_ptr<BaseTimeScheme> scalScheme(new TimeSchemeGLM(scheme3,0));
        feFunctions_[ACOU_PMLAUXSCALAR]->SetTimeScheme(scalScheme);
      }
    }
  }

  void AcousticPDE::ComputeSOS(PtrCoefFct& c0, PtrCoefFct dens, PtrCoefFct blk,
		  PtrCoefFct regionTemp, std::string tempId) {
//TODO 1 erase inconsistent coding style 1/c0^2
	  if ( complexFluidFormulation_ ) {
		  //here c0 is actually 1/c0^2!!
		  c0 = CoefFunction::Generate( mp_,  Global::COMPLEX,
				  CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV ) );
	  }
	  else {
		  if (tempId != "") {
			  if( formulation_ == ACOU_POTENTIAL )
				  EXCEPTION("We need for temperature dependent speed of sound a pressure formulation");

			  // gasR=287.058 J/kg K   ... universal gas constant
			  // kappa=1.402, adabatic exponent for air
        shared_ptr< CoefFunctionConst<Double> > constVal(new CoefFunctionConst<Double>());
        constVal->SetScalar(402.4553160);

			  c0 = CoefFunction::Generate( mp_,  Global::REAL,
					  CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT),
							  CoefXpr::OP_SQRT) );
		  }
		  else {
			  // c0 = sqrt(bulk_modulus / density)
			  c0 = CoefFunction::Generate( mp_,  Global::REAL,
					  CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV),
							  CoefXpr::OP_SQRT) );
		  }
	  }
  }

  void AcousticPDE::ComputeSOS_SQR(PtrCoefFct& cSQR, PtrCoefFct dens, PtrCoefFct blk,
		                           PtrCoefFct regionTemp, std::string tempId) {
    if (tempId != "") {
      if( formulation_ == ACOU_POTENTIAL )
          EXCEPTION("We need for temperature dependent speed of sound a pressure formulation");

      // gasR=287.058 J/kg K   ... universal gas constant
      // kappa=1.402, adabatic exponent for air
      shared_ptr< CoefFunctionConst<Double> > constVal(new CoefFunctionConst<Double>());
      constVal->SetScalar(402.4553160);

      cSQR = CoefFunction::Generate( mp_,  Global::REAL,
                                    CoefXprBinOp(mp_, constVal, regionTemp, CoefXpr::OP_MULT) );
    }
    else {
      // c^2 = bulk_modulus / density
      cSQR = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV) );
    }
  }

  template <UInt DIM, bool IS_COMPLEX>
  void AcousticPDE::DefineConvectiveIntegrators(RegionIdType actRegion, PtrParamNode curRegNode, 
                                                shared_ptr<ElemList> actSDList, PtrCoefFct coeffM)
  {
    std::string flowId = curRegNode->Get("flowId")->As<std::string>();
    if (flowId != "") 
    {
      std::cout << "Assigning convective integrators for AcousticPDE" << std::endl;
      if (complexFluidFormulation_) 
      {
        EXCEPTION("Complex fluid and flow currently not allowed");
      }

      // Get result info object for flow
      shared_ptr<ResultInfo> flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);

      // Add the region information
      PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow", "name", flowId.c_str());

      bool fullForm = false;

      if (myParam_->Get("flowFormulation")->As<std::string>() == "withDivergence") 
      {
        fullForm = true;
      }

      // Read coefficient flow coefficient function for this region and add it to flow functor
      PtrCoefFct regionFlow;
      std::set<UInt> definedDofs;
      bool coefUpdateGeo;
      ReadUserFieldValues(actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType,
                          IS_COMPLEX, regionFlow, definedDofs, coefUpdateGeo);
      meanFlowCoef_->AddRegion(actRegion, regionFlow);

      PtrCoefFct divRegionFlow;
      PtrCoefFct divUFactors;

      if (fullForm) 
      {
        ReadUserFieldValues(actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType,
                            IS_COMPLEX, divRegionFlow, definedDofs, coefUpdateGeo);
        divRegionFlow->SetDerivativeOperation(CoefFunction::VECTOR_DIVERGENCE);
        divMeanFlowCoef_->AddRegion(actRegion, divRegionFlow);
      }

      // Create integrators
      BaseBDBInt* convectiveStiff = nullptr;
      BiLinearForm* convectiveStiffDivU = nullptr;
      BiLinearForm* convectiveDamp = nullptr;
      BiLinearForm* convectiveDampT = nullptr;
      BiLinearForm* convectiveDampDivU = nullptr;

      // Define aliases depending on IS_COMPLEX (template parameter)
      using ScalarType = typename std::conditional<IS_COMPLEX, Complex, Double>::type;
      using ConvectiveOperatorType = typename std::conditional<IS_COMPLEX, ConvectiveOperator<FeH1, DIM, 1, Complex>, ConvectiveOperator<FeH1, DIM, 1>>::type;
      using IdentityOperatorType = typename std::conditional<IS_COMPLEX, IdentityOperator<FeH1, DIM, 1, Complex>, IdentityOperator<FeH1, DIM, 1>>::type;

      convectiveDamp = new ABInt<ScalarType>(new IdentityOperatorType(), new ConvectiveOperatorType(), 
                                              coeffM, 1.0, coefUpdateGeo);
      convectiveDampT = new ABInt<ScalarType>(new ConvectiveOperatorType(), new IdentityOperatorType(), 
                                              coeffM, -1.0, coefUpdateGeo);
      convectiveStiff = new BBInt<ScalarType>(new ConvectiveOperatorType(), coeffM, -1.0, coefUpdateGeo);

      if (fullForm) 
      {
        divUFactors = CoefFunction::Generate(mp_, IS_COMPLEX ? Global::COMPLEX : Global::REAL,
                                              CoefXprBinOp(mp_, coeffM, divRegionFlow, CoefXpr::OP_MULT));
        convectiveDampDivU = new BBInt<ScalarType>(new IdentityOperatorType(), divUFactors, 1.0, coefUpdateGeo);
        convectiveStiffDivU = new ABInt<ScalarType>(new IdentityOperatorType(), new ConvectiveOperatorType(), 
                                                    divUFactors, 1.0, coefUpdateGeo);
      }

      convectiveStiff->SetBCoefFunctionOpB(meanFlowCoef_);
      convectiveStiff->SetName("convectiveStiffPierce");
      convectiveDamp->SetBCoefFunctionOpB(meanFlowCoef_);
      convectiveDamp->SetName("convectiveDampPierce");
      convectiveDampT->SetBCoefFunctionOpA(meanFlowCoef_);
      convectiveDampT->SetName("convectiveDampPierceTransposed");

      convectiveInts_[actRegion] = convectiveStiff;

      BiLinFormContext* convectiveContextStiff = new BiLinFormContext(convectiveStiff, STIFFNESS);
      BiLinFormContext* convectiveContextDamp = new BiLinFormContext(convectiveDamp, DAMPING);
      BiLinFormContext* convectiveContextDampT = new BiLinFormContext(convectiveDampT, DAMPING);

      convectiveContextDamp->SetEntities(actSDList, actSDList);
      convectiveContextDamp->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      convectiveContextDampT->SetEntities(actSDList, actSDList);
      convectiveContextDampT->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      convectiveContextStiff->SetEntities(actSDList, actSDList);
      convectiveContextStiff->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
      assemble_->AddBiLinearForm(convectiveContextDamp);
      assemble_->AddBiLinearForm(convectiveContextDampT);
      assemble_->AddBiLinearForm(convectiveContextStiff);

      if (fullForm) 
      {
        convectiveStiffDivU->SetBCoefFunctionOpB(meanFlowCoef_);
        convectiveStiffDivU->SetName("convectiveStiffPierceDivU");
        convectiveDampDivU->SetName("convectiveDampPierceDivU");
        BiLinFormContext* convectiveContextStiffDivU = new BiLinFormContext(convectiveStiffDivU, STIFFNESS);
        BiLinFormContext* convectiveContextDampDivU = new BiLinFormContext(convectiveDampDivU, DAMPING);

        convectiveContextDampDivU->SetEntities(actSDList, actSDList);
        convectiveContextDampDivU->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);
        convectiveContextStiffDivU->SetEntities(actSDList, actSDList);
        convectiveContextStiffDivU->SetFeFunctions(feFunctions_[formulation_], feFunctions_[formulation_]);

        assemble_->AddBiLinearForm(convectiveContextDampDivU);
        assemble_->AddBiLinearForm(convectiveContextStiffDivU);
      }
    }
  }
}

template void AcousticPDE::DefineTransientPMLInts<2>(shared_ptr<ElemList>, std::string,
    RegionIdType actRegion, std::string tempId);
template void AcousticPDE::DefineTransientPMLInts<3>(shared_ptr<ElemList>, std::string,
    RegionIdType actRegion, std::string tempId);
template void AcousticPDE::DefineConvectiveIntegrators<2, true>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);
template void AcousticPDE::DefineConvectiveIntegrators<2, false>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);
template void AcousticPDE::DefineConvectiveIntegrators<3, true>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);
template void AcousticPDE::DefineConvectiveIntegrators<3, false>(RegionIdType actRegion, PtrParamNode curRegNode, 
    shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);