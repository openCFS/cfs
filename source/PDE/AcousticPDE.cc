#include "AcousticPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
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
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/Mesh/NcInterfaces/BaseNcInterface.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <def_expl_templ_inst.hh>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

  DECLARE_LOG(acousticpde)
   DEFINE_LOG(acousticpde, "pde.acoustic")


  AcousticPDE::AcousticPDE( Grid* aGrid, PtrParamNode paramNode,
                            PtrParamNode infoNode,
                            shared_ptr<SimState> simState, Domain* domain)
              : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

    pdename_           = "acoustic";
    pdematerialclass_  = FLUID;
    nonLin_            = false;
    isMechCoupled_     = false;
    
    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;

    isTimeDomPML_      = false;

    isAPML_ = false;

    //check for pressure or potential formulation
    std::string pdeFormulation = myParam_->Get("formulation")->As<std::string>();
    if(pdeFormulation == "default"){
      formulation_ = ACOU_PRESSURE;
    }else{
      formulation_ = SolutionTypeEnum.Parse(pdeFormulation);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticPDE::CreateFeSpaces( const std::string&  formulation,
                  PtrParamNode infoNode ){

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
      //now define the additional uknowns

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

    //       if( regionNodes.GetSize() > 0 ) {
    //         Info->PrintF( pdename_, "Damping in following region(s)\n" );
    //       }

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
      

      //         // Log to info file
      //         std::string dampString;
      //         Enum2String( dampingList_[actRegionId], dampString );
      //
      //         if( dampingList_[actRegionId] == FRACTIONAL_GL ) {
      //           dampString += "( Gruenwald-Letnikov algorithm )";
      //         }
      //         if( dampingList_[actRegionId] == FRACTIONAL_BLANK ) {
      //           dampString += "( Blanks algorithm )";
      //         }
      //         if( dampingList_[actRegionId] == FRACTIONAL_GL_INT ||
      //             dampingList_[actRegionId] == FRACTIONAL_BLANK_INT ) {
      //           dampString += "(linear interpol. of single past values)";
      //         }
      //         Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(),
      //                       dampString.c_str() );
    }
    //       Info->PrintF( pdename_, "\n" );  
  }

  void AcousticPDE::DefineIntegrators(){

    RegionIdType actRegion;
    // BaseMaterial * actSDMat = NULL;

    //type of geometry
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    //flag indicating frequency PML formulation
    bool harmonicPML = false;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      // actSDMat = it->second;

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
      // Generate coefficient functions
      //=======================================================================

      PtrCoefFct dens = materials_[actRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct blk = materials_[actRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
      
      // c0 = sqrt(bulk_modulus / density)
      PtrCoefFct c0 = 
          CoefFunction::Generate( mp_,  Global::REAL,
                                  CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV),
                                  CoefXpr::OP_SQRT) );
      
      // store coefficient functions
      matCoefs_[ELEM_DENSITY]->AddRegion(actRegion, dens);
      matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->AddRegion( actRegion, c0);

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
          
      // build coefficient for mass matrix as (factor / (c0*c0))
      PtrCoefFct coeffc =
          CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV) );
      PtrCoefFct coeffM =
                CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, factor, coeffc, CoefXpr::OP_MULT) );


      // ====================================================================
      // Take account for pml (frequency domain only)
      // ====================================================================
      shared_ptr<CoefFunction> coeffPMLScal, coeffPMLVec;
      //just coeffunctions for the transformation of jacobians
      shared_ptr<CoefFunction> coeffPMLfactor;
      shared_ptr<CoefFunction> coeffPMLMass;


      if( dampingList_[actRegion] == PML ) {
        std::string dampId;
        curRegNode->GetValue("dampingId",dampId);
        if(analysistype_ == HARMONIC){
          PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",dampId.c_str());
          coeffPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode,c0,actSDList,regions_,true));
          coeffPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode,c0,actSDList,regions_,false));
          // store pml factor
          matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLVec);
          coeffPMLfactor = CoefFunction::Generate( mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_, factor,coeffPMLScal,CoefXpr::OP_MULT));
          coeffPMLMass = CoefFunction::Generate( mp_, Global::COMPLEX,
                                            CoefXprBinOp(mp_,coeffM,coeffPMLScal,CoefXpr::OP_MULT));
          harmonicPML = true;
        }else{
          if(dim_==2)
            DefineTransientPMLInts<2>(actSDList,dampId);
          else
            DefineTransientPMLInts<3>(actSDList,dampId);
        }
      }else{
        harmonicPML = false;
      } // damping == PML


      // ====================================================================
      // standard stiffness integrator
      // ====================================================================
      BaseBDBInt * stiffInt = NULL;
      if( dim_ == 2 ) {
        if(harmonicPML){
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,2,Complex>(),
                                        coeffPMLfactor, 1.0, updatedGeo_ );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), factor, 
                                       1.0, updatedGeo_ );
        }
      }
      else{
        if(harmonicPML){
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,3,Complex>(),
                                        coeffPMLfactor, 1.0, updatedGeo_ );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), factor,
                                       1.0, updatedGeo_ );
        }
      }
      
      stiffInt->SetName("LaplaceIntegrator");

      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );

      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH ) {
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
      bdbInts_[actRegion] = stiffInt;

      // ====================================================================
      // standard mass integrator
      // ====================================================================

      BaseBDBInt *massInt = NULL;
      
      if(dim_==2){
        if(harmonicPML)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1,Complex>(), 
                                       coeffPMLMass, 1.0, updatedGeo_ );
        else
          massInt = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>,coeffM, 
                                      1.0, updatedGeo_ );
      }else{
        if(harmonicPML)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1,Complex>(), coeffPMLMass, 
                                       1.0, updatedGeo_ );
        else
          massInt = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>, coeffM, 
                                      1.0, updatedGeo_ );
      }

      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace() );

      BiLinFormContext *massContext =  new BiLinFormContext(massInt, MASS );

      // Check for damping (mass part)
      if ( dampingList_[actRegion] == RAYLEIGH ) {
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
      // check for flow (Pierce equation)
      // ====================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if(flowId != "") {
        if( dampingList_[actRegion] == PML ) {
          EXCEPTION("PML not available for flow domains!");
        }
        if ( formulation_ != ACOU_POTENTIAL )
          EXCEPTION("Pierce-Equation just possible in velocity potential formulation" );

        // Get result info object for flow
        shared_ptr<ResultInfo> flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY); 
        
        //Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow","name",flowId.c_str());

        bool fullForm = false;
        if(myParam_->Get("flowFormulation")->As<std::string>() == "withDivergence")
          fullForm = true;

        // Read coefficient flow coefficient function for this region and add it to flow functor
        PtrCoefFct regionFlow;
        std::set<UInt> definedDofs;
        bool coefUpdateGeo;
        ReadUserFieldValues( actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType, 
                             isComplex_, regionFlow, definedDofs, coefUpdateGeo );
        meanFlowCoef_->AddRegion( actRegion, regionFlow );


        PtrCoefFct divRegionFlow;
        PtrCoefFct divUFactors;
        if(fullForm){
          ReadUserFieldValues( actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType,
                               isComplex_, divRegionFlow, definedDofs, coefUpdateGeo );
          divRegionFlow->SetDerivativeOperation(CoefFunction::VECTOR_DIVERGENCE);
          divMeanFlowCoef_->AddRegion( actRegion, divRegionFlow );
        }



        //now create the integrators
        BaseBDBInt   *convectiveStiff = NULL;
        BiLinearForm *convectiveStiffDivU = NULL;
        BiLinearForm *convectiveDamp = NULL;
        BiLinearForm *convectiveDampT = NULL;
        BiLinearForm *convectiveDampDivU = NULL;
        if( dim_ == 2 ) {
          if( isComplex_ ) {


            convectiveDamp  = new ABInt<Complex>(new IdentityOperator<FeH1,2,1>(),
                                                 new ConvectiveOperator<FeH1,2,1,Complex>(),
                                                 coeffM, 1.0, coefUpdateGeo);
            convectiveDampT  = new ABInt<Complex>(new ConvectiveOperator<FeH1,2,1,Complex>(),
                                                  new IdentityOperator<FeH1,2,1>(),
                                                  coeffM, -1.0, coefUpdateGeo);

            convectiveStiff = new BBInt<Complex>(new ConvectiveOperator<FeH1,2,1,Complex>(),
                                                 coeffM, -1.0, coefUpdateGeo);
            if(fullForm){
              divUFactors = CoefFunction::Generate( mp_, Global::COMPLEX,
                                                CoefXprBinOp(mp_,coeffM,divRegionFlow,CoefXpr::OP_MULT));
              convectiveDampDivU = new BBInt<Complex>(new IdentityOperator<FeH1,2,1>(),
                                                      divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<Complex>( new IdentityOperator<FeH1,2,1>(),
                                                        new ConvectiveOperator<FeH1,2,1,Complex>(),
                                                        divUFactors, 1.0, coefUpdateGeo);
            }

          } else {

            convectiveDamp  = new ABInt<>(new IdentityOperator<FeH1,2,1>(),
                                          new ConvectiveOperator<FeH1,2,1>(),
                                          coeffM, 1.0, coefUpdateGeo);
            convectiveStiff = new BBInt<>(new ConvectiveOperator<FeH1,2,1>(),
                                          coeffM, -1.0, coefUpdateGeo);

            convectiveDampT  = new ABInt<>(new ConvectiveOperator<FeH1,2,1>(),
                                           new IdentityOperator<FeH1,2,1>(),
                                           coeffM, -1.0, coefUpdateGeo);
            if(fullForm){
              divUFactors = CoefFunction::Generate( mp_, Global::REAL,
                                                CoefXprBinOp(mp_,divRegionFlow,coeffM,CoefXpr::OP_MULT));
              convectiveDampDivU = new BBInt<>(new IdentityOperator<FeH1,2,1>(),
                                                      divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<>( new IdentityOperator<FeH1,2,1>(),
                                                        new ConvectiveOperator<FeH1,2,1>(),
                                                        divUFactors, 1.0, coefUpdateGeo);
            }
          }        
        } else {
          if( isComplex_ ) {

            convectiveDamp  = new ABInt<Complex>(new IdentityOperator<FeH1,3,1>(),
                                                 new ConvectiveOperator<FeH1,3,1,Complex>(),
                                                 coeffM, 1.0, coefUpdateGeo);
            convectiveDampT  = new ABInt<Complex>(new ConvectiveOperator<FeH1,3,1,Complex>(),
                                                  new IdentityOperator<FeH1,3,1>(),
                                                  coeffM, -1.0, coefUpdateGeo);

            convectiveStiff = new BBInt<Complex>(new ConvectiveOperator<FeH1,3,1,Complex>(),
                                                 coeffM, -1.0, coefUpdateGeo);
            if(fullForm){
              divUFactors = CoefFunction::Generate( mp_, Global::COMPLEX,
                                                CoefXprBinOp(mp_,coeffM,divRegionFlow,CoefXpr::OP_MULT));

              convectiveDampDivU = new BBInt<Complex>(new IdentityOperator<FeH1,3,1>(),
                                                      divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<Complex>( new IdentityOperator<FeH1,3,1>(),
                                                        new ConvectiveOperator<FeH1,3,1,Complex>(),
                                                        divUFactors, 1.0, coefUpdateGeo);
            }
          } else {            

            convectiveDamp  = new ABInt<>(new IdentityOperator<FeH1,3,1>(),
                                          new ConvectiveOperator<FeH1,3,1>(),
                                          coeffM, 1.0, coefUpdateGeo);
            convectiveDampT  = new ABInt<>(new ConvectiveOperator<FeH1,3,1>(),
                                           new IdentityOperator<FeH1,3,1>(),
                                           coeffM, -1.0, coefUpdateGeo);

            convectiveStiff = new BBInt<>(new ConvectiveOperator<FeH1,3,1>(),
                                          coeffM, -1.0, coefUpdateGeo);

            if(fullForm){
              divUFactors = CoefFunction::Generate( mp_, Global::REAL,
                                                CoefXprBinOp(mp_,coeffM,divRegionFlow,CoefXpr::OP_MULT));

              convectiveDampDivU = new BBInt<>(new IdentityOperator<FeH1,3,1>(),
                                                      divUFactors, 1.0, coefUpdateGeo);

              convectiveStiffDivU = new ABInt<>( new IdentityOperator<FeH1,3,1>(),
                                                        new ConvectiveOperator<FeH1,3,1>(),
                                                        divUFactors, 1.0, coefUpdateGeo);
            }
          }          
        }
        convectiveStiff->SetBCoefFunctionOpB(meanFlowCoef_);
        convectiveStiff->SetName("convectiveStiffPierce");
        convectiveDamp->SetBCoefFunctionOpB(meanFlowCoef_);
        convectiveDamp->SetName("convectiveDampPierce");
        convectiveDampT->SetBCoefFunctionOpA(meanFlowCoef_);
        convectiveDampT->SetName("convectiveDampPierceTransposed");

        convectiveInts_[actRegion] = convectiveStiff;

        BiLinFormContext *convectiveContextStiff =  new BiLinFormContext(convectiveStiff, STIFFNESS );
        BiLinFormContext *convectiveContextDamp  =  new BiLinFormContext(convectiveDamp, DAMPING );
        BiLinFormContext *convectiveContextDampT  =  new BiLinFormContext(convectiveDampT, DAMPING );

        convectiveContextDamp->SetEntities( actSDList, actSDList );
        convectiveContextDamp->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
        convectiveContextDampT->SetEntities( actSDList, actSDList );
        convectiveContextDampT->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
        convectiveContextStiff->SetEntities( actSDList, actSDList );
        convectiveContextStiff->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
        assemble_->AddBiLinearForm( convectiveContextDamp );
        assemble_->AddBiLinearForm( convectiveContextDampT );
        assemble_->AddBiLinearForm( convectiveContextStiff );

        if(fullForm){
          convectiveStiffDivU->SetBCoefFunctionOpB(meanFlowCoef_);
          convectiveStiffDivU->SetName("convectiveStiffPierceDivU");
          convectiveDampDivU->SetName("convectiveDampPierceDivU");
          BiLinFormContext *convectiveContextStiffDivU =  new BiLinFormContext(convectiveStiffDivU, STIFFNESS );
          BiLinFormContext *convectiveContextDampDivU  =  new BiLinFormContext(convectiveDampDivU, DAMPING );


          convectiveContextDampDivU->SetEntities( actSDList, actSDList );
          convectiveContextDampDivU->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
          convectiveContextStiffDivU->SetEntities( actSDList, actSDList );
          convectiveContextStiffDivU->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);

          assemble_->AddBiLinearForm( convectiveContextDampDivU );
          assemble_->AddBiLinearForm( convectiveContextStiffDivU );
        }

      }
    }
  }

  template<UInt DIM>
  void AcousticPDE::DefineTransientPMLInts(shared_ptr<ElemList> eList, std::string id){

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

    // c0 = sqrt(bulk_modulus / density)
    PtrCoefFct c0 =
        CoefFunction::Generate( mp_,  Global::REAL,
                                CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV),
                                CoefXpr::OP_SQRT) );
    PtrCoefFct coeffc =
        CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, dens, blk, CoefXpr::OP_DIV) );

  // PtrCoefFct coeffc =
  //      CoefFunction::Generate( mp_, Global::REAL,
  //                             CoefXprBinOp( mp_, one,
  //                                           CoefXprBinOp(mp_, c0,c0,CoefXpr::OP_MULT),
  //                                           CoefXpr::OP_DIV ) );



    PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",id.c_str());
    shared_ptr<CoefFunction> coeffPMLVec;
    coeffPMLVec.reset(new CoefFunctionPML<Double>(pmlNode,c0,eList,regions_,true));
    // store pml factor
    matCoefs_[PML_DAMP_FACTOR]->AddRegion(eList->GetRegion(), coeffPMLVec);

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
    vars["a"]  = coeffPMLVec;
    vars["b"] = coeffc;
    var["a"] = coeffPMLVec;

    StdVector<std::string> matAReal;
    StdVector<std::string> matBReal;
    std::string alpha;
    std::string beta;

    if(DIM == 3 ){
      const std::string Amat[] =  { "a_0_R", "0.0" , "0.0" , "0.0" , "a_1_R" , "0.0" , "0.0" , "0.0" , "a_2_R" };
      const std::string Bmat[] = {"( a_0_R - a_1_R - a_2_R )", "0.0","0.0","0.0",  "( a_1_R - a_0_R - a_2_R)","0.0","0.0","0.0","( a_2_R - a_1_R - a_0_R)"};
      alpha = "( (a_0_R + a_1_R + a_2_R ) * b_R )";
      beta = " ( ( (a_0_R * a_1_R) + (a_0_R * a_2_R) + (a_1_R * a_2_R) ) * b_R )";
      matAReal.Import(Amat,9);
      matBReal.Import(Bmat,9);
    }else{
      const std::string Amat[] =  { "a_0_R", "0.0" , "0.0" , "a_1_R" };
      const std::string Bmat[] = {"( a_0_R - a_1_R)", "0.0" , "0.0",  "( a_1_R - a_0_R )"};
      alpha = "( (a_0_R + a_1_R ) * b_R )";
      beta = " ( (a_0_R * a_1_R ) * b_R )";
      matAReal.Import(Amat,4);
      matBReal.Import(Bmat,4);
    }

    coefA->SetTensor(matAReal,DIM,DIM,var);
    coefB->SetTensor(matBReal,DIM,DIM,var);
    coefAlpha->SetScalar(alpha,vars);
    coefBeta->SetScalar(beta,vars);

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
    BaseBDBInt *     divU     = new ABInt<>(new GradientOperator<FeH1,DIM>(), new IdentityOperator<FeH1,DIM,DIM>(),mechAcouFactor,1.0,updatedGeo_);

    BaseBDBInt *     dUdt     = new BBInt<>(new IdentityOperator<FeH1,DIM,DIM>(), one, 1.0, updatedGeo_ );
    BaseBDBInt *     AU       = new BDBInt<>(new IdentityOperator<FeH1,DIM,DIM>(), coefA,1.0, updatedGeo_ );
    BaseBDBInt *     BgradP   = new ADBInt<>(new IdentityOperator<FeH1,DIM,DIM>(),new GradientOperator<FeH1,DIM>(),coefB,1.0,updatedGeo_);

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
      const std::string Cmat[]=  { " (a_1_R * a_2_R) " , "0.0" , "0.0" , "0.0" , " (a_0_R * a_2_R) " , "0.0", "0.0", "0.0", " (a_0_R * a_1_R) "};
      matCReal.Import(Cmat,9);
      std::string gamma = "( (a_0_R * a_1_R * a_2_R) * b_R )";
      coefC->SetTensor(matCReal,3,3,var);
      coefGamma->SetScalar(gamma,vars);

      PtrCoefFct mAcouCorrect_CoefGamma =
          CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp( mp_, mechAcouFactor, coefGamma, CoefXpr::OP_MULT ) );

      //the scalar auxiliary variable is called Nu...
      BaseBDBInt *     dampNu   = new BBInt<>(new IdentityOperator<FeH1,DIM>(), mAcouCorrect_CoefGamma, 1.0, updatedGeo_ );
      BaseBDBInt *     CgradNu  = new ADBInt<>(new IdentityOperator<FeH1,DIM,DIM>(),new GradientOperator<FeH1,DIM>(),coefC,-1.0,updatedGeo_);
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
    if( bcNode ) {
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );

      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        // c0 = sqrt(bulk_modulus / density)
        PtrCoefFct dens = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
        PtrCoefFct blk = materials_[aRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
        PtrCoefFct c0 = 
            CoefFunction::Generate( mp_,  Global::REAL,
                                    CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, blk, dens, 
                                                                 CoefXpr::OP_DIV),
                                                    CoefXpr::OP_SQRT) );
        
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

/*
        // In the case of acou-mech coupling we have to multiply the
        // abc-Integrator matrix with -1
        std::string factor = "1.0";

        if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
          factor = "-1.0";
        }
*/                
        
        // factor for damping matrix: factor / c0
        PtrCoefFct coeffDamp = 
        CoefFunction::Generate( mp_, Global::REAL,
                               CoefXprBinOp(mp_, factor, c0, CoefXpr::OP_DIV ) );
        BiLinearForm * abcInt = NULL;
        if( dim_ == 2 ) {
          abcInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), coeffDamp, 1.0, updatedGeo_ );
        } else {
          abcInt = new BBInt<>(new IdentityOperator<FeH1,3,1>(), coeffDamp, 1.0, updatedGeo_ );
        }

        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, DAMPING );

        abcContext->SetEntities( actSDList, actSDList );
        abcContext->SetFeFunctions( feFunctions_[formulation_] , feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }
    }
  }

  void AcousticPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(acousticpde) << "Defining rhs load integrators for acoustic PDE";
    
    // Get FESpace and FeFunction of mechanical displacement
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
    // integrators by -densiy
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
        PtrCoefFct tmp2 = 
            CoefFunction::Generate( mp_, part, 
                                   CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
        exValue = CoefFunction::Generate( mp_, part, 
                                         CoefXprBinOp(mp_, tmp2,surfDens, 
                                                      CoefXpr::OP_MULT) );
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
        PtrCoefFct tmp2 = 
            CoefFunction::Generate( mp_, part, 
                                   CoefXprBinOp(mp_, tmp, exValue, CoefXpr::OP_MULT) );
        exValue = CoefFunction::Generate( mp_, part, 
                                         CoefXprBinOp(mp_, surfDens, tmp2, 
                                                      CoefXpr::OP_MULT) );
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
    //  PRESSURE (surface) 
    // ===========================
    // Only possible for harmonic potential formulation
    ReadRhsExcitation( "pressure", empty, ResultInfo::SCALAR, isComplex_,
                       ent, coef, coefUpdateGeo );
    if( formulation_ == ACOU_POTENTIAL && isComplex_ ) { 
      for( UInt i = 0; i < ent.GetSize(); ++i ) {

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

    // =====================================
    //  rhsValues for e.g. for aeroacoustics
    // =====================================
    ReadRhsExcitation( "rhsValues", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      coef[i]->SetConservative(true);
      this->rhsFeFunctions_[formulation_]->AddLoadCoefFunction(coef[i], ent[i]);
    }

  }

  void AcousticPDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticPDE::DefinePrimaryResults(){

    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    if ( formulation_ ==  ACOU_PRESSURE) {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "";
      res1->unit = "Pa";
    } else {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "";
      res1->unit = "m^2/s";
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

    CreateMeanFlowFunction(vecDofNames);
    
    // ==============================================
    //  Define CoefFunctions for material parameters
    // ==============================================
    
    // === DENSITY ===
    shared_ptr<ResultInfo> density ( new ResultInfo );
    density->resultType = ELEM_DENSITY;
    density->dofNames = "";
    density->unit = "kg/m^3";
    density->definedOn = ResultInfo::ELEMENT;
    density->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> densFct(new CoefFunctionMulti(CoefFunction::SCALAR,1,1,
                                                                false));
    matCoefs_[ELEM_DENSITY] = densFct;
    DefineFieldResult(densFct, density);
    
    // === SPEED OF SOUND ===
    shared_ptr<ResultInfo> sos ( new ResultInfo );
    sos->resultType = ACOU_ELEM_SPEED_OF_SOUND;
    sos->dofNames = "";
    sos->unit = "m/s";
    sos->definedOn = ResultInfo::ELEMENT;
    sos->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> sosFct(new CoefFunctionMulti(CoefFunction::SCALAR,1,1,
                                                               false));
    matCoefs_[ACOU_ELEM_SPEED_OF_SOUND] = sosFct;
    DefineFieldResult(sosFct, sos);
    
    // === PML DAMPING FACTORS ===
    //if( matCoefs_.find(PML_DAMP_FACTOR) != matCoefs_.end() ) {
    shared_ptr<ResultInfo> pml ( new ResultInfo );
    pml->resultType = PML_DAMP_FACTOR;
    pml->dofNames = vecDofNames;
    //pml->dofNames = "";
    pml->unit = "";
    pml->definedOn = ResultInfo::ELEMENT;
    pml->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> pmlFct(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, 
                                                               isComplex_));
    matCoefs_[PML_DAMP_FACTOR] = pmlFct;
    DefineFieldResult(pmlFct, pml);
    //}

    // === PML AUX Variables ===
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
  
  void AcousticPDE::FinalizePostProcResults(){
    //first call base class method
    SinglePDE::FinalizePostProcResults();

    ////now we additonally add the convective stuff
    //std::map<RegionIdType, BaseBDBInt*>::iterator stiffIt = this->convectiveInts_.begin();
    //for(; stiffIt != convectiveInts_.end(); ++stiffIt ) {
    //  RegionIdType region = stiffIt->first;
    //  BaseBDBInt* bdb = stiffIt->second;
    //  if( !bdb)
    //    continue;
    //
    //  convectiveCoef_->AddIntegrator(bdb, region);
    //}
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
    
    // === PRESSURE / POTENTIAL - 1.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv1(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv1->resultType = ACOU_POTENTIAL_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "m^2/s^2";
    } else {
      deriv1->resultType = ACOU_PRESSURE_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "Pa/s";
    }
    deriv1->entryType = res1->entryType;
    deriv1->definedOn = res1->definedOn;
    availResults_.insert( deriv1 );
    DefineTimeDerivResult( deriv1->resultType, 1, formulation_ );

    // === PRESSURE / POTENTIAL - 2.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv2(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv2->resultType = ACOU_POTENTIAL_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "m^2/s^3";
    } else {
      deriv2->resultType = ACOU_PRESSURE_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "Pa/s^2";
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
      pres->unit = "Pa";
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
    
    
    shared_ptr<ResultInfo> vel, velNormal, intensity, intensNormal, power, pres;
    PtrCoefFct intensFct, velFct;
    shared_ptr<CoefFunctionSurf> sNormIntens, velFctNormal;
    shared_ptr<CoefFunctionFormBased>  presGradFct, velFctPot;
    shared_ptr<ResultFunctor> powerFct;
    
    // some results are only available in potential and / or
    // harmonic simulation
    if( formulation_ == ACOU_POTENTIAL ||
        isComplex_ ){
      // === ACOU_PRESSURE ===
      pres.reset(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = "Pa";
      pres->entryType = ResultInfo::SCALAR;
      pres->definedOn = ResultInfo::ELEMENT;
      
      // === ACOU_VELOCITY ===
      vel.reset(new ResultInfo);
      vel->resultType = ACOU_VELOCITY;
      vel->dofNames = vecDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::ELEMENT;
      
      // === ACOU_NORMAL_VELOCITY ===
      velNormal.reset(new ResultInfo);
      velNormal->resultType = ACOU_NORMAL_VELOCITY;
      velNormal->dofNames = "";
      velNormal->unit = "m/s";
      velNormal->entryType = ResultInfo::SCALAR;
      velNormal->definedOn = ResultInfo::SURF_ELEM;
      
      // === ACOU_INTENSITY ===
      intensity.reset(new ResultInfo);
      intensity->resultType = ACOU_INTENSITY;
      intensity->dofNames = vecDofNames;
      intensity->unit = "W/m^2";
      intensity->entryType = ResultInfo::VECTOR;
      intensity->definedOn = ResultInfo::ELEMENT;
      
      // === ACOU_NORMAL_INTENSITY ===
      intensNormal.reset(new ResultInfo);
      intensNormal->resultType = ACOU_NORMAL_INTENSITY;
      intensNormal->dofNames = "";
      intensNormal->unit = "W/m^2";
      intensNormal->entryType = ResultInfo::SCALAR;
      intensNormal->definedOn = ResultInfo::ELEMENT;
      
      // === ACOU_POWER ===
      power.reset(new ResultInfo);
      power->resultType = ACOU_POWER;
      power->dofNames = "";
      power->unit = "W";
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
    }
    
    if( formulation_ == ACOU_POTENTIAL ) {
      // --------------------------
      //  POTENTENTIAL FORMULATION
      // --------------------------
      // === ACOU_VELOCITY ===
      // Velocity v = - grad Psi
      if( isComplex_ ) {
        velFctPot.reset(new CoefFunctionBOp<Complex>(feFct, vel, -1.0));
      } else {
        velFctPot.reset(new CoefFunctionBOp<Double>(feFct, vel, -1.0));
      }
      velFct = velFctPot;
      stiffFormCoefs_.insert(velFctPot);
      DefineFieldResult( velFct, vel );
      
      // === ACOU_NORMAL_VELOCITY ===
      velFctNormal.reset(new CoefFunctionSurf(true, 1.0, velNormal));
      DefineFieldResult(velFctNormal, velNormal);
      surfCoefFcts_[velFctNormal] = velFctPot;
      
      // === ACOU_INTENSITY ===
      // Intensity I = p * conj(v)
      intensFct = 
          CoefFunction::Generate( mp_, part,
                                 CoefXprBinOp(mp_, presFct, velFct, CoefXpr::OP_MULT_CONJ ) );
      DefineFieldResult(intensFct, intensity);
      
      // === ACOU_NORMAL_INTENSITY ===
      sNormIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult( sNormIntens, intensNormal );
      surfCoefFcts_[sNormIntens] = intensFct;
      
      // === ACOU_POWER ===
      // Power p = \int_Gamma I *n dGamma
      // Integrate normal intensity
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
    
    if( formulation_ == ACOU_PRESSURE && isComplex_ ) {
      // --------------------------------
      //  COMPLEX & PRESSURE FORMULATION
      // --------------------------------
      
      // === ACOU_VELOCITY ===
      
      // Velocity v = j* 1/(omega*rho) * grad(p)
      // a) define gradient of pressure
      if( isComplex_ ) {
        presGradFct.reset(new CoefFunctionBOp<Complex>(feFct, vel, 1.0));  
      } else {
        presGradFct.reset(new CoefFunctionBOp<Double>(feFct, vel, 1.0));
      }
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

      // === ACOU_NORMAL_VELOCITY ===
      velFctNormal.reset(new CoefFunctionSurf(true, 1.0, velNormal));
      DefineFieldResult(velFctNormal, velNormal);
      surfCoefFcts_[velFctNormal] = velFct;
      
      // === ACOU_INTENSITY ===
      // Intensity I = p * v
      intensFct = 
          CoefFunction::Generate( mp_, part,
                                 CoefXprBinOp(mp_, presFct, velFct, CoefXpr::OP_MULT_CONJ ) );
      DefineFieldResult(intensFct, intensity);

      // === ACOU_NORMAL_INTENSITY ===
      sNormIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult( sNormIntens, intensNormal );
      surfCoefFcts_[sNormIntens] = intensFct;
      
      // === ACOU_POWER ===
      // Power p = \int_Gamma I *n dGamma
      //  Integrate normal intensity
      powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens, 
                                                         feFct, power ) );
      resultFunctors_[ACOU_POWER] = powerFct;
      availResults_.insert(power);
    }
    
    // === ACOUSTIC KINETIC ENERGY ===
    shared_ptr<ResultFunctor> keFunc;
    shared_ptr<ResultInfo> kinEnergy(new ResultInfo);
    kinEnergy->resultType = ACOU_KIN_ENERGY;
    kinEnergy->dofNames = "";
    kinEnergy->unit = "Ws";
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
    shared_ptr<ResultFunctor> keFuncPot;
    shared_ptr<ResultInfo> potEnergy(new ResultInfo);
    potEnergy->resultType = ACOU_POT_ENERGY;
    potEnergy->dofNames = "";
    potEnergy->unit = "Ws";
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

  }

   void AcousticPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
     //// === MEAN FLUIDMECH VELOCITY ===
     shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
     flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
     flowvelocity->dofNames = dofNames;
     flowvelocity->unit = "m/s";

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
       divflowvelocity->unit = "1/s";

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

    if(this->isTimeDomPML_){
      //basically the choice for alpha scheme needs to be done everytime we have
      //a damping matrix not just for PML

      //scheme for main unknown
      GLMScheme * scheme1 = new Newmark(0.5,0.25,alpha);
      GLMScheme * scheme2 = new Newmark(0.5,0.25,alpha);
      shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(scheme1,0));
      shared_ptr<BaseTimeScheme> vecScheme(new TimeSchemeGLM(scheme2,0));

      feFunctions_[ACOU_PMLAUXVEC]->SetTimeScheme(vecScheme);
      feFunctions_[formulation_]->SetTimeScheme(acouScheme);

      if(!this->isAPML_ && dim_ == 3){
        GLMScheme * scheme3 = new Newmark(0.5,0.25,alpha);
        shared_ptr<BaseTimeScheme> scalScheme(new TimeSchemeGLM(scheme3,0));
        feFunctions_[ACOU_PMLAUXSCALAR]->SetTimeScheme(scalScheme);
      }
    }else{
      //GLMScheme * scheme1 = new Newmark(0.8,0.4225,-0.3);
      //GLMScheme * scheme1 = new Newmark(0.6,0.3025,alpha);
      GLMScheme * scheme1 = new Newmark(0.5,0.25,alpha);
      shared_ptr<BaseTimeScheme> acouScheme(new TimeSchemeGLM(scheme1,0));
      feFunctions_[formulation_]->SetTimeScheme(acouScheme);
    }
  }

}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template void AcousticPDE::DefineTransientPMLInts<2>(shared_ptr<ElemList>, std::string);
  template void AcousticPDE::DefineTransientPMLInts<3>(shared_ptr<ElemList>, std::string);
#endif
