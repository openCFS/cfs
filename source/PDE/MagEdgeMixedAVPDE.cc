#include <fstream>

#include "MagEdgePDE.hh"
#include "MagEdgeMixedAVPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "FeBasis/H1/FeSpaceH1Hi.hh"
#include "FeBasis/H1/H1Elems.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSUPG.hh"

// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"


// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(magEdgeMixedAVPde, "magEdgeMixedAVPde")


  // **************
  //  Constructor
  // **************
  MagEdgeMixedAVPDE::MagEdgeMixedAVPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagEdgePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdgeMixedAV";
    pdematerialclass_ = ELECTROMAGNETIC;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("MagEdgeMixedAVPDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));

  }


  // *************
  //  Destructor
  // *************
  MagEdgeMixedAVPDE::~MagEdgeMixedAVPDE() {
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagEdgeMixedAVPDE::InitNonLin() {

    SinglePDE::InitNonLin();
  }
  void MagEdgeMixedAVPDE::ReadSpecialBCs() {
    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------

    ReadCoils();

  }

  void MagEdgeMixedAVPDE::DefineIntegrators() {

    this->DefineAVIntegrators();


    DefineCoilIntegrators(1.0);
  }


  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgeMixedAVPDE::DefineAVIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<FeSpace> magVecPotFeSpace = magVecPotFeFunc->GetFeSpace();
    shared_ptr<FeSpace> elecScalPotFeSpace = elecScalPotFeFunc->GetFeSpace();


    PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());

      // Get flag if we need to consider an electric scalar potential in this region
      bool isConducRegion = curRegNode->Get("isConducRegion")->As<bool>();

      // Get polynomial and integration order for magnetic vector potential
      std::string magVecPolyId = curRegNode->Get("magVecPolyId")->As<std::string>();
      std::string magVecIntegId = curRegNode->Get("magVecIntegId")->As<std::string>();
      magVecPotFeSpace->SetRegionApproximation(actRegion, magVecPolyId, magVecIntegId);

      if( isConducRegion ){
        // Get polynomial and integration order for electric scalar potential
        std::string elecScalPolyId = curRegNode->Get("elecScalPolyId")->As<std::string>();
        std::string elecScalIntegId = curRegNode->Get("elecScalIntegId")->As<std::string>();
        elecScalPotFeSpace->SetRegionApproximation(actRegion, elecScalPolyId, elecScalIntegId);
      }

      // Get possible nonlinearities defined in this region
      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // Create new entity list, based on the region
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // Pass entitylist to fespace / fefunction for magnetic vector and electric scalar potential
      magVecPotFeFunc->AddEntityList( actSDList );
      if( isConducRegion ){
        elecScalPotFeFunc->AddEntityList( actSDList );
      }

      if(matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1){
        EXCEPTION("MagEdgeMixedAVPDE does not support nonlinear"
            "(temperatur dependent) electric conductivity yet!\n"
            "But the implementation is not big deal...I promise.\n"
            "Mostly c&p form MagEdgePDE.");
      }

      //check if a veloctyId is assigned 
      std::string velocityId = curRegNode->Get("velocityId")->As<std::string>();

      if ( nonLinTypes.GetSize() > 0 ){
        // =================================================================================
        //  NONLINEAR SECTION
        // =================================================================================
        EXCEPTION("MagEdgeMixedAVPDE does not support nonlinear reluctivity yet!";)
      } else {
        // =================================================================================
        //  LINEAR STIFFNESS SECTION
        // =================================================================================

        /* ==============================================
         * Handling of material parameters
           ============================================== */
        // Magnetic Reluctivity
        PtrCoefFct nuNl = NULL;
        nuNl = actMat->GetScalCoefFnc( MAG_RELUCTIVITY_SCALAR, Global::REAL);
        // Add material to global, distributed reluctivity coefficient function
        reluc_->AddRegion(actRegion, nuNl);

        // Magnetic Permeability
        PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        PtrCoefFct permeability = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, nuNl, CoefXpr::OP_DIV ) );
        matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

        // Electric Conductivity
        Double conductivity;
        materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
        PtrCoefFct conducCoef = materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL);


        /* ==============================================
         * Upper left STIFFNESS part:
         * curl(A) \cdot curl(A’)
           ============================================== */
        BaseBDBInt* stiffUpperLeft = NULL;
        stiffUpperLeft = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), nuNl, 1.0, updatedGeo_) ;
        stiffUpperLeft->SetName("CurlACurlAIntegratorUpperLeft");

        BiLinFormContext * stiffUpperLeftContext = new BiLinFormContext(stiffUpperLeft, STIFFNESS );
        stiffUpperLeftContext->SetEntities( actSDList, actSDList );
        stiffUpperLeftContext->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc );
        assemble_->AddBiLinearForm( stiffUpperLeftContext );
        // Add bdb-integrator to global list, needed for flux density evaluation
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffUpperLeft) );

        /* ==============================================
         * Upper right STIFFNESS part:
         * \sigma grad(V) \cdot A’
           ============================================== */
        if( isConducRegion ){
          BiLinearForm* stiffUpperRight = NULL;
          stiffUpperRight = new ABInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                                        new GradientOperator<FeH1,3,1,Double>(),
                                        conducCoef, 1.0, updatedGeo_);
          stiffUpperRight->SetName("GradVIdentityAIntegratorUpperRight");

          BiLinFormContext * stiffUpperRightContext = new BiLinFormContext(stiffUpperRight, STIFFNESS );
          stiffUpperRightContext->SetEntities( actSDList, actSDList );
          stiffUpperRightContext->SetFeFunctions( magVecPotFeFunc, elecScalPotFeFunc );
          assemble_->AddBiLinearForm( stiffUpperRightContext );
        }


        /* ==============================================
         * Lower right STIFFNESS part:
         * \sigma grad(V’) \cdot grad(V)
           ============================================== */
        if( isConducRegion ){
          BaseBDBInt* stiffLowerRight = NULL;
          stiffLowerRight = new BBInt<>(
              new  GradientOperator<FeH1,3,1,Double>(), conducCoef, 1.0, updatedGeo_) ;
          stiffLowerRight->SetName("GradVGradVIntegratorLowerRight");

          BiLinFormContext * stiffLowerRightContext = new BiLinFormContext(stiffLowerRight, STIFFNESS );
          stiffLowerRightContext->SetEntities( actSDList, actSDList );
          stiffLowerRightContext->SetFeFunctions( elecScalPotFeFunc, elecScalPotFeFunc );
          assemble_->AddBiLinearForm( stiffLowerRightContext );
          // Add bdb-integrator to global list, needed for gradient evaluation
          bdbIntsAux1_[actRegion] = stiffLowerRight;
        }


        /* ==============================================
         * Lower left STIFFNESS part:
         * \sigma grad(V) \cdot A
         * This term is only needed because of symmetry reason!
         * If one uses this formulation for a case with velocity, this term is not allowed to be used!
           ============================================== */ 
        /*   
        if( isConducRegion){
          BiLinearForm* stiffLowerLeft = NULL;
          stiffLowerLeft = new ABInt<>(new GradientOperator<FeH1,3,1,Double>(),
                                        new IdentityOperator<FeHCurl,3,1,Double>(),
                                        conducCoef, 1.0, updatedGeo_);
          stiffLowerLeft->SetName("GradVIdentityAIntegratorLowerLeft");

          BiLinFormContext * stiffLowerLeftContext = new BiLinFormContext(stiffLowerLeft, STIFFNESS );
          stiffLowerLeftContext->SetEntities( actSDList, actSDList );
          stiffLowerLeftContext->SetFeFunctions(elecScalPotFeFunc, magVecPotFeFunc  );
          assemble_->AddBiLinearForm( stiffLowerLeftContext );
        }
      */

      // ====================================================================
      // check for velocity
      // ====================================================================
      if(velocityId != "")
      { 
        bool coefUpdateGeo;
        ReadRegionVelocityField(velocityId,actSDList,actRegion,coefUpdateGeo);
        // Get the stablisation method
        PtrParamNode velNode = myParam_->Get("velocityList")->GetByVal("velocity","name",velocityId.c_str());
        StabilisationType stabilisation = BasePDE::stabilisationType.Parse(velNode->Get("stabilisation")->As<std::string>());


        /* ==============================================
         * Upper VELOCITY part:
         * \gamma \cdot (A’ \cdot (v x curl(A)))
           ============================================== 
        */
        // Create the integrators
        BaseBDBInt   *velocityStiff = NULL;

        // CurlOperatorMag doesn't work with FeH1
        if( isComplex_ )
        {
          EXCEPTION("MagEdgePDE: VelocityStiffInt not implemented for complex case!");
        } else {
          if(dim_ == 2)
            velocityStiff  = new ABInt<>(new IdentityOperator<FeHCurl,2,1>(),new CurlOperatorMag<FeHCurl,2,Double>(),conducCoef, -1.0, coefUpdateGeo);
          else
            velocityStiff  = new ABInt<>(new IdentityOperator<FeHCurl,3,1>(),new CurlOperatorMag<FeHCurl,3,Double>(),conducCoef, -1.0, coefUpdateGeo);
        }
        assert(velocityStiff != NULL);

        //Gets the Velocity Vector, and using this Vector in CurlOperatorMag
        velocityStiff->SetBCoefFunctionOpB(VelocityCoef_);
        velocityStiff->SetName("VelocityStiffInt");
        velocityInts_[actRegion] = velocityStiff;

        BiLinFormContext *VelocityContextStiff =  new BiLinFormContext(velocityStiff, STIFFNESS );
        VelocityContextStiff->SetEntities( actSDList, actSDList );
        VelocityContextStiff->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc);
        assemble_->AddBiLinearForm( VelocityContextStiff );


        /* ==============================================
         * Lower VELOCITY part:
         * \gamma \cdot (grad(V’) \cdot (v x curl(A)))
           ============================================== 
        */
        BaseBDBInt   *velocityStiff2 = NULL;

        // CurlOperatorMag doesn't work with FeH1
        if( isComplex_ )
        {
          EXCEPTION("MagEdgePDE: VelocityStiffInt2 not implemented for complex case!");
        } else {
          if(dim_ == 2)
            velocityStiff2  = new ABInt<>(new GradientOperator<FeH1,2,1>(),new CurlOperatorMag<FeHCurl,2,Double>(),conducCoef, -1.0, coefUpdateGeo);
          else
            velocityStiff2  = new ABInt<>(new GradientOperator<FeH1,3,1>(),new CurlOperatorMag<FeHCurl,3,Double>(),conducCoef, -1.0, coefUpdateGeo);
        }
        assert(velocityStiff2 != NULL);

        // Gets the Velocity Vector, and using this Vector in CurlOperatorMag
        velocityStiff2->SetBCoefFunctionOpB(VelocityCoef_);
        velocityStiff2->SetName("VelocityStiffInt2");
        velocityInts_2[actRegion] = velocityStiff2;

        BiLinFormContext *VelocityContextStiff2 =  new BiLinFormContext(velocityStiff2, STIFFNESS );
        VelocityContextStiff2->SetEntities( actSDList, actSDList );
        VelocityContextStiff2->SetFeFunctions( elecScalPotFeFunc, magVecPotFeFunc);
        assemble_->AddBiLinearForm( VelocityContextStiff2 );

      // Upwinding/Stabilisation Terms
      switch (stabilisation)
      { 
        case StabilisationType::NO_STABILISATION:
        {
          break;
        }
        case StabilisationType::SUPG:
        {
          PtrCoefFct StabilizationFactor;
          // Get the Stabilization factor
          // Magnetic Peclet number is P = conductivity*permeability*abs(velocity)*lElm/2
          // The material coefficient for the peclet number is conductivity*permeability, but in CoefFunctionSUPG
          // the Peclet number is defined as lElm*abs(velocity)/(2*material_coefficient). Therefore we need the 1/(conductivity*permeability)
          PtrCoefFct one  = CoefFunction::Generate(mp_, Global::REAL, "1.0");
          PtrCoefFct materialUpwind = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_,permeability, conducCoef ,CoefXpr::OP_MULT));
          PtrCoefFct peclet_material_coef = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_,one, materialUpwind , CoefXpr::OP_DIV));
          StabilizationFactor.reset(new CoefFunctionSUPG(VelocityCoef_, peclet_material_coef, elecScalPotFeFunc));

          // scalar for Upwind-Terms
          PtrCoefFct coeffUpwindingFactor = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_,StabilizationFactor,conducCoef,CoefXpr::OP_MULT));
          
          // Create the first upwind term
          BaseBDBInt   *UpwindingInt1 = NULL;

          // CurlOperatorMagUpwind doesn't work with FeH1
          if( isComplex_ )
          {
            EXCEPTION("MagEdgePDE: The stablilisation with SUPG is not implemented for the complex case!");
          } else {
            if(dim_ == 2)
              UpwindingInt1  = new ABInt<>(new CurlOperatorMag<FeHCurl,2,Double>(), new CurlOperatorMag<FeHCurl,2,Double>(),coeffUpwindingFactor, 1.0, coefUpdateGeo);
            else
              UpwindingInt1  = new ABInt<>(new CurlOperatorMag<FeHCurl,3,Double>(), new CurlOperatorMag<FeHCurl,3,Double>(),coeffUpwindingFactor, 1.0, coefUpdateGeo);
          }
          assert(UpwindingInt1 != NULL);
          UpwindingInt1->SetBCoefFunctionOpA(VelocityCoef_);
          UpwindingInt1->SetBCoefFunctionOpB(VelocityCoef_);
          UpwindingInt1->SetName("UpwindingInt1");
          UpwindingInt1_1[actRegion] = UpwindingInt1;

          BiLinFormContext *UpwindingContextStiff1 =  new BiLinFormContext(UpwindingInt1, STIFFNESS );
          UpwindingContextStiff1->SetEntities( actSDList, actSDList );
          UpwindingContextStiff1->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc);
          assemble_->AddBiLinearForm( UpwindingContextStiff1 );

          // Create the second upwind term
          BaseBDBInt   *UpwindingInt2 = NULL;
          // CurlOperatorMagUpwind doesn't work with FeH1
          if( isComplex_ )
          {
            EXCEPTION("MagEdgePDE: The stablilisation with SUPG is not implemented for the complex case!");
          } else {
            if(dim_ == 2)
              UpwindingInt2  = new ABInt<>(new CurlOperatorMag<FeHCurl,2,Double>(),new GradientOperator<FeH1,2,1>(),coeffUpwindingFactor, -1.0, coefUpdateGeo);
            else
              UpwindingInt2  = new ABInt<>(new CurlOperatorMag<FeHCurl,3,Double>(),new GradientOperator<FeH1,3,1>(),coeffUpwindingFactor, -1.0, coefUpdateGeo);
          }
          assert(UpwindingInt2 != NULL);

          UpwindingInt2->SetBCoefFunctionOpA(VelocityCoef_);
          UpwindingInt2->SetName("UpwindingInt2");
          UpwindingInt2_1[actRegion] = UpwindingInt2;

          BiLinFormContext *UpwindingContextStiff2 =  new BiLinFormContext(UpwindingInt2, STIFFNESS );
          UpwindingContextStiff2->SetEntities( actSDList, actSDList );
          UpwindingContextStiff2->SetFeFunctions( magVecPotFeFunc, elecScalPotFeFunc);
          assemble_->AddBiLinearForm( UpwindingContextStiff2 );
          break;
        }
        default:
          {
            EXCEPTION("MagEdgeMixedAVPDE: Stabilisation method is not implemented!");
          }
       } // end Stabilisation
      } // end velocityId


        // =================================================================================
        //  LINEAR MASS SECTION
        // =================================================================================
        bool scaleByEdgeSize = false;
        if ( conductivity < 1e-10 || analysistype_ == STATIC ) {
          Matrix<Double> reluc;
          // Get tensor of permeability and determine max. value
          materials_[actRegion]->GetTensor( reluc, MAG_RELUCTIVITY_TENSOR, Global::REAL );
          conductivity =  regularizationFactor * reluc[0][0];
          scaleByEdgeSize = true;
          // Add region to set of "regularized" regions
          regularizedRegions_.insert(actRegion);
        }

        PtrCoefFct conducCoefReg = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
        // Also add material to global, distributed reluctivity coefficient function
        conduc_->AddRegion(actRegion, conducCoefReg);


        /* ==============================================
         * Upper left MASS part:
         * \sigma \frac{dA}{dt} \cdot A’
           ============================================== */
        BaseBDBInt *massUpperLeftInt;
        BiLinFormContext * massUpperLeftContext;

        if ( analysistype_ == STATIC) {
          // we have to guarantee, that we add some mass to curl-curl integrator.
          // Additionally, the integrator gets scaled by the edge size for a uniform
          // conditioning
          massUpperLeftInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),
              conducCoefReg,1.0);
          massUpperLeftInt->SetName("MassIntegratorUpperLeft");
          massUpperLeftContext =  new BiLinFormContext(massUpperLeftInt, STIFFNESS );
        } else {
          // here we add the "normal" mass integrator, which gets not scaled by the
          // edge size
          if( scaleByEdgeSize ) {
            massUpperLeftInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(),
                conducCoefReg,1.0, updatedGeo_);
            massUpperLeftContext = new BiLinFormContext(massUpperLeftInt, STIFFNESS );
          } else {
            massUpperLeftInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(),
                conducCoefReg,1.0, updatedGeo_);
            massUpperLeftContext = new BiLinFormContext(massUpperLeftInt, DAMPING );
          }
          massUpperLeftInt->SetName("MassIntegratorUpperLeft");
        }
        massUpperLeftContext->SetEntities( actSDList, actSDList );
        massUpperLeftContext->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc );
        assemble_->AddBiLinearForm( massUpperLeftContext );
     } // END OF NONLIN/LIN PART
    } // end for regions
  } // end DefineIntegrators



  void MagEdgeMixedAVPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(), endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        EXCEPTION("MagEdgeMixedAVPDE: Mortar NC interfaces not tested!");
      case NC_NITSCHE:
      {
        /*
         * that's kind of a dirty hack because for Nitsche NC, we need to access the
         * electric conductivity as MAG_CONDUCTIVITY_SCALAR. But this should only be done in
         * the MagEdgeMixedAVPDE
         */
        shared_ptr<CoefFunctionMulti> identifier = NULL;
        identifier.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, true));
        if (dim_ == 2)
          EXCEPTION("MagEdgeMixedAVPDE possible only in 3D!")
        else
          //DefineNitscheCoupling<3,1>(ELEC_POTENTIAL, *ncIt, identifier );
          DefineNitscheCoupling<3,1>(MAG_POTENTIAL, *ncIt );
        break;
      }
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }


  void MagEdgeMixedAVPDE::DefineSurfaceIntegrators( ){
    /*
    ========================================================================================
     E x n boundary
     Since we use A-V formulation, E x n = -dA/dt - \nabla V
     usually E x n = 0 and this can be split into two BC's:
      1) A x n = 0 (which is probably already the case, because it's the classical mag BC)
      2) \nabla V x n = 0, can be accomplished by setting V=V0=const. on boundary
     BUT we have to set E x n = -U(t) gradsurf(V), which I also split into two parts:
      1) A x n is already zero, therefore the second term has to do the work
      2) \nabla V x n = -U(t) gradsurf(V)
    ========================================================================================
        */
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {

      ParamNodeList eNodes = bcNode->GetList( "voltOnEfield" );

      for( UInt i = 0; i < eNodes.GetSize(); i++ ) {
        std::string regionName = eNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = eNodes[i]->Get("volumeRegion")->As<std::string>();
        std::string voltage = eNodes[i]->Get("value")->As<std::string>();

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        std::set<RegionIdType> volRegions;
        volRegions.insert(aRegion);

        PtrCoefFct voltCoef = CoefFunction::Generate( mp_, Global::REAL, voltage);
        
        BiLinearForm * voltInt = NULL;
        if( dim_ == 2 ) {
          voltInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,1,Complex>(),
                                       new GradientOperator<FeH1,3,1,Complex>(),
                                       voltCoef, -1.0, volRegions, updatedGeo_);
        }
        else {
          voltInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,1,Double>(),
                                       new GradientOperator<FeH1,3,1,Double>(),
                                       voltCoef, -1.0, volRegions, updatedGeo_);
        }


        voltInt->SetName("voltIntegrator");
        BiLinFormContext *voltContext = new BiLinFormContext(voltInt, STIFFNESS );

        voltContext->SetEntities( actSDList, actSDList );
        voltContext->SetFeFunctions( feFunctions_[ELEC_POTENTIAL] , feFunctions_[ELEC_POTENTIAL]);
        feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( voltContext );

      }
    }
  }

  void MagEdgeMixedAVPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(magEdgeMixedAVPde) << "Defining rhs load integrators for MagEdgeMixedAVPDE";

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    //LinearForm * upperRHSInt = NULL;

    StdVector<std::string> vecDofNames = magVecPotFeFunc->GetResultInfo()->dofNames;
    StdVector<std::string> scalDofNames = elecScalPotFeFunc->GetResultInfo()->dofNames;

    bool coefUpdateGeo = true;
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgeMixedAVPde) << "Reading prescribed flux density";

    ReadRhsExcitation( "fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      EXCEPTION("Currently no rhs for fluxDensity possible in MagEdgeMexedAVPDE");
    }

    // ==================
    //  FIELD INTENSITY
    // ==================
    LOG_DBG(magEdgeMixedAVPde) << "Reading prescribed field intensity";

    ReadRhsExcitation( "fieldIntensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      EXCEPTION("Currently no rhs for fieldIntensity possible in MagEdgeMexedAVPDE");
    }
  }


  void MagEdgeMixedAVPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }

  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgeMixedAVPDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    // Important: Create a new time scheme just for the elec potential unknowns, as otherwise the
    // size of the vectors does not match!
    GLMScheme * scheme2 = new Trapezoidal(gamma);
    shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(scheme2, 0, nlType) );
    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme2);
  }

  /*
  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagEdgeMixedAVPDE::ReadCoils() {
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    PtrParamNode coilInfoNode = myInfo_->Get( "coilList", ParamNode::PASS );
    if ( coilNode ){
      EXCEPTION("Currently no coils are supported for MagEdgeMixedAVPDE");
    }
  }
  */

  void MagEdgeMixedAVPDE::DefinePrimaryResults() {
    MagEdgePDE::DefinePrimaryResults();

    // === ELECTRIC SCALAR POTENTIAL ===
    shared_ptr<BaseFeFunction> scalFct = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<ResultInfo> res2(new ResultInfo);
    res2->resultType = ELEC_POTENTIAL;
    res2->dofNames = "";
    res2->unit = "V";
    res2->definedOn = ResultInfo::NODE;
    res2->entryType = ResultInfo::SCALAR;
    res2->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
    results_.Push_back( res2 );
    availResults_.insert( res2 );
    scalFct->SetResultInfo(res2);
    DefineFieldResult( scalFct, res2 );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    //hdbcSolNameMap_[ELEC_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[ELEC_POTENTIAL] = "elecPotential";

  }

  void MagEdgeMixedAVPDE::DefinePostProcResults() {
    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> magVecPotFeFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFct = feFunctions_[ELEC_POTENTIAL];


    // === TIME DERIVATIVES OF PRIMARY RESULTS ===
    if( analysistype_ != STATIC ) {
      // === MAGNETIC VECTOR POTENTIAL - 1ST DERIVATIVE ===
      shared_ptr<ResultInfo> aDot(new ResultInfo);
      aDot->resultType = MAG_POTENTIAL_DERIV1;
      aDot->dofNames = vecComponents;
      aDot->unit = "V/m";
      aDot->definedOn = ResultInfo::ELEMENT;
      aDot->entryType = ResultInfo::VECTOR;
      aDot->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
      availResults_.insert( aDot );
      DefineTimeDerivResult( MAG_POTENTIAL_DERIV1, 1, MAG_POTENTIAL );

    }

      // === GRADIENT ELEC SCALAR POTENTIAL ===
      shared_ptr<ResultInfo> gradV(new ResultInfo);
      gradV->resultType = GRAD_ELEC_POTENTIAL;
      gradV->dofNames = vecComponents;
      gradV->unit = "V/m";
      gradV->definedOn = ResultInfo::ELEMENT;
      gradV->entryType = ResultInfo::VECTOR;
      gradV->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
      availResults_.insert( gradV );
      shared_ptr<CoefFunctionFormBased> gradVFunc;
      if( isComplex_ ) {
        gradVFunc.reset(new CoefFunctionBOp<Complex>(elecScalPotFeFct, gradV));
      } else {
        gradVFunc.reset(new CoefFunctionBOp<Double>(elecScalPotFeFct, gradV));
      }
      DefineFieldResult( gradVFunc, gradV );
      stiffFormCoefsAux1_.insert(gradVFunc);


      // === ELECTRIC FIELD INTENSITY ===
      shared_ptr<ResultInfo> elecIntens(new ResultInfo);
      elecIntens->resultType = ELEC_FIELD_INTENSITY;
      elecIntens->SetVectorDOFs(dim_, isaxi_);
      elecIntens->dofNames = vecComponents;
      elecIntens->unit = "V/m";
      elecIntens->definedOn = ResultInfo::ELEMENT;
      elecIntens->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> elecIntensFunc(
          new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( elecIntensFunc, elecIntens );


      // === EDDY CURRENT DENSITY ===
      shared_ptr<ResultInfo> eddyJ(new ResultInfo);
      eddyJ->resultType = MAG_EDDY_CURRENT_DENSITY;
      eddyJ->dofNames = vecComponents;
      eddyJ->unit = "A/m^2";
      eddyJ->definedOn = ResultInfo::ELEMENT;
      eddyJ->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> eddyJFunc(
          new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( eddyJFunc, eddyJ );


      /* For integrating eddy current density over a surface
       * we have two possibilities:
       * 1) integration in Hcurl FE space MAG_EDDY_CURRENT1
       * 2) integration in H1 FE space MAG_EDDY_CURRENT2
       * Both integrations give the same result
       * but we'll keep it for testing purposes.
       */
      // === EDDY CURRENT (SURFACE RESULT) ===
      shared_ptr<ResultInfo> ec1(new ResultInfo());
      ec1->resultType = MAG_EDDY_CURRENT1;
      ec1->dofNames = "";
      ec1->unit = "A";
      ec1->definedOn = ResultInfo::SURF_REGION;
      ec1->entryType = ResultInfo::SCALAR;
      availResults_.insert( ec1 );

      // first, create normal mapping, -1.0 because we want the inward pointing normal vector
      shared_ptr<CoefFunctionSurf> ncd(new CoefFunctionSurf(true, -1.0, ec1));
      surfCoefFcts_[ncd] = eddyJFunc;

      // then, integrate values
      shared_ptr<ResultFunctor> eddyCurrentFuncMagVecPot;
      if( isComplex_ ) {
        eddyCurrentFuncMagVecPot.reset(new ResultFunctorIntegrate<Complex>(ncd, magVecPotFeFct, ec1 ) );
      } else {
        eddyCurrentFuncMagVecPot.reset(new ResultFunctorIntegrate<Double>(ncd, magVecPotFeFct, ec1 ) );
      }
      resultFunctors_[MAG_EDDY_CURRENT1] = eddyCurrentFuncMagVecPot;



      // === EDDY CURRENT (SURFACE RESULT) ===
      shared_ptr<ResultInfo> ec2(new ResultInfo());
      ec2->resultType = MAG_EDDY_CURRENT2;
      ec2->dofNames = "";
      ec2->unit = "A";
      ec2->definedOn = ResultInfo::SURF_REGION;
      ec2->entryType = ResultInfo::SCALAR;
      availResults_.insert( ec2 );
      // first, create normal mapping, -1.0 because we want the inward pointing normal vector
      shared_ptr<CoefFunctionSurf> ncd2(new CoefFunctionSurf(true, -1.0, ec2));
      surfCoefFctsAux1_[ncd2] = eddyJFunc;
      // then, integrate values
      shared_ptr<ResultFunctor> eddyCurrentFuncElecScalPot;
      if( isComplex_ ) {
        eddyCurrentFuncElecScalPot.reset(new ResultFunctorIntegrate<Complex>(ncd2, elecScalPotFeFct, ec2 ) );
      } else {
        eddyCurrentFuncElecScalPot.reset(new ResultFunctorIntegrate<Double>(ncd2, elecScalPotFeFct, ec2 ) );
      }
      resultFunctors_[MAG_EDDY_CURRENT2] = eddyCurrentFuncElecScalPot;


    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
    fluxDens->resultType = MAG_FLUX_DENSITY;
    fluxDens->dofNames = vecComponents;
    fluxDens->unit = "Vs/m^2";
    fluxDens->definedOn = ResultInfo::ELEMENT;
    fluxDens->entryType = ResultInfo::VECTOR;
    fluxDens->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    shared_ptr<CoefFunctionFormBased> bFunc;
    if( isComplex_ ) {
      bFunc.reset(new CoefFunctionBOp<Complex>(magVecPotFeFct, fluxDens));
    } else {
      bFunc.reset(new CoefFunctionBOp<Double>(magVecPotFeFct, fluxDens));
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
    normFlux->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
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
      fluxFct.reset(new ResultFunctorIntegrate<Complex>(sNormFDens, magVecPotFeFct, flux ) );
    } else {
      fluxFct.reset(new ResultFunctorIntegrate<Double>(sNormFDens, magVecPotFeFct, flux ) );
    }
    resultFunctors_[MAG_FLUX] = fluxFct;
    availResults_.insert(flux);

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


    /* For integrating total current density over a surface
      * integration in Hcurl FE space MAG_EDDY_CURRENT
    */
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

    // === TOTAL CURRENT (SURFACE RESULT) ===
    shared_ptr<ResultInfo> ect(new ResultInfo());
    ect->resultType = MAG_TOTAL_CURRENT;
    ect->dofNames = "";
    ect->unit = "A";
    ect->definedOn = ResultInfo::SURF_REGION;
    ect->entryType = ResultInfo::SCALAR;
    availResults_.insert( ect );

    // first, create normal mapping, -1.0 because we want the inward pointing normal vector
    shared_ptr<CoefFunctionSurf> ncdt(new CoefFunctionSurf(true, -1.0, ect));
    surfCoefFcts_[ncdt] = tcdCoef;

    // then, integrate values
    shared_ptr<ResultFunctor> totalCurrentFuncMagVecPot;
    if( isComplex_ ) {
      totalCurrentFuncMagVecPot.reset(new ResultFunctorIntegrate<Complex>(ncdt, magVecPotFeFct, ect ) );
    } else {
      totalCurrentFuncMagVecPot.reset(new ResultFunctorIntegrate<Double>(ncdt, magVecPotFeFct, ect ) );
    }
    resultFunctors_[MAG_TOTAL_CURRENT] = totalCurrentFuncMagVecPot;

    // === RELUCTIVITY  ===
    shared_ptr<ResultInfo> reluc(new ResultInfo);
    reluc->resultType = MAG_ELEM_RELUCTIVITY;
    reluc->dofNames = "";
    reluc->unit = "Am/Vs";
    reluc->definedOn = ResultInfo::ELEMENT;
    reluc->entryType = ResultInfo::SCALAR;
    DefineFieldResult( reluc_, reluc );

  }

  void MagEdgeMixedAVPDE::FinalizePostProcResults() {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;


    shared_ptr<BaseFeFunction> magVecPotFeFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFct = feFunctions_[ELEC_POTENTIAL];

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // === ELECTRIC FIELD INTENSITY ===
    // Assemble coefficient function for
    // E = -\frac{\partial A}{\partial t} - \nabla V
    shared_ptr<CoefFunctionMulti> elecIntensCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "-1.0");
    regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ){
      std::string regionName = ptGrid_->GetRegion().ToString(*regIt);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      // Get flag if we need to consider an electric scalar potential in this region
      bool isConducRegion = curRegNode->Get("isConducRegion")->As<bool>();
      if( isConducRegion ){
        if( analysistype_ != STATIC ) {
        PtrCoefFct h = CoefFunction::Generate( mp_, part,
            CoefXprBinOp( mp_, GetCoefFct( GRAD_ELEC_POTENTIAL ),
                GetCoefFct( MAG_POTENTIAL_DERIV1 ), CoefXpr::OP_ADD ) );
        PtrCoefFct h2 = CoefFunction::Generate( mp_, part,
            CoefXprVecScalOp(mp_, h, constOne, CoefXpr::OP_MULT));
            elecIntensCoef->AddRegion(*regIt, h2);
        }
        else{
          PtrCoefFct h2 = 
          CoefFunction::Generate( mp_, part, CoefXprVecScalOp(mp_, GetCoefFct( GRAD_ELEC_POTENTIAL ), constOne, CoefXpr::OP_MULT));
          elecIntensCoef->AddRegion(*regIt, h2);
        }
      }
    }


    // === EDDY CURRENT DENSITY ===
    // J = \sigma * E
    shared_ptr<CoefFunctionMulti> eddyJCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_EDDY_CURRENT_DENSITY]);
    regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ){
      Double conductivity;
      materials_[*regIt]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
      PtrCoefFct conducCoef = materials_[*regIt]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL);
      PtrCoefFct jE = CoefFunction::Generate( mp_, part,
                  CoefXprVecScalOp(mp_, elecIntensCoef, conducCoef, CoefXpr::OP_MULT));
      eddyJCoef->AddRegion(*regIt, jE);
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
    PtrCoefFct velocity = GetCoefFct(MEAN_FLUIDMECH_VELOCITY);
    PtrCoefFct fluxDensity = GetCoefFct(MAG_FLUX_DENSITY);
    shared_ptr<CoefFunctionMulti> tcdCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
    // loop over all regions and assemble total current density:
    //  - if region is coil -> take coil current
    //  - if region is no coil and analyis is transient/harmonic -> eddy

    StdVector<RegionIdType>::iterator regIt2 = regions_.Begin();
    for( ; regIt2 != regions_.End(); ++regIt2 ) {
      RegionIdType actRegion = *regIt2;
      if( coilCurrentDens_.find(actRegion) != coilCurrentDens_.end() ) {
        // region is a coil
        tcdCoef->AddRegion( actRegion, coilCurrentDens_[actRegion] );
      } else {
        // region is no coil
        if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC || analysistype_ == STATIC) {
          Double conductivity;
          materials_[*regIt2]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
          PtrCoefFct conducCoef = materials_[*regIt2]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL);
          
          PtrCoefFct F_L = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, velocity, fluxDensity, CoefXpr::OP_CROSS));
          PtrCoefFct jvel = CoefFunction::Generate( mp_, part,
                  CoefXprVecScalOp(mp_, F_L, conducCoef, CoefXpr::OP_MULT));
          PtrCoefFct jTotal = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, jvel, jEddy, CoefXpr::OP_ADD));
          tcdCoef->AddRegion( actRegion, jTotal );
        }
      }
    }


  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeMixedAVPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    PtrParamNode magVecPpotSpaceNode = infoNode->Get("magPotential");
    crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, magVecPpotSpaceNode, FeSpace::HCURL, ptGrid_ );
    crSpaces[MAG_POTENTIAL]->Init(solStrat_);

    PtrParamNode elecScalPotSpaceNode = infoNode->Get("elecPotential");
    crSpaces[ELEC_POTENTIAL] = FeSpace::CreateInstance(myParam_, elecScalPotSpaceNode, FeSpace::H1, ptGrid_);
    crSpaces[ELEC_POTENTIAL]->Init(solStrat_);


    return crSpaces;
  }

} // end of namespace

