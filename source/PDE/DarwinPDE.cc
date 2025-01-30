#include <fstream>

#include "DarwinPDE.hh"

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
DEFINE_LOG(darwinPDE, "darwinPDE")


  // **************
  //  Constructor
  // **************
  DarwinPDE::DarwinPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "Darwin";
    pdematerialclass_ = ELECTROMAGNETIC_DARWIN;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true; //true;

    formulationType_ = paramNode->Get("formulation")->As<std::string>();

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if ( !is3d )
      EXCEPTION("DarwinPDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    eps_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
  }


  // *************
  //  Destructor
  // *************
  DarwinPDE::~DarwinPDE() {
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void DarwinPDE::InitNonLin() {

    SinglePDE::InitNonLin();
  }


  // *****************************
  //  Definition of Integrators
  // *****************************
  void DarwinPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<BaseFeFunction> lagrangeMultFeFunc = feFunctions_[LAGRANGE_MULT]; // b Lagrange multiplier
    shared_ptr<FeSpace> magVecPotFeSpace = magVecPotFeFunc->GetFeSpace();
    shared_ptr<FeSpace> elecScalPotFeSpace = elecScalPotFeFunc->GetFeSpace();
    shared_ptr<FeSpace> lagrangeMultFeSpace = lagrangeMultFeFunc->GetFeSpace();

    shared_ptr<BaseFeFunction> lagrangeMult1FeFunc = NULL;
    shared_ptr<FeSpace> lagrangeMult1FeSpace = NULL;
    if(formulationType_ == "Darwin_doubleLagrange"){
        lagrangeMult1FeFunc = feFunctions_[LAGRANGE_MULT_1]; // c Lagrange multiplier
        lagrangeMult1FeSpace = lagrangeMult1FeFunc->GetFeSpace();
    }


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


      // Get polynomial and integration order for electric scalar potential
      std::string elecScalPolyId = curRegNode->Get("elecScalPolyId")->As<std::string>();
      std::string elecScalIntegId = curRegNode->Get("elecScalIntegId")->As<std::string>();
      elecScalPotFeSpace->SetRegionApproximation(actRegion, elecScalPolyId, elecScalIntegId);


      // Get polynomial and integration order for Lagrange multiplier
      std::string lagrangeMultPolyId = curRegNode->Get("lagrangeMultPolyId")->As<std::string>();
      std::string lagrangeMultIntegId = curRegNode->Get("lagrangeMultIntegId")->As<std::string>();
      lagrangeMultFeSpace->SetRegionApproximation(actRegion, lagrangeMultPolyId, lagrangeMultIntegId);


      if(formulationType_ == "Darwin_doubleLagrange"){
		  // Get polynomial and integration order for the second Lagrange multiplier
		  std::string lagrangeMult1PolyId = curRegNode->Get("lagrangeMultPolyId")->As<std::string>();
		  std::string lagrangeMult1IntegId = curRegNode->Get("lagrangeMultIntegId")->As<std::string>();
		  lagrangeMult1FeSpace->SetRegionApproximation(actRegion, lagrangeMultPolyId, lagrangeMultIntegId);
      }


      // Get possible nonlinearities defined in this region
      const StdVector<NonLinType>& matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      const StdVector<NonLinType>& nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // Create new entity list, based on the region
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // Pass entitylist to fespace / fefunction for magnetic vector, scalar potential and Lagrange multiplier(s)
      magVecPotFeFunc->AddEntityList( actSDList );
      elecScalPotFeFunc->AddEntityList( actSDList );
      lagrangeMultFeFunc->AddEntityList( actSDList );
      if(formulationType_ == "Darwin_doubleLagrange"){
    	  lagrangeMult1FeFunc->AddEntityList( actSDList );
      }


      if(matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1){
        EXCEPTION("DarwinPDE does not support nonlinear"
            "(temperatur dependent) electric conductivity yet!\n"
            "But the implementation is not big deal...I promise.\n"
            "Mostly c&p form MagEdgePDE.");
      }



      if ( nonLinTypes.GetSize() > 0 ){
        // =================================================================================
        //  NONLINEAR SECTION
        // =================================================================================
        EXCEPTION("DarwinPDE does not support nonlinear reluctivity yet!";)
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

        // Electric Permittivity
        PtrCoefFct eps = actMat->GetScalCoefFnc( MAG_PERMITTIVITY_SCALAR, Global::REAL);
        // Add material to global, distributed coefficient function
        eps_->AddRegion(actRegion, eps);

        /* ==============================================
         * K_A_A^(nu):
         * nu curl(A) \cdot curl(A)
           ============================================== */
        BaseBDBInt* K_A_A_nu = NULL;
        K_A_A_nu = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), nuNl, 1.0, updatedGeo_) ;
        K_A_A_nu->SetName("K_A_A_nu");

        BiLinFormContext * K_A_A_nuContext = new BiLinFormContext(K_A_A_nu, STIFFNESS );
        K_A_A_nuContext->SetEntities( actSDList, actSDList );
        K_A_A_nuContext->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc );
        assemble_->AddBiLinearForm( K_A_A_nuContext );
        // Add bdb-integrator to global list, needed for flux density evaluation
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,K_A_A_nu) );


        /* ==============================================
         * K_phi_phi^(sigma): JUST AUXILIARY INTEGRATOR FOR CORRECT POSTPROCESSING
         * grad(V) \cdot grad(V)
           ============================================== */
        BaseBDBInt* AUX_phi_phi = NULL;
        AUX_phi_phi = new BBInt<>(new GradientOperator<FeH1,3,1,Double>(), constOne, 1.0, updatedGeo_) ;
        AUX_phi_phi->SetName("AUX_phi_phi");

        BiLinFormContext * AUX_phi_phiContext = new BiLinFormContext(AUX_phi_phi, DAMPING_AUX );
        AUX_phi_phiContext->SetEntities( actSDList, actSDList );
        AUX_phi_phiContext->SetFeFunctions( elecScalPotFeFunc, elecScalPotFeFunc );
        //assemble_->AddBiLinearForm( AUX_phi_phiContext );
        // Add bdb-integrator to global list, needed for gradient evaluation
        bdbIntsAux1_[actRegion] = AUX_phi_phi;



        if( isConducRegion ){

        	/* ==============================================
			 * K_A_phi^(sigma):
			 * \sigma A \cdot grad(V)
			   ============================================== */
        	BiLinearForm* K_A_phi_sigma = NULL;
        	K_A_phi_sigma = new ABInt<>( new IdentityOperator<FeHCurl,3,1,Double>(), new GradientOperator<FeH1,3,1,Double>(),
					conducCoef, 1.0, updatedGeo_);
        	K_A_phi_sigma->SetName("K_A_phi_sigma");

        	BiLinFormContext * K_A_phi_sigmaContext = new BiLinFormContext(K_A_phi_sigma, STIFFNESS );
        	K_A_phi_sigmaContext->SetEntities( actSDList, actSDList );
        	K_A_phi_sigmaContext->SetFeFunctions( magVecPotFeFunc, elecScalPotFeFunc);
        	assemble_->AddBiLinearForm( K_A_phi_sigmaContext );



        	/* ==============================================
			 * K_phi_A^(sigma):
			 * \sigma grad(V) \cdot  A
			   ============================================== */
        	BiLinearForm* K_phi_A_sigma = NULL;
        	K_phi_A_sigma = new ABInt<>( new GradientOperator<FeH1,3,1,Double>(), new IdentityOperator<FeHCurl,3,1,Double>(),
					conducCoef, 1.0, updatedGeo_);
        	K_phi_A_sigma->SetName("K_phi_A_sigma");

        	BiLinFormContext * K_phi_A_sigmaContext = new BiLinFormContext(K_phi_A_sigma, STIFFNESS );
        	K_phi_A_sigmaContext->SetEntities( actSDList, actSDList );
        	K_phi_A_sigmaContext->SetFeFunctions( elecScalPotFeFunc, magVecPotFeFunc);
        	assemble_->AddBiLinearForm( K_phi_A_sigmaContext );





        	/* ==============================================
        	 * K_phi_phi^(sigma):
        	 * -j/omega * \sigma grad(V) \cdot grad(V)
			   ============================================== */
        	BaseBDBInt* K_phi_phi_sigma = NULL;
        	K_phi_phi_sigma = new BBInt<>(new GradientOperator<FeH1,3,1,Double>(), conducCoef, -1.0, updatedGeo_) ;
        	K_phi_phi_sigma->SetName("K_phi_phi_sigma");

        	BiLinFormContext * K_phi_phi_sigmaContext = new BiLinFormContext(K_phi_phi_sigma, DAMPING_AUX );
        	K_phi_phi_sigmaContext->SetEntities( actSDList, actSDList );
        	K_phi_phi_sigmaContext->SetFeFunctions( elecScalPotFeFunc, elecScalPotFeFunc );
        	assemble_->AddBiLinearForm( K_phi_phi_sigmaContext );
        	// Add bdb-integrator to global list, needed for gradient evaluation
        	//bdbIntsAux1_[actRegion] = K_phi_phi_sigma;



        	/* ==============================================
        	 * K_A_A^(sigma):
        	 * \sigma (A) \cdot (A)
           ============================================== */
        	BaseBDBInt *K_A_A_sigma;
//        	K_A_A_sigma = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(), conducCoef,1.0, updatedGeo_);
        	K_A_A_sigma = new BBIntMassEdge<>(new IdentityOperator<FeHCurl,3,1,Double>(), conducCoef,1.0, updatedGeo_);
        	K_A_A_sigma->SetName("K_A_A_sigma");

        	BiLinFormContext * K_A_A_sigmaContext;
        	K_A_A_sigmaContext = new BiLinFormContext(K_A_A_sigma, DAMPING );
        	K_A_A_sigmaContext->SetEntities( actSDList, actSDList );
        	K_A_A_sigmaContext->SetFeFunctions( magVecPotFeFunc, magVecPotFeFunc );
        	assemble_->AddBiLinearForm( K_A_A_sigmaContext );
        }


          /* ==============================================
           * K_phi_phi^(epsilon):
           * \epsilon grad(V) \cdot grad(V)
             ============================================== */
          BaseBDBInt *K_phi_phi_epsilon;

          K_phi_phi_epsilon = new BBInt<>(new GradientOperator<FeH1,3,1,Double>(), eps, 1.0, updatedGeo_) ;
          K_phi_phi_epsilon->SetName("K_phi_phi_epsilon");

          BiLinFormContext * K_phi_phi_epsilonContext = new BiLinFormContext(K_phi_phi_epsilon, STIFFNESS );
          K_phi_phi_epsilonContext->SetEntities( actSDList, actSDList );
          K_phi_phi_epsilonContext->SetFeFunctions( elecScalPotFeFunc, elecScalPotFeFunc );
          assemble_->AddBiLinearForm( K_phi_phi_epsilonContext );



          /* ==============================================
           * K_A_phi^(epsilon):
           * \epsilon (A) \cdot grad(V)
             ============================================== */
          BiLinearForm* K_A_phi_epsilon = NULL;
          K_A_phi_epsilon = new ABInt<>( new IdentityOperator<FeHCurl,3,1,Double>(), new GradientOperator<FeH1,3,1,Double>(),
                                   eps, 1.0, updatedGeo_);
          K_A_phi_epsilon->SetName("K_A_phi_epsilon");

          BiLinFormContext * K_A_phi_epsilonContext = new BiLinFormContext(K_A_phi_epsilon, DAMPING );
          K_A_phi_epsilonContext->SetEntities( actSDList, actSDList );
          K_A_phi_epsilonContext->SetFeFunctions( magVecPotFeFunc, elecScalPotFeFunc );
          assemble_->AddBiLinearForm( K_A_phi_epsilonContext );



          /* ==============================================
           * K_phi_A^(epsilon):
           * \epsilon grad(V) \cdot (A)
             ============================================== */
          BiLinearForm* K_phi_A_epsilon = NULL;
          K_phi_A_epsilon = new ABInt<>( new GradientOperator<FeH1,3,1,Double>(), new IdentityOperator<FeHCurl,3,1,Double>(),
                                   eps, 1.0, updatedGeo_);
          K_phi_A_epsilon->SetName("K_phi_A_epsilon");

          BiLinFormContext * K_phi_A_epsilonContext = new BiLinFormContext(K_phi_A_epsilon, DAMPING );
          K_phi_A_epsilonContext->SetEntities( actSDList, actSDList );
          K_phi_A_epsilonContext->SetFeFunctions( elecScalPotFeFunc, magVecPotFeFunc );
          assemble_->AddBiLinearForm( K_phi_A_epsilonContext );




          /* ==============================================
           * K_A_p^(epsilon):
           * \epsilon (A) \cdot grad(p)
             ============================================== */
          BiLinearForm* K_A_p_epsilon = NULL;
          K_A_p_epsilon = new ABInt<>( new IdentityOperator<FeHCurl,3,1,Double>(), new GradientOperator<FeH1,3,1,Double>(),
                                   eps, -1.0, updatedGeo_);
          K_A_p_epsilon->SetName("K_A_p_epsilon");

          BiLinFormContext * K_A_p_epsilonContext = new BiLinFormContext(K_A_p_epsilon, DAMPING );
          K_A_p_epsilonContext->SetEntities( actSDList, actSDList );
          K_A_p_epsilonContext->SetFeFunctions( magVecPotFeFunc, lagrangeMultFeFunc  );
          assemble_->AddBiLinearForm( K_A_p_epsilonContext );



          /* ==============================================
           * K_p_A^(epsilon):
           * \epsilon  grad(p) \cdot (A)
             ============================================== */
          BiLinearForm* K_p_A_epsilon = NULL;
          K_p_A_epsilon = new ABInt<>( new GradientOperator<FeH1,3,1,Double>(), new IdentityOperator<FeHCurl,3,1,Double>(),
                                   eps, -1.0, updatedGeo_);
          K_p_A_epsilon->SetName("K_p_A_epsilon");

          BiLinFormContext * K_p_A_epsilonContext = new BiLinFormContext(K_p_A_epsilon, DAMPING );
          K_p_A_epsilonContext->SetEntities( actSDList, actSDList );
          K_p_A_epsilonContext->SetFeFunctions( lagrangeMultFeFunc, magVecPotFeFunc  );
          assemble_->AddBiLinearForm( K_p_A_epsilonContext );


          if(formulationType_ == "Darwin_doubleLagrange"){
              /* ==============================================
               * K_A_c^(epsilon):
               * \epsilon (A) \cdot grad(c)
                 ============================================== */
              BiLinearForm* K_A_c_epsilon = NULL;
              K_A_c_epsilon = new ABInt<>( new IdentityOperator<FeHCurl,3,1,Double>(), new GradientOperator<FeH1,3,1,Double>(),
                                       eps, -1.0, updatedGeo_);
              K_A_c_epsilon->SetName("K_A_p_epsilon");

              BiLinFormContext * K_A_c_epsilonContext = new BiLinFormContext(K_A_c_epsilon, DAMPING );
              K_A_c_epsilonContext->SetEntities( actSDList, actSDList );
              K_A_c_epsilonContext->SetFeFunctions( magVecPotFeFunc, lagrangeMult1FeFunc  );
              assemble_->AddBiLinearForm( K_A_c_epsilonContext );


              /* ==============================================
               * K_c_A^(epsilon):
               * \epsilon  grad(c) \cdot (A)
                 ============================================== */
              BiLinearForm* K_c_A_epsilon = NULL;
              K_c_A_epsilon = new ABInt<>( new GradientOperator<FeH1,3,1,Double>(), new IdentityOperator<FeHCurl,3,1,Double>(),
                                       eps, -1.0, updatedGeo_);
              K_c_A_epsilon->SetName("K_c_A_epsilon");

              BiLinFormContext * K_c_A_epsilonContext = new BiLinFormContext(K_c_A_epsilon, DAMPING );
              K_c_A_epsilonContext->SetEntities( actSDList, actSDList );
              K_c_A_epsilonContext->SetFeFunctions( lagrangeMult1FeFunc, magVecPotFeFunc  );
              assemble_->AddBiLinearForm( K_c_A_epsilonContext );



              // Reading the alpha parameter from the xml file
              std::string alpha = myParam_->Get("alpha")->As<std::string>();
              PtrCoefFct alphaCoef = CoefFunction::Generate(mp_, Global::REAL, alpha);

              /* ==============================================
               * K_b_b^(alpha):
               * \alpha  b \cdot b
                 ============================================== */
              BiLinearForm* K_b_b_alpha = NULL;
              K_b_b_alpha = new BBInt<>( new IdentityOperator<FeH1,3,1,Double>(),
                                       alphaCoef, 1.0, updatedGeo_);
              K_b_b_alpha->SetName("K_b_b_alpha");

              BiLinFormContext * K_b_b_alphaContext = new BiLinFormContext(K_b_b_alpha, STIFFNESS );
              K_b_b_alphaContext->SetEntities( actSDList, actSDList );
              K_b_b_alphaContext->SetFeFunctions( lagrangeMultFeFunc, lagrangeMultFeFunc  );
              assemble_->AddBiLinearForm( K_b_b_alphaContext );

              /* ==============================================
               * K_c_c^(alpha):
               * \alpha  c \cdot c
                 ============================================== */
              BiLinearForm* K_c_c_alpha = NULL;
              K_c_c_alpha = new BBInt<>( new IdentityOperator<FeH1,3,1,Double>(),
                                       alphaCoef, 1.0, updatedGeo_);
              K_c_c_alpha->SetName("K_c_c_alpha");

              BiLinFormContext * K_c_c_alphaContext = new BiLinFormContext(K_c_c_alpha, STIFFNESS );
              K_c_c_alphaContext->SetEntities( actSDList, actSDList );
              K_c_c_alphaContext->SetFeFunctions( lagrangeMult1FeFunc, lagrangeMult1FeFunc  );
              assemble_->AddBiLinearForm( K_c_c_alphaContext );


              /* ==============================================
               * K_b_c^(alpha):
               * \alpha  b \cdot c
                 ============================================== */
              BiLinearForm* K_b_c_alpha = NULL;
              K_b_c_alpha = new BBInt<>( new IdentityOperator<FeH1,3,1,Double>(),
                                       alphaCoef, -1.0, updatedGeo_);
              K_b_c_alpha->SetName("K_b_c_alpha");

              BiLinFormContext * K_b_c_alphaContext = new BiLinFormContext(K_b_c_alpha, STIFFNESS );
              K_b_c_alphaContext->SetEntities( actSDList, actSDList );
              K_b_c_alphaContext->SetFeFunctions( lagrangeMultFeFunc, lagrangeMult1FeFunc  );
              K_b_c_alphaContext->SetCounterPart(true);
              assemble_->AddBiLinearForm( K_b_c_alphaContext );

          }



     } // END OF NONLIN/LIN PART
    } // end for regions
  } // end DefineIntegrators



  void DarwinPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(), endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        EXCEPTION("DarwinPDE: Mortar NC interfaces not tested!");
      case NC_NITSCHE:
      {
        /*
         * that's kind of a dirty hack because for Nitsche NC, we need to access the
         * electric conductivity as MAG_CONDUCTIVITY_SCALAR. But this should only be done in
         * the DarwinPDE
         */
        shared_ptr<CoefFunctionMulti> identifier = NULL;
        identifier.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, true));
        if (dim_ == 2)
          EXCEPTION("DarwinPDE possible only in 3D!")
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


  void DarwinPDE::DefineSurfaceIntegrators( ){
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
//    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
//    if( bcNode ) {
//
//      ParamNodeList eNodes = bcNode->GetList( "voltOnEfield" );
//
//      for( UInt i = 0; i < eNodes.GetSize(); i++ ) {
//        std::string regionName = eNodes[i]->Get("name")->As<std::string>();
//        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
//        std::string volRegName = eNodes[i]->Get("volumeRegion")->As<std::string>();
//        std::string voltage = eNodes[i]->Get("value")->As<std::string>();
//
//        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
//        std::set<RegionIdType> volRegions;
//        volRegions.insert(aRegion);
//
//
//        PtrCoefFct voltCoef = CoefFunction::Generate( mp_, Global::REAL, voltage);
//
//
//
//
//
//
//
//        BiLinearForm * voltInt = NULL;
//
//
//
//
//        if( dim_ == 2 ) {
//          voltInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,1,Complex>(),
//                                       new GradientOperator<FeH1,3,1,Complex>(),
//                                       voltCoef, -1.0, volRegions, updatedGeo_);
//        }
//        else {
//          voltInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,1,Double>(),
//                                       new GradientOperator<FeH1,3,1,Double>(),
//                                       voltCoef, -1.0, volRegions, updatedGeo_);
//        }
//
//
//
//        voltInt->SetName("voltIntegrator");
//        BiLinFormContext *voltContext = new BiLinFormContext(voltInt, STIFFNESS );
//
//        voltContext->SetEntities( actSDList, actSDList );
//        voltContext->SetFeFunctions( feFunctions_[ELEC_POTENTIAL] , feFunctions_[ELEC_POTENTIAL]);
//        feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );
//        assemble_->AddBiLinearForm( voltContext );
//
//      }
//    }
  }




  void DarwinPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(darwinPDE) << "Defining rhs load integrators for DarwinPDE";

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;

    StdVector<std::string> vecDofNames = magVecPotFeFunc->GetResultInfo()->dofNames;
    StdVector<std::string> scalDofNames = elecScalPotFeFunc->GetResultInfo()->dofNames;

    bool coefUpdateGeo = true;
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(darwinPDE) << "Reading prescribed flux density";

    ReadRhsExcitation( "fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      EXCEPTION("Currently no rhs for fluxDensity possible in DarwinPDE");
    }

    // ==================
    //  FIELD INTENSITY
    // ==================
    LOG_DBG(darwinPDE) << "Reading prescribed field intensity";

    ReadRhsExcitation( "fieldIntensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      EXCEPTION("Currently no rhs for fieldIntensity possible in DarwinPDE");
    }
  }


  void DarwinPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void DarwinPDE::InitTimeStepping()
  {
    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    // Helper lambda to create the appropriate scheme
    auto makeScheme = [&]() -> GLMScheme* {
      if (integrationScheme)
        return GetXmlDefinedScheme(integrationScheme);
      else
        return new Trapezoidal(1.0);
    };

    // Important: each field gets its own time scheme instance, as otherwise the
    // size of the vectors does not match!
    TimeSchemeGLM::NonLinType nlType = (nonLin_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    feFunctions_[MAG_POTENTIAL] ->SetTimeScheme(shared_ptr<BaseTimeScheme>(new TimeSchemeGLM(makeScheme(), 0, nlType)));
    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(shared_ptr<BaseTimeScheme>(new TimeSchemeGLM(makeScheme(), 0, nlType)));
    feFunctions_[LAGRANGE_MULT] ->SetTimeScheme(shared_ptr<BaseTimeScheme>(new TimeSchemeGLM(makeScheme(), 0, nlType)));
  }

  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void DarwinPDE::ReadCoils() {
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    PtrParamNode coilInfoNode = myInfo_->Get( "coilList", ParamNode::PASS );
    if ( coilNode ){
      EXCEPTION("Currently no coils are supported for DarwinPDE");
    }
  }


  void DarwinPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;
    potInfo->SetFeFunction(feFunctions_[MAG_POTENTIAL]);

    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL], potInfo );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";


    // === ELECTRIC SCALAR POTENTIAL ===
    shared_ptr<BaseFeFunction> scalFct = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<ResultInfo> res2(new ResultInfo);
    res2->resultType = ELEC_POTENTIAL;
    res2->dofNames = "";
    res2->unit = "V";
    res2->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_POTENTIAL);
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


    // === LAGRANGE MULTIPLIER ===
    shared_ptr<BaseFeFunction> lagMplFct = feFunctions_[LAGRANGE_MULT];
    shared_ptr<ResultInfo> res3(new ResultInfo);
    res3->resultType = LAGRANGE_MULT;
    res3->dofNames = "";
    res3->unit = "";
    res3->definedOn = ResultInfo::MapSolTypeToDefinedOn(LAGRANGE_MULT);
    res3->entryType = ResultInfo::SCALAR;
    res3->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
    results_.Push_back( res3);
    availResults_.insert( res3 );
    lagMplFct->SetResultInfo(res3);
    DefineFieldResult( lagMplFct, res3 );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    //hdbcSolNameMap_[ELEC_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[LAGRANGE_MULT] = "lagrangeMultiplier";



    if(formulationType_ == "Darwin_doubleLagrange"){
        // === LAGRANGE MULTIPLIER 1 ===
        shared_ptr<BaseFeFunction> lagMpl1Fct = feFunctions_[LAGRANGE_MULT_1];
        shared_ptr<ResultInfo> res4(new ResultInfo);
        res4->resultType = LAGRANGE_MULT_1;
        res4->dofNames = "";
        res4->unit = "";
        res4->definedOn = ResultInfo::MapSolTypeToDefinedOn(LAGRANGE_MULT_1);
        res4->entryType = ResultInfo::SCALAR;
        res4->SetFeFunction(feFunctions_[LAGRANGE_MULT_1]);
        results_.Push_back( res4);
        availResults_.insert( res4 );
        lagMpl1Fct->SetResultInfo(res4);
        DefineFieldResult( lagMpl1Fct, res4 );

        idbcSolNameMap_[LAGRANGE_MULT_1] = "lagrangeMultiplier1";
    }



    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability ( new ResultInfo );
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_ELEM_PERMEABILITY);
    permeability->entryType = ResultInfo::SCALAR;
    permeability->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult(permFct, permeability);

  }

  void DarwinPDE::DefinePostProcResults() {
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

      // === GRADIENT ELEC SCALAR POTENTIAL ===
      shared_ptr<ResultInfo> gradV(new ResultInfo);
      gradV->resultType = GRAD_ELEC_POTENTIAL;
      gradV->dofNames = vecComponents;
      gradV->unit = "V/m";
      gradV->definedOn = ResultInfo::MapSolTypeToDefinedOn(GRAD_ELEC_POTENTIAL);
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


      // === ELECTRIC FIELD INTENSITY TRANSVERSAL ===
      shared_ptr<ResultInfo> elecIntensT(new ResultInfo);
      elecIntensT->resultType = ELEC_FIELD_INTENSITY_TRANSVERSAL;
      elecIntensT->SetVectorDOFs(dim_, isaxi_);
      elecIntensT->dofNames = vecComponents;
      elecIntensT->unit = "V/m";
      elecIntensT->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_FIELD_INTENSITY_TRANSVERSAL);
      elecIntensT->entryType = ResultInfo::VECTOR;
      elecIntensT->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
      shared_ptr<CoefFunctionMulti> elecIntensTFunc(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( elecIntensTFunc, elecIntensT );


      // === ELECTRIC FIELD INTENSITY LONGITUDINAL ===
      shared_ptr<ResultInfo> elecIntensL(new ResultInfo);
      elecIntensL->resultType = ELEC_FIELD_INTENSITY_LONGITUDINAL;
      elecIntensL->SetVectorDOFs(dim_, isaxi_);
      elecIntensL->dofNames = vecComponents;
      elecIntensL->unit = "V/m";
      elecIntensL->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_FIELD_INTENSITY_LONGITUDINAL);
      elecIntensL->entryType = ResultInfo::VECTOR;
      elecIntensL->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
      shared_ptr<CoefFunctionMulti> elecIntensLFunc(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( elecIntensLFunc, elecIntensL );


      // === ELECTRIC FIELD INTENSITY ===
      shared_ptr<ResultInfo> elecIntens(new ResultInfo);
      elecIntens->resultType = ELEC_FIELD_INTENSITY;
      elecIntens->SetVectorDOFs(dim_, isaxi_);
      elecIntens->dofNames = vecComponents;
      elecIntens->unit = "V/m";
      elecIntens->definedOn = ResultInfo::MapSolTypeToDefinedOn(ELEC_FIELD_INTENSITY);
      elecIntens->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> elecIntensFunc(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( elecIntensFunc, elecIntens );

      // === DISPLACEMENT CURRENT DENSITY ===
      shared_ptr<ResultInfo> displcurrIntens(new ResultInfo);
      displcurrIntens->resultType = DISPLACEMENT_CURRENT_FIELD_INTENSITY;
      displcurrIntens->SetVectorDOFs(dim_, isaxi_);
      displcurrIntens->dofNames = vecComponents;
      displcurrIntens->unit = "";
      displcurrIntens->definedOn = ResultInfo::MapSolTypeToDefinedOn(DISPLACEMENT_CURRENT_FIELD_INTENSITY);
      displcurrIntens->entryType = ResultInfo::VECTOR;
      displcurrIntens->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
      shared_ptr<CoefFunctionMulti> displcurrIntensFunc(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( displcurrIntensFunc, displcurrIntens );

      // === DISPLACEMENT CURRENT ===
      shared_ptr<ResultInfo> displcurr(new ResultInfo);
      displcurr->resultType = DISPLACEMENT_CURRENT_SURF;
      displcurr->SetVectorDOFs(dim_, isaxi_);
      displcurr->dofNames = "";
      displcurr->unit = "";
      displcurr->definedOn = ResultInfo::MapSolTypeToDefinedOn(DISPLACEMENT_CURRENT_SURF);
      displcurr->entryType = ResultInfo::SCALAR;
      displcurr->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
      availResults_.insert( displcurr );
      // first, create normal mapping
      shared_ptr<CoefFunctionSurf> displcurrSurf(new CoefFunctionSurf(true, 1.0, displcurr));
      surfCoefFcts_[displcurrSurf] = displcurrIntensFunc;
      // then, integrate values
      shared_ptr<ResultFunctor> displcurrSurfFunc;
      if( isComplex_ ) {
    	  displcurrSurfFunc.reset(new ResultFunctorIntegrate<Complex>(displcurrSurf, magVecPotFeFct, displcurr ) );
      } else {
    	  displcurrSurfFunc.reset(new ResultFunctorIntegrate<Double>(displcurrSurf, magVecPotFeFct, displcurr ) );
      }
      resultFunctors_[DISPLACEMENT_CURRENT_SURF] = displcurrSurfFunc;

    }

    // === EDDY CURRENT DENSITY ===
    shared_ptr<ResultInfo> eddyJ(new ResultInfo);
    eddyJ->resultType = MAG_EDDY_CURRENT_DENSITY;
    eddyJ->dofNames = vecComponents;
    eddyJ->unit = "A/m^2";
    eddyJ->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_EDDY_CURRENT_DENSITY);
    eddyJ->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> eddyJFunc(
        new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
    DefineFieldResult( eddyJFunc, eddyJ );


    // === EDDY CURRENT (SURFACE RESULT) ===
    shared_ptr<ResultInfo> ec(new ResultInfo());
    ec->resultType = MAG_EDDY_CURRENT;
    ec->dofNames = "";
    ec->unit = "A";
    ec->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_EDDY_CURRENT);
    ec->entryType = ResultInfo::SCALAR;
    availResults_.insert( ec );
    // first, create normal mapping
    shared_ptr<CoefFunctionSurf> ncd(new CoefFunctionSurf(true, 1.0, ec));
    surfCoefFcts_[ncd] = eddyJFunc;
    // then, integrate values
    shared_ptr<ResultFunctor> eddyCurrentFunc;
    if( isComplex_ ) {
      eddyCurrentFunc.reset(new ResultFunctorIntegrate<Complex>(ncd, magVecPotFeFct, ec ) );
    } else {
      eddyCurrentFunc.reset(new ResultFunctorIntegrate<Double>(ncd, magVecPotFeFct, ec ) );
    }
    resultFunctors_[MAG_EDDY_CURRENT] = eddyCurrentFunc;



    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
    fluxDens->resultType = MAG_FLUX_DENSITY;
    fluxDens->dofNames = vecComponents;
    fluxDens->unit = "Vs/m^2";
    fluxDens->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_FLUX_DENSITY);
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
    normFlux->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_NORMAL_FLUX_DENSITY);
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
    flux->definedOn = ResultInfo::MapSolTypeToDefinedOn(MAG_FLUX);
    flux->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    shared_ptr<ResultFunctor> fluxFct;
    if( isComplex_ ) {
      fluxFct.reset(new ResultFunctorIntegrate<Complex>(sNormFDens, magVecPotFeFct, flux ) );
    } else {
      fluxFct.reset(new ResultFunctorIntegrate<Double>(sNormFDens, magVecPotFeFct, flux ) );
    }
    resultFunctors_[MAG_FLUX] = fluxFct;
    availResults_.insert(flux);

  }

  void DarwinPDE::FinalizePostProcResults() {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;


    shared_ptr<BaseFeFunction> magVecPotFeFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFct = feFunctions_[ELEC_POTENTIAL];

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // === ELECTRIC FIELD INTENSITY ===
    // Assemble coefficient function for
    // E_T = - \nabla V
    // E_L = -\frac{\partial A}{\partial t}
    // E = E_T + E_L

    // === EDDY CURRENT DENSITY ===
    // Assemble coefficient function for
    // J = \gamma (E_T + E_L)

    shared_ptr<CoefFunctionMulti> elecIntensTCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY_TRANSVERSAL]);
    shared_ptr<CoefFunctionMulti> elecIntensLCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY_LONGITUDINAL]);
    shared_ptr<CoefFunctionMulti> elecIntensCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);
    shared_ptr<CoefFunctionMulti> jDensLCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_EDDY_CURRENT_DENSITY]);
    shared_ptr<CoefFunctionMulti> dispCurrDensCoef =
            dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[DISPLACEMENT_CURRENT_FIELD_INTENSITY]);



    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "-1.0");
    regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ){
      std::string regionName = ptGrid_->GetRegion().ToString(*regIt);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());

      // Material parameters
      Double conductivity, permittivity;
      materials_[*regIt]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
      PtrCoefFct conducCoef = materials_[*regIt]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR,Global::REAL);
      materials_[*regIt]->GetScalar(permittivity,MAG_PERMITTIVITY_SCALAR,Global::REAL);
      PtrCoefFct permittivityCoef = materials_[*regIt]->GetScalCoefFnc(MAG_PERMITTIVITY_SCALAR,Global::REAL);


      // Electric Field
      PtrCoefFct longitudinal = CoefFunction::Generate( mp_,
          part, CoefXprVecScalOp(mp_, GetCoefFct( GRAD_ELEC_POTENTIAL ), constOne, CoefXpr::OP_MULT));
      PtrCoefFct transversal = CoefFunction::Generate( mp_,
          part, CoefXprVecScalOp(mp_, GetCoefFct( MAG_POTENTIAL_DERIV1 ), constOne, CoefXpr::OP_MULT));
      PtrCoefFct sum = CoefFunction::Generate( mp_,
                part, CoefXprBinOp(mp_, transversal, longitudinal, CoefXpr::OP_ADD));

      elecIntensTCoef->AddRegion(*regIt, transversal);
      elecIntensLCoef->AddRegion(*regIt, longitudinal);
      elecIntensCoef->AddRegion(*regIt, sum);

      // Current Density
      PtrCoefFct curr_tmp = CoefFunction::Generate( mp_,
          part, CoefXprBinOp(mp_, longitudinal, transversal, CoefXpr::OP_ADD));
      PtrCoefFct curr = CoefFunction::Generate( mp_,
          part, CoefXprVecScalOp(mp_, curr_tmp, conducCoef, CoefXpr::OP_MULT));

      jDensLCoef->AddRegion(*regIt, curr);


      // Displacement Current Density dD/dt = epsilon * dE/dt
      // and in frequency-domain hat(D) = j * 2*pi*f * hat(E)
      // but remember: due to the Darwin approximation, dE/dt only consists of the time-derivative of the scalar potential
      PtrCoefFct tmp = CoefFunction::Generate( mp_,
          part, CoefXprVecScalOp(mp_, longitudinal, permittivityCoef, CoefXpr::OP_MULT));
      PtrCoefFct coefFreqFactor = NULL;
      coefFreqFactor = CoefFunction::Generate( mp_, Global::COMPLEX, "-2*pi*f");
      PtrCoefFct dispcurr = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, coefFreqFactor, tmp, CoefXpr::OP_MULT) );

      dispCurrDensCoef->AddRegion(*regIt, dispcurr);


    }

  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  DarwinPDE::CreateFeSpaces(const std::string& formulation,
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

    PtrParamNode lagrangeMultSpaceNode = infoNode->Get("lagrangeMultiplier");
	crSpaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, lagrangeMultSpaceNode, FeSpace::H1, ptGrid_);
	crSpaces[LAGRANGE_MULT]->Init(solStrat_);

	if(formulationType_ == "Darwin_doubleLagrange"){
	    PtrParamNode lagrangeMultSpaceNode = infoNode->Get("lagrangeMultiplier1");
		crSpaces[LAGRANGE_MULT_1] = FeSpace::CreateInstance(myParam_, lagrangeMultSpaceNode, FeSpace::H1, ptGrid_);
		crSpaces[LAGRANGE_MULT_1]->Init(solStrat_);
	}

    return crSpaces;
  }

} // end of namespace

