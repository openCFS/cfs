// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>

#include "magneticPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/stdSolveStep.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
#include "Forms/curlfieldop.hh"
#include "Forms/forms_header.hh"
#include "Forms/linElecInt.hh"
#include "trapezoidal.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/assemble.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {

  // **************
  //  Constructor
  // **************
  MagPDE::MagPDE( Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    ENTER_FCN( "MagPDE::MagPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magnetic";
    pdematerialclass_ = ELECTROMAGNETIC;
    maxTimeDerivOrder_ = 1;

    // check if we have a 3d setup
    is3d_ = false;
    is3d_ = param->Get("domain")->Get("geometryType")->AsString() == "3d";
  }


  // *************
  //  Destructor
  // *************
  MagPDE::~MagPDE() {
    ENTER_FCN( "MagPDE::~MagPDE" );
  }


  // *********************************
  //  Read special boundary conitions
  // *********************************

  void MagPDE::ReadSpecialBCs() {

    ENTER_FCN( "MagPDE::ReadSpecialBCs" );

    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();
 
    // -----------------------------
    // Check for permanent magnets
    // -----------------------------
    ReadMagnets();
 
  }
  

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagPDE::InitNonLin() {

    ENTER_FCN( "MagPDE::InitNonLin" );

    // Check, if "nonLinList" is present
    ParamNode * nonLinListNode = myParam_->Get("nonLinList", false );
    if( !nonLinListNode) 
      return;

    // Get nonlinear types
    StdVector<ParamNode*> nonLinNodes = nonLinListNode->GetChildren();
    for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {
      
      std::string actTypeString = nonLinNodes[i]->GetName();
      std::string actId = nonLinNodes[i]->Get("id")->AsString();

      NonLinType actType;
      String2Enum( actTypeString, actType );
      
      // check type
      if( actType == PERMEABILITY ) {
        nonLin_ = TRUE;
      }

      nonLinIdType_[actId] = actType;
    }

    // Run over all region and set entry in "regionNonLinId"
    StdVector<ParamNode*> regionNodes = 
      myParam_->Get("regionList")->GetChildren();
    
    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;
    
    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
      
      // get data
      regionNodes[i]->Get( "name", actRegionName );
      regionNodes[i]->Get( "nonLinId", actNonLinId );
      
      if( actNonLinId == "" )
        continue;
      
      actRegionId = ptgrid_->RegionNameToId( actRegionName );
      
      // Check nonLinId was already registerd
      if( nonLinIdType_.find( actNonLinId) == nonLinIdType_.end() ) {
        EXCEPTION( "NonLinearity with id '" << actNonLinId 
                   << "' was not defined in 'nonLinList'" );
      }
      
      regionNonLinId_[actRegionId] = actNonLinId;
      regionNonLinType_[actRegionId] = nonLinIdType_[actNonLinId];
      
      // Log to info file
      std::string nonLinString;
      Enum2String( nonLinIdType_[actNonLinId], nonLinString );
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(), 
                    nonLinString.c_str() );
      
    }

    // Here we need in addition the nonLinMethod_ for the definition
    // of the integrators
    ParamNode * nonLinNode = myParam_->Get("nonLinear", false );
    nonLinMethod_ = "fixPoint";
    if( nonLinNode ) {
      nonLinNode->Get(  "method", nonLinMethod_, false );
    }
    
  }


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagPDE::DefineIntegrators() {

    ENTER_FCN( "MagPDE::DefineIntegerators" );


    // flag for updatedLagrange formulation
    bool upLagrangeForm = true;

    // Loop over all regions this PDE lives on
    for ( UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {

      // create new entity list
      RegionIdType actRegionId = subdoms_[actSD];
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegionId );

      BaseMaterial* actMat = materials_[subdoms_[actSD]];

      if ( regionNonLinType_[actRegionId] != NO_NONLINEARITY ) {

        //read in the BH-curve data and compute the approximation
//         std::string nlfnc = materials_[actRegionId]->GetNonlinFileName();
//         ApproxData *nlinFnc = new SmoothSpline(nlfnc, BH);
//         //ApproxData *nlinFnc = new LinInterpolate(nlfnc);
//         nlinFnc->SetAccuracy( 0.05 );
//         nlinFnc->SetMaxY( 3.0 );   //maximal value of magnetic induction B
//         nlinFnc->CalcBestParameter();
//         nlinFnc->CalcApproximation();
//         nlinFnc->Print();

        BaseForm *curlcurlNL; 
        if (is3d_ ) {
          curlcurlNL = new nLinCurlCurlNode3DInt( actMat, upLagrangeForm );
        }
        else {
          curlcurlNL = new nLinCurlCurlNode2DInt( actMat, isaxi_, upLagrangeForm );
        }

        curlcurlNL->SetNonLinMethod( nonLinMethod_ );      
        curlcurlNL->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      
        BiLinFormContext * stiffContext = 
          new BiLinFormContext( curlcurlNL, STIFFNESS );
        stiffContext->SetPtPdes(this, this);   
        stiffContext->SetResults( results_[0], results_[0],
                                  actSDList, actSDList );     
        assemble_->AddBiLinearForm( stiffContext);
        

        //save bilinearForm
        pdeBilinearForms_[actRegionId][curlcurlNL->GetName()] = curlcurlNL;

        // nonlinear RHS linearform!!
        LinearForm * rhsSource;
        if ( is3d_ ) {
          rhsSource = new nLinMagNode3D_linFormInt( actMat, upLagrangeForm);
        }
        else {
          rhsSource = new nLinMagNode2D_linFormInt( actMat, isaxi_, 
                                                    upLagrangeForm);
        }

        rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        LinearFormContext * rhsContext = 
          new LinearFormContext( rhsSource );
        rhsContext->SetPtPde( this );
        rhsContext->SetResult( results_[0], actSDList );
        assemble_->AddLinearForm( rhsContext );
      }
      else {
        BaseForm *curlcurl;
        if (is3d_ ) {
          curlcurl = new CurlCurlNode3DInt( actMat, upLagrangeForm);
        }
        else {
          curlcurl = new CurlCurlNode2DInt( actMat, isaxi_, upLagrangeForm);
        }

        BiLinFormContext * stiffContext = 
          new BiLinFormContext(curlcurl, STIFFNESS );
	stiffContext->SetPtPdes(this, this);  
        stiffContext->SetResults( results_[0], results_[0],
                                  actSDList, actSDList );      
	assemble_->AddBiLinearForm( stiffContext );

        //save bilinearForm
        pdeBilinearForms_[actRegionId][curlcurl->GetName()] = curlcurl;


        if ( nonLin_ == true ) {
          // for nonlinear RHS linearform we need linear and nonlinear
          // subdomains
          LinearForm * rhsSource;

          if ( is3d_ ) {
            rhsSource = new nLinMagNode3D_linFormInt( actMat, upLagrangeForm );
          }
          else { 
            rhsSource = new nLinMagNode2D_linFormInt( actMat, isaxi_, upLagrangeForm );
          }

          rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
          LinearFormContext * rhsContext = 
            new LinearFormContext( rhsSource );
          rhsContext->SetPtPde( this );
          rhsContext->SetResult( results_[0], actSDList );
          assemble_->AddLinearForm( rhsContext );
        }
      }

      // Mass matrix
      Double conductivity;
      materials_[actRegionId]->GetScalar(conductivity,MAG_CONDUCTIVITY,REAL);

      BaseForm *bilinear_mass;
      if ( is3d_ ) {
        bilinear_mass = new MassInt(conductivity, dim_, isaxi_, upLagrangeForm );
      }
      else {
        bilinear_mass = new MassInt(conductivity, 1, isaxi_,  upLagrangeForm );
      }

      BiLinFormContext * massContext = 
	new BiLinFormContext(bilinear_mass, MASS );
      massContext->SetPtPdes(this, this);
      massContext->SetResults( results_[0], results_[0],
                               actSDList, actSDList );        
      assemble_->AddBiLinearForm( massContext );

      //add scalar potential bilinear-forms and coupling to vector potential in 3d 
      //case, if region is conductive

      if ( conductivity > EPS && is3d_ &&
           analysistype_ != STATIC ) {
        //define bilinear form for scalar potential
        LaplaceInt* bilinear_Scalar = new LaplaceInt(conductivity, isaxi_, 
                                                     upLagrangeForm ); 
        BiLinFormContext * scalarContext = 
          new BiLinFormContext(bilinear_Scalar, MASS );
        scalarContext->SetPtPdes(this, this);
        scalarContext->SetResults( results_[1], results_[1],
                                   actSDList, actSDList );        
        assemble_->AddBiLinearForm( scalarContext );
        
        //define coupling bilinearform
        MagCoupVectorScalarPotentialInt* coupBilinear = 
          new MagCoupVectorScalarPotentialInt(conductivity, upLagrangeForm );
        BiLinFormContext * coupContext = 
          new BiLinFormContext(coupBilinear, MASS );
        coupContext->SetPtPdes(this, this);
        coupContext->SetResults( results_[0], results_[1],
                                 actSDList, actSDList ); 
        coupContext->SetCounterPart( true );       
        assemble_->AddBiLinearForm( coupContext );
        
        // Give result to equation numbering class
        if ( analysistype_ != STATIC ) 
          eqnMap_->AddResult( *results_[1], actSDList );

        // Define additional bilinearform for calculating the gradient of 
        // the scalar potential
        // Note: this integrator is not passed to the assemble class.
        // It is only needed later for calculating the eddy current density
        linElecInt * elecBiLin = new linElecInt( materials_[actRegionId],
                                                 FULL, true );
        pdeBilinearForms_[actRegionId][elecBiLin->GetName()] = elecBiLin;
      }


      // If this subdomain is a coil we have to do special things
      for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
        if ( actRegionId == coilRegionId_[coil] ) {
          std::string factor = coilDef_[coil]->value_ + "/" +
            GenStr(coilDef_[coil]->windingCrossSection_);

          if ( is3d_ ) {
            MechVolForceInt *coilSource3d = 
              new MechVolForceInt( dim_, coilDef_[coil]->phase_,
                                   isaxi_ );


            StdVector<std::string> currDensity(3);
            currDensity[0] = factor + "*" + GenStr(coilDef_[coil]->locFlowDir_[0]);
            currDensity[1] = factor + "*" + GenStr(coilDef_[coil]->locFlowDir_[1]);
            currDensity[2] = factor + "*" + GenStr(coilDef_[coil]->locFlowDir_[2]);
            coilSource3d->SetVolForceVector( currDensity,
                                             coilDef_[coil]->flowCoordSys_,
                                             true, 1.0 );
            LinearFormContext * coilContext = 
              new LinearFormContext( coilSource3d );
            coilContext->SetPtPde( this );
            coilContext->SetResult( results_[0], actSDList );
            assemble_->AddLinearForm( coilContext );
          }
          else {
            LinearForm *coilSource = 
              new VolumeSrcInt( factor, isaxi_, upLagrangeForm  );  
            LinearFormContext * coilContext = 
              new LinearFormContext( coilSource );
            coilContext->SetPtPde( this );
            coilContext->SetResult( results_[0], actSDList );
            assemble_->AddLinearForm( coilContext );
          }
        }
      }

      // check, if this subdomain is a permanent magnet
      for ( UInt perm = 0; perm < magnetsDomain_.GetSize(); perm++ ) {
        if ( actRegionId == magnetsDomain_[perm] ) {

          //change direction of magnetization for 2D, so that we can use the
          //standard GlobalDerivatives and obtain: (curl w) . M
          Vector<Double> magnetization(dim_);
          if ( is3d_ ) {
            magnetization[0] = magnetsOriX_[perm];
            magnetization[1] = magnetsOriY_[perm];
            magnetization[2] = magnetsOriZ_[perm];
          }
          else if (isaxi_) {
            magnetization[0] = magnetsOriY_[perm];
            magnetization[1] = -magnetsOriX_[perm];
          }
          else {
            magnetization[0] = -magnetsOriY_[perm];
            magnetization[1] = magnetsOriX_[perm];
          }
          
          // Get reluctivity for this domain and perform consistency check
          Double reluctivity;
          materials_[subdoms_[actSD]]->GetScalar(reluctivity,MAG_RELUCTIVITY,REAL);

          std::string fncname = "none";
          LinearForm *permSource;
          if ( is3d_ ) {
            permSource = new MagPerm3DInt(magnetization, reluctivity, 
                                          isaxi_, upLagrangeForm );
          }
          else {
            permSource = new MagPerm2DInt(magnetization, reluctivity, 
                                          isaxi_, upLagrangeForm );
          }

          LinearFormContext * permContext = 
            new LinearFormContext( permSource );
          permContext->SetPtPde( this );
          permContext->SetResult( results_[0], actSDList );
          assemble_->AddLinearForm( permContext );
        }
      }

      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );

    }
  }

  void MagPDE::DefineSolveStep()
  {
    ENTER_FCN( "MagPDE::DefineSolveStep" );

    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagPDE :: InitTimeStepping() {
    ENTER_FCN( "MagPDE::InitTimeStepping" );

    TS_alg_ = new Trapezoidal( algsys_ );
  }




  // ***************
  //   CalcResults
  // ***************
  void MagPDE::CalcResults ( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "MagPDE::CalcResults" );

    switch (res->GetResultInfo()->resultType ) {

    case MAG_POTENTIAL:
      if( isComplex_ ) {
        ExtractResult<Complex>( res, sol_ );
      } else {
        ExtractResult<Double>( res, sol_ );
      }
      break;

    case ELEC_POTENTIAL:
      if( isComplex_ ) {
        ExtractResult<Complex>( res, sol_ );
      } else {
        ExtractResult<Double>( res, sol_ );
      }
      break;

    case MAG_FLUX_DENSITY:
      if( isComplex_ ) {
        CalcFluxDensity<Complex>( res );
      } else {
        CalcFluxDensity<Double>( res );
      }
      break;

    case MAG_EDDY_CURRENT:
      if( isComplex_ ) {
        CalcEddyCurrent<Complex>( res );
      } else {
        CalcEddyCurrent<Double>( res );
      }
      break;
      
    case MAG_ENERGY:
      if( isComplex_ ) {
        CalcEnergy<Complex>( res );
      } else {
        CalcEnergy<Double>( res );
      }
      break;

    case MAG_EDDY_POWER:
      if( isComplex_ ) {
        CalcEddyPower<Complex>( res );
      } else {
        CalcEddyPower<Double>( res );
      }
      break;

    case MAG_FORCE_VWP:
      if( !isComplex_ ) {
        CalcForceVWP( res );
      }
      break;

    default:
      Warning( "Resulttype not computable by magnetic PDE",
               __FILE__, __LINE__ );
    }
  }

  

  //  void MagPDE::PostProcess() {

  //     ENTER_FCN( "MagPDE::PostProcess" );

  //     NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  
  //     if (calcEnergy_.GetSize() !=0 )
  //       CalcEnergy();

  //     if (calcBfield_.GetSize() !=0 ) {






  //     if (UIfilename_.size() > 0 && analysistype_ == TRANSIENT) {
  //       Vector<Double> uiSD;
  //       ComputeUI(uiSD);
  //       WriteUI2File(uiSD);
  //     }

  //     //check for force computation
  //     if ( calcForceVWP_.GetSize() > 0 ) {
  //       Vector<Double> totalForce;
  //       Vector<Double> ForceValues(dim_*ForceNodes_.GetSize());

  //       ForceOpVWP_->CalcNodeForce(ForceValues, totalForce);

  //       // put information into storesol

  //       UInt actPos =  0;
  //       for (UInt iNode = 0; iNode < ForceNodes_.GetSize(); iNode++) {
  //         for (UInt iDof = 0; iDof < dim_; iDof++) {
  //           Force_(ForceNodes_[iNode]-1,iDof) = ForceValues[actPos++];
  //         }
  //       }


  //       // write information in .info-file
  //       Info->PrintF(pdename_, "Sum of magnetic force (VWM):\n");
  //       Info->PrintVec(totalForce);
  //     }

  //     // Last but no least trigger postprocessing fromt within script-file
  //   }

  void MagPDE::CalcForceVWP( shared_ptr<BaseResult> result ) {
    ENTER_FCN( "MagPDE::CalcForceVWP" );

    // fetch result object and convert data
    Result<Double> &  actSol = 
      dynamic_cast<Result<Double>&>(*result);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    
    Vector<Double> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );

    Vector<Double> totalForce;
    ForceOpVWP_->CalcNodeForce(actVal, totalForce);
    std::cerr << "totalForce (VWP) is:" << totalForce.Serialize() << std::endl;

  }

  template<class TYPE>
  void MagPDE::CalcFluxDensity( shared_ptr<BaseResult> result ) {
    ENTER_FCN( "MagPDE::CalcFluxDensity" );

     Vector<TYPE> elemFlux(dim_);
    
    // fetch result object and convert data
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*result);
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );

    // loop over elements
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      CalcFluxDensityAtIP( it, 0, elemFlux );
      for( UInt iDof = 0; iDof < dim_; iDof++ ) {
        actVal[it.GetPos()*dim_ + iDof] = elemFlux[iDof];
      }
    }
  }
  
  template<class TYPE>
  void MagPDE::CalcEddyCurrent( shared_ptr<BaseResult> result ) {
    ENTER_FCN( "MagPDE::CalcEddyCurrent" );
    
    
    // fetch result object and convert data
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*result);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    Vector<TYPE> & actVal = actSol.GetVector();
    UInt jEddyDofs = is3d_ ? 3 : 1;
    actVal.Resize( actSol.GetEntityList()->GetSize() * jEddyDofs );
  
    Vector<TYPE> jEddyElem;
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      CalcEddyCurrentAtIP( it, 0, jEddyElem );
      for( UInt iDof = 0; iDof < jEddyDofs; iDof++ ) {
        actVal[it.GetPos()*jEddyDofs + iDof] = jEddyElem[iDof];
      }
    }
  }

  template<class TYPE>
  void MagPDE::CalcEddyPower( shared_ptr<BaseResult> result ) {
    ENTER_FCN( "MagPDE::CalcEddyPower" );

    // fetch result object and convert data
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*result);
    EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
    Vector<TYPE> & actVal = actSol.GetVector();
    
    UInt jEddyDofs = is3d_ ? 3 : 1;
    actVal.Resize( actSol.GetEntityList()->GetSize() * jEddyDofs );
    
    Vector<TYPE> jEddyElem;
    TYPE regionEddyPower, elemEddyPower;

    // iterate over the regions
    for ( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {


      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();
      regionEddyPower = 0.0;
      
      // get conductivity of region
      Double conductivity;
      materials_[regionIt.GetRegion()]
        ->GetScalar(conductivity,MAG_CONDUCTIVITY,REAL);
      
      // iterate over elements
      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        BaseFE * ptElem = elemIt.GetElem()->ptElem;
        ptElem->SetAnsatzFct( results_[0]->fctType );
        const UInt nrIntPts = ptElem->GetNumIntPoints();
        const Vector<Double> & intWeights = ptElem->GetIntWeights();  
        Vector<Double> shapeFnc;
        
        // iterate over all integration points
        TYPE temp;
        elemEddyPower = 0.0;
        Matrix<Double> ptCoord;
        ptgrid_->GetElemNodesCoord( ptCoord, 
                                    elemIt.GetElem()->connect,
                                    true );
        
        for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {  
          
          CalcEddyCurrentAtIP( elemIt, actIntPt, jEddyElem );
          temp = jEddyElem * jEddyElem;
          temp /= conductivity;
          ptElem->GetShFncAtIp(shapeFnc,actIntPt, elemIt.GetElem() );
          Double jacDet = ptElem->CalcJacobianDetAtIp(actIntPt, ptCoord, 
                                                      elemIt.GetElem() );
          if (isaxi_) {
            Vector<Double> CoordAtIP;
            CoordAtIP = ptCoord * shapeFnc;
            jacDet *=  2 * PI * CoordAtIP[0];
          }
          temp *= jacDet * intWeights[actIntPt-1];
          elemEddyPower += temp;
        }
        regionEddyPower += elemEddyPower;
      }
      // store value in vector
      actVal[regionIt.GetPos()] = regionEddyPower;
    }

  }
  
  

  // **************
  //   CalcEnergy
  // **************
  template<class TYPE>
  void MagPDE::CalcEnergy( shared_ptr<BaseResult> result ) {
    
    ENTER_FCN( "MagPDE::CalcEnergy" );
    
    Matrix<Double> elemmat, ptCoord;
    
    StdVector<UInt> connecth;
    StdVector<Integer> Eqns;  
    Vector<TYPE> help;
    BaseForm *bilinear_stiff = NULL; 
    
    // get result
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*result);      
    EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
    
    // resize vector
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() );
    
    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
   
      
      Double reluctivity;
      materials_[regionIt.GetRegion()]
        ->GetScalar(reluctivity,MAG_RELUCTIVITY,REAL);
      
      // Create stiffness integrator
      if ( regionNonLinType_[regionIt.GetRegion()] != NO_NONLINEARITY ) {
        
        //read in the BH-curve data and compute the approximation
        std::string nlfnc = materials_[regionIt.GetRegion()]->GetNonlinFileName();
        ApproxData *nlinFnc = new SmoothSpline(nlfnc);
        nlinFnc->CalcBestParameter();
        nlinFnc->CalcApproximation();
        bilinear_stiff = new nLinCurlCurlNode2DInt( materials_[regionIt.GetRegion()],
                                                    isaxi_);
        
        bilinear_stiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        
        // VERY IMPORTANT: Set nonlinear-method "fixpoint", as otherwise also
        // the Frechet part of the stiffness is calculated!
        bilinear_stiff->SetNonLinMethod( "fixPoint" );
      } else {
        bilinear_stiff = new CurlCurlNode2DInt( materials_[regionIt.GetRegion()], isaxi_);
      }
      
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();
      
      TYPE energy = 0.0;
      Vector<TYPE> magvecpot;         
      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        sol_->GetElemSolution(magvecpot, elemIt);
        bilinear_stiff->CalcElementMatrix(elemmat, elemIt, elemIt);
        
        help = elemmat * magvecpot;
        energy += 0.5 * (help * magvecpot);
        
      }  
      actVal[regionIt.GetPos()] = energy;
      delete bilinear_stiff;    
    }


  }

  // *************
  //   ComputeUI
  // *************
  template<class TYPE>
  void MagPDE::CalcFlux( shared_ptr<BaseResult> result,
                         bool timeDeriv ) {

    ENTER_FCN( "MagPDE::ComputeUI" );

    Warning( "Do implementation! -> Also for the coils and UI!",
             __FILE__, __LINE__ );

//     uiSD.Resize(coilDef_.GetSize());
//     uiSD.Init();
  
//     // loop over all subdomains
//     for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {

//       for (UInt dom=0; dom<coilDef_.GetSize(); dom++) {
//         if (subdoms_[actSD] == coilRegionId_[dom]) {
           
//           ElemList actSDList(ptgrid_ );
//           actSDList.SetRegion( subdoms_[actSD] );
          
//           EntityIterator it = actSDList.GetIterator();
//           for ( it.Begin(); !it.IsEnd(); it++ ) {
//             BaseFE * ptEl = it.GetElem()->ptElem;
                
//             const UInt nrIntPts= ptEl->GetNumIntPoints();
//             const Vector<Double> & intWeights = ptEl->GetIntWeights();  
//             Double jacDet;
                
//             StdVector<UInt> connect;
//             connect = it.GetElem()->connect;

//             Matrix<Double> ptCoord;
//             ptgrid_->GetElemNodesCoord( ptCoord, connect, true );

//             Vector<Double> magVecDeriv1Elem;
//             if( timeDeriv ) {
//               GetDerivSolVecOfElement(magVecDeriv1Elem,it);
//             } else {
//               GetSolVecOfElement(magVecDeriv1Elem,it);              
//             }
//             Double uiElem=0;
                
//             for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {
//               Vector<Double> shapeFnc;
//               jacDet = ptEl->CalcJacobianDetAtIp(actIntPt, ptCoord, 
//                                                  it.GetElem());    
//               ptEl -> GetShFncAtIp(shapeFnc, actIntPt, it.GetElem() );

//               uiElem += shapeFnc * magVecDeriv1Elem;
                    
//               if (isaxi_) {
//                 Vector<Double> coordAtIP;
//                 coordAtIP = ptCoord * shapeFnc;
//                 uiElem += shapeFnc * magVecDeriv1Elem * 2 * PI * coordAtIP[0]
//                   * intWeights[actIntPt-1] * jacDet;
//               }
//               else {
//                 uiElem += shapeFnc * magVecDeriv1Elem 
//                   * intWeights[actIntPt-1] * jacDet;
//               }
//             }
            
//             uiSD[dom] += uiElem;
//           }    
//         }
//       }
//     }
  }


  // ***************************************************
  //   Store currents/voltages and inductivity in file
  // ***************************************************
  void MagPDE::WriteUI2File(Vector<Double>& uiSD) {

    ENTER_FCN( "MagPDE::WriteUI2File" );

    Vector<UInt> coilIDs;   // just positive ids
    coilIDs.Push_back(abs(coilDef_[0]->id_));

    Integer maxID = coilDef_[0]->id_;

    for (UInt dom=1; dom < coilDef_.GetSize(); dom++) {

      bool isInVec = false;
      
      for (UInt dom2=0; dom2 < coilIDs.GetSize(); dom2++)    
        if (abs(coilDef_[dom]->id_) == coilIDs[dom2])
          isInVec = true;
      
      if (!isInVec)
        coilIDs.Push_back(abs(coilDef_[dom]->id_));
      
      if (abs(coilDef_[dom]->id_) > maxID)
        maxID = coilDef_[dom]->id_;
    }

    Vector<Double> uiID(maxID);
    uiID.Init();

    for (UInt dom=0; dom < coilDef_.GetSize(); dom++) {
      Integer actCoilID = coilDef_[dom]->id_;
      uiID[abs(actCoilID)-1] += uiSD[dom] * actCoilID/abs(actCoilID);
    }

    *UIfile_ << solveStep_->GetActTime() << " \t";
    for (UInt actID=0; actID < coilIDs.GetSize(); actID++) {
      if ( coilDef_[coilIDs[actID]-1]->coilType_ == Coil::MEASUREMENT2D )
        *UIfile_ << uiID[coilIDs[actID]-1] *
          coilDef_[coilIDs[actID]-1]->windingCrossSection_ << " \t";
    }
  
    *UIfile_ << myEndl;
  
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagPDE::ReadCoils() {

    ENTER_FCN( "MagEdgePDE::ReadCoils" );

    // Check if the element "coils" is present at all.
    // Otherwise leave
    ParamNode * coilNode = myParam_->Get( "coils", false );
    if ( !coilNode )
      return;
    
    // Get single coil nodes
    StdVector<ParamNode*> coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    if( coilNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Using the following coils:\n" );
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {
        
        // get region name of actual coil
        std::string regionName = coilNodes[i]->Get("name")->AsString();
        RegionIdType regionId = ptgrid_->RegionNameToId( regionName );

        coilRegionId_.Push_back( regionId );
        coilDef_.Push_back( shared_ptr<Coil>( new Coil( regionId,
                                                        coilNodes[i], ptgrid_) ) );
        Info->PrintCoil( *coilDef_.Last(), analysistype_ );
      }
    }
  }


  // ********************************************************
  //   Query parameter object for information about magnets
  // ********************************************************
  void MagPDE::ReadMagnets() {

    ENTER_FCN( "MagEdgePDE::ReadMagnets" );


    // Check if the element "magnets" is present at all.
    // Otherwise leave
    ParamNode * magnetNode = myParam_->Get( "magnets", false );
    if ( !magnetNode )
      return;

    // Get single magnet nodes
    StdVector<ParamNode*> magnetNodes = magnetNode->GetChildren();

    // trigger definition of magnets
    if( magnetNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_,
              "Found permanent magnets in the following regions:\n" );
      
      Double tmpDir;
      for( UInt i = 0; i < magnetNodes.GetSize(); i++ ) {
        
        // get region name of actual magnet
        std::string regionName = magnetNodes[i]->Get("name")->AsString();
        RegionIdType regionId = ptgrid_->RegionNameToId( regionName );
        
        magnetsDomain_.Push_back( regionId );
        
        // read orientation
        magnetNodes[i]->Get( "orientX", tmpDir );
        magnetsOriX_.Push_back( tmpDir );
          
        magnetNodes[i]->Get( "orientY", tmpDir );
        magnetsOriY_.Push_back( tmpDir );
          
        magnetNodes[i]->Get( "orientZ", tmpDir );
        magnetsOriZ_.Push_back( tmpDir );

        // report name to logfile
        Info->PrintF( pdename_, " %s\n", regionName.c_str());
      }
      Info->PrintF( "", "\n" );
    }
  }


  void MagPDE::DefineAvailResults() {
    ENTER_FCN( "MagPDE::DefineAvailResults" );
    
    StdVector<std::string> vecComponents;
    if( is3d_ ) {
      vecComponents = "x", "y", "z";
    }
    else if( isaxi_ ) {
      vecComponents = "r", "z";
    } 
    else {
      vecComponents = "x", "y";
    }
    
    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = MAG_POTENTIAL;
    if ( is3d_ ) {
      res1->dofNames = vecComponents;
    }
    else {
      res1->dofNames = "";
    }
    res1->unit = "Vs/m";
    res1->definedOn = ResultInfo::NODE;
    if ( is3d_ ) {
      res1->entryType = ResultInfo::VECTOR;
    }
    else {
      res1->entryType = ResultInfo::SCALAR;
    }
    res1->fctType = fct;
    
    results_.Push_back( res1 );
    availResults_.insert( res1 );

    if (is3d_ & analysistype_ != STATIC ) {
      // === MAGNETIC SCALAR POTENTIAL ===
      shared_ptr<ResultInfo> res2(new ResultInfo);
      res2->resultType = ELEC_POTENTIAL;
      res2->dofNames = "";
      res2->unit = "V";
      res2->definedOn = ResultInfo::NODE;
      res2->entryType = ResultInfo::SCALAR;
      res2->fctType = res1->fctType;
      results_.Push_back( res2 );
      availResults_.insert( res2 );
    }

    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX_DENSITY;
    flux->dofNames = vecComponents;
    flux->unit = "Vs/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    flux->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( flux ); 
    
    // === EDDY CURRENT DENSITY ===
    shared_ptr<ResultInfo> eddy(new ResultInfo);
    eddy->resultType = MAG_EDDY_CURRENT;
    if ( is3d_ ) {
      eddy->dofNames = vecComponents;
    }
    else {
      eddy->dofNames = "";
    }
    eddy->unit = "A/m^2";
    eddy->definedOn = ResultInfo::ELEMENT;
    eddy->entryType = ResultInfo::VECTOR;
    eddy->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( eddy ); 

    // === EDDY CURRENT POWER ===
    shared_ptr<ResultInfo> eddyPower(new ResultInfo);
    eddyPower->resultType = MAG_EDDY_POWER;
    eddyPower->dofNames = "";
    eddyPower->unit = "W";
    eddyPower->definedOn = ResultInfo::REGION;
    eddyPower->entryType = ResultInfo::SCALAR;
    eddyPower->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( eddyPower );            

    // === MAGNETIC ENERGY ===
    shared_ptr<ResultInfo> energy(new ResultInfo);
    energy->resultType = MAG_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( energy );            
    
  }

  void MagPDE::ReadSpecialResults() {
    ENTER_FCN( "MagPDE::ReadSpecialResults" );
    
    // There is a small problem with the magnetic force VWP:
    // The force itself is primarily calculated on nodes,
    // which makes it dependent on the discretization.

    // Thereofore we sould choose between two approaches:
    // a) Either we calculate the force-density or
    // b) We calculate only the sum of the nodes. However
    //    since we have no possibility to associate a 
    //    result on "nameNodesName", we should consider
    //    using a surface mesh for the determination of the
    //    virtual work principle
    


  //   StdVector<std::string> vecComponents;
//     if( isaxi_ ) {
//       vecComponents = "r", "z";
//     } else {
//       vecComponents = "x", "y";
//     }

//     StdVector<std::string> keyVec, attrVec, valVec;
    
//     // === MAGNETIC_FORCE (Virtual Work Principle, VWP ) ===
//     StdVector<std::string> forceNodeNames, forceRegions;
//     StdVector<UInt> forceNodes;
//     StdVector<RegionIdType> forceRegionIds;
//     std::string quantity;
//     Enum2String(MAG_FORCE_VWP, quantity);

//     // try to find nodes 
//     keyVec  = pdename_, "storeResults", "nodeResult", "nodes";
//     attrVec = "", "", "type";
//     valVec = "", "",quantity;
//     params->GetList( keyVec, attrVec, valVec, forceNodeNames);

//     keyVec = pdename_, "storeResults", "nodeResults", "region";
//     params->GetList( keyVec, attrVec, valVec, forceRegions);
    
//     if ( forceNodeNames.GetSize() > 0 ) {

//       // Force only computable for static/transient case
//       if( isComplex_) {
//         EXCEPTION( "Magnetic Force only computable for static/transient "
//                    << "case!" );
//       }

//       UInt numNodes = 0;
//       shared_ptr<ResultInfo> force(new ResultInfo);
//       force->resultType = MAG_FORCE_VWP;;
//       force->definedOn = ResultInfo::NODE;
//       force->entryType = ResultInfo::VECTOR;
//       force->dofNames = vecComponents;
//       force->unit = "N";


//       // count complete number of nodes
//       for ( UInt i=0; i<forceNodeNames.GetSize(); i++ ) {
//         numNodes+= ptgrid_->GetNumNodes(forceNodeNames[i]);
//       }
//       forceNodes.Resize(numNodes);
//       forceNodes.Init();
      
//       UInt inode =0;
//       StdVector<UInt> nodesConverted;
      
//       Info->PrintF( pdename_, " Computing '%s' for nodes:\n", quantity.c_str() );   
//       shared_ptr<EntityList> actList;
//       // get for each nodeslist all nodes
//       for (UInt i=0; i<forceNodeNames.GetSize(); i++) {
        
//         Info->PrintF( pdename_, " %s\n", forceNodeNames[i].c_str() );
        
//         // Get pure vector of nodes for force operator
//         ptgrid_->GetNodesByName( nodesConverted, forceNodeNames[i]);
        
//         for ( UInt j=0; j<nodesConverted.GetSize(); j++, inode++)
//           forceNodes[inode] = nodesConverted[j];
        
        
//         // Create Result object for postprocessing
//         actList = ptgrid_->GetEntityList( EntityList::NODE_LIST,
//                                           forceNodeNames[i], 
//                                           EntityList::NAMED_NODES );
//         shared_ptr<BaseResult> actSol;
//         actSol = shared_ptr<BaseResult>(new Result<Double>);
//         actSol->SetResultInfo( force );
//         actSol->SetEntityList( actList );
//         resultLists_[force].Push_back( actSol );
//       }
      
//       //initialize the force operator
//       NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);
//       ForceOpVWP_ = new  MagForceOp(ptgrid_, this, eqnMap_, *solhelp, dim_, 
//                                     materials_,  isaxi_, true );
      
//       ptgrid_->RegionNameToId( forceRegionIds, forceRegions ); 
//       ForceOpVWP_->Setup( forceRegionIds, forceNodes );
//     }
  }
    
    
  // =======================================================================
  //   HELPER METHODS FOR CALCULATING AUXILIARY QUANTITIES 
  // =======================================================================
  
  template<class TYPE>
  void MagPDE::CalcFluxDensityAtIP( EntityIterator it,
                                    UInt ip,
                                    Vector<TYPE>& field ) {
    ENTER_FCN( "MagPDE::CalcFluxDensityAtIP" );

    // get element solution
    Vector<TYPE> elSol;
    sol_->GetElemSolution(elSol, it);
    
    if( is3d_ ) {
      field.Resize(3);
      RegionIdType actRegionId = it.GetElem()->regionId;

      CurlCurlNode3DInt* curlOp = NULL;
      if ( regionNonLinType_[actRegionId] != NO_NONLINEARITY ) {
        std::string bilinearName = "nLinCurlCurlNode3DInt";
        curlOp = dynamic_cast<nLinCurlCurlNode3DInt*>(pdeBilinearForms_[actRegionId][bilinearName]);
      }
      else {
        std::string bilinearName = "CurlCurlNode3DInt";
        curlOp = dynamic_cast<CurlCurlNode3DInt*>(pdeBilinearForms_[actRegionId][bilinearName]);
      }
      
      //set element info
      curlOp->ExtractElemInfo( it );
      Matrix<Double> CornerCoords, bMatCurl, bMatDiv;
      ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect, 
                                  true );
      // case 1: element midpoint
      if( ip == 0 ) {
        Vector<Double> intPoint;
        it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
        curlOp->SetIntPoint(intPoint);
        curlOp->calcBMat(bMatCurl, bMatDiv, 0, CornerCoords);
        curlOp->UnsetIntPoint();     
      } else {
        // case2: real integration point
        curlOp->calcBMat(bMatCurl, bMatDiv, ip, CornerCoords);
      }
      field = bMatCurl * elSol;
      
      
    } else {
      field.Resize(2);
      RegionIdType actRegionId = it.GetElem()->regionId;
      
      CurlCurlNode2DInt* curlOp = NULL;
      
      if ( regionNonLinType_[actRegionId] != NO_NONLINEARITY ) {
        std::string bilinearName = "nLinCurlCurlNode2DInt";
        curlOp = dynamic_cast<nLinCurlCurlNode2DInt*>(pdeBilinearForms_[actRegionId][bilinearName]);
      } else {
        curlOp = dynamic_cast<CurlCurlNode2DInt*>
          (pdeBilinearForms_[actRegionId]["CurlCurlNode2DInt"]);
      }
      
      //set element info
      curlOp->ExtractElemInfo( it );
      Matrix<Double> CornerCoords, bMat;
      ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect, 
                                  true );
      // case 1: element midpoint
      if( ip == 0 ) {
        Vector<Double> intPoint;
        it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
        curlOp->SetIntPoint(intPoint);
        curlOp->calcBMat(bMat, 0, CornerCoords);
        curlOp->UnsetIntPoint();     
      } else {
        // case2: real integration point
        curlOp->calcBMat(bMat, ip, CornerCoords);
      }      
      field = bMat * elSol;

      // Account for curl 
      TYPE temp = field[0];
      if (isaxi_) {
        field[0] = -field[1];
        field[1] = temp;
      } else {
        field[0] = field[1];
        field[1] = -temp;
      }
    }
  }
  
  template<class TYPE>
  void MagPDE::CalcEddyCurrentAtIP( EntityIterator it,
                                    UInt ip,
                                    Vector<TYPE>& jEddy ) {
    ENTER_FCN( "MagPDE::CalcEddyCurrentAtIP" );


    Vector<Double> lCoord, shpFnc;
    Vector<TYPE> magVecDeriv1Elem, elecPotElem;
    Double conductivity = 0.0;
    
     // Get the right material parameter for current element
    RegionIdType actRegionId = it.GetElem()->regionId;
    materials_[actRegionId]
      ->GetScalar(conductivity,MAG_CONDUCTIVITY,REAL);
    BaseFE * ptEl = it.GetElem()->ptElem;
    
    
    if (is3d_ ) {
      jEddy.Resize(3);
    
      // Get electric bilinear form for regions with non-zero
      // conductivity
      linElecInt* elecBiLin = 
        dynamic_cast<linElecInt*>
        (pdeBilinearForms_[actRegionId]["linElecInt"]);
      
         // case 1: dummy integration point
      if( ip == 0 ) {
        it.GetElem()->ptElem->GetCoordMidPoint( lCoord );
        ptEl->GetShFnc(shpFnc,lCoord,it.GetElem());
      } else {
        ptEl->GetShFncAtIp(shpFnc,ip,it.GetElem());
      }
         
      // 1) get part coming from vector potential
      GetDerivSolVecOfElement(magVecDeriv1Elem,it,results_[0]);
      jEddy.Init();
      for( UInt iDof = 0; iDof < 3; iDof++ ) {
        for( UInt i = 0; i < shpFnc.GetSize(); i++ ) {
          jEddy[iDof] += shpFnc[i] * magVecDeriv1Elem[i*3+iDof];
        }
      }
      
      // 2) get part from scalar vector potential (only if region is conductive)
      if ( conductivity > EPS ) {
        GetDerivSolVecOfElement(elecPotElem,it,results_[1]);
        Matrix<Double> bMat;
        if( ip == 0 ) {
          elecBiLin->BDBInt::calcBMat( it, bMat);
        } else {
          Matrix<Double> ptCoord;
          ptgrid_->GetElemNodesCoord( ptCoord, it.GetElem()->connect, true );
          elecBiLin->ExtractElemInfo( it );
          elecBiLin->calcBMat( bMat, ip, ptCoord );
        }
        jEddy += bMat*elecPotElem;
      }
      jEddy *= -conductivity;
      
    } else {
      // **** 2D case ****
      jEddy.Resize(1);
      TYPE jEddyTemp;
      if( ip == 0 ) {
        it.GetElem()->ptElem->GetCoordMidPoint( lCoord );
        ptEl->GetShFnc(shpFnc,lCoord,it.GetElem());
      } else { 
        ptEl->GetShFncAtIp(shpFnc,ip,it.GetElem());
      }
      GetDerivSolVecOfElement(magVecDeriv1Elem,it,results_[0]);
      
      jEddyTemp = magVecDeriv1Elem * shpFnc;
      jEddyTemp *= -conductivity;
      jEddy[0] = jEddyTemp;
    }
  }
    

  


  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MagPDE::InitCoupling(PDECoupling * Coupling) {

    ENTER_FCN( "MagPDE::InitCoupling" );
  
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;

    // Enable update of geometry
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();  

    StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
    elemNodeToCouplingNode_.Resize(numCouplings);
    elemNodeToCouplingNode_.Init();


    for ( UInt actCoupling = 0; actCoupling < numCouplings; actCoupling++ ){

      if (ptCoupling_->GetOutputQuantity(actCoupling) == MAG_FORCE_LORENTZ) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

        //get the element-node to coupling node matching
        StdVector<std::string> couplRegions;
        StdVector<RegionIdType> regionIds;
        ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
        ptgrid_->RegionNameToId( regionIds, couplRegions );

        //Get total number of coupling elements
        UInt totalCouplingElems = ptgrid_->GetNumElems( regionIds );
        
        elemNodeToCouplingNode_tmp.Resize(totalCouplingElems);
        elemNodeToCouplingNode_tmp.Init();

        UInt offset = 0;
        for ( UInt reg = 0; reg < couplRegions.GetSize(); reg++ ) {

          // find subdomain index
          Integer SDidx=-1; for (UInt sd=0; sd<subdoms_.GetSize(); sd++) {
            if (regionIds[reg] == subdoms_[sd]) {
              SDidx = sd;
              break;
            }
          }
              
          if (SDidx==-1) {
            EXCEPTION( "magneticPDE: Coupling Region is not within the "
                       << "subdomains of the PDE" );
          }

          StdVector<Elem*> elemssd;
          ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);

          StdVector<UInt> * couplingnodes = NULL;
          ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
          if ( couplingnodes == NULL ) {
            EXCEPTION( "Pointer couplingnodes = NULL!" );
          }

          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            StdVector<UInt> & connecth = elemssd[actEl]->connect;
            elemNodeToCouplingNode_tmp[offset+actEl].Resize(connecth.GetSize());
            elemNodeToCouplingNode_tmp[offset+actEl].Init();

            for ( UInt ielemnode = 0; ielemnode < connecth.GetSize();
                  ielemnode++ ) {
              for ( UInt cnode = 0; cnode < (*couplingnodes).GetSize();
                    cnode++ ) {
                if (connecth[ielemnode] == (*couplingnodes)[cnode] ) {
                  elemNodeToCouplingNode_tmp[offset+actEl][ielemnode] = cnode;
                  break;
                }
              }
            }
          }

          //in the case, that we have more than one coupling region!
          offset = elemssd.GetSize();
        }

        elemNodeToCouplingNode_[actCoupling]  = elemNodeToCouplingNode_tmp;
      }
    }
  }

  void MagPDE::CalcOutputCoupling() {

    ENTER_FCN( "MagPDE::CalcOutputCoupling" );

    SolutionType quantity;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * values = NULL;
    UInt forcesCount = 0;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for ( UInt actCoupling = 0;
          actCoupling < ptCoupling_->GetNumOutputCouplings();
          actCoupling++ ) {

      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      ptCoupling_->GetOutputValues(actCoupling, values);

      Vector<Double> *temp = dynamic_cast<Vector<Double> *>(values);
      
      switch(ptCoupling_->GetOutputType(actCoupling)) {
          
      case NODE:
        ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
          
        if (quantity == MAG_FORCE_LORENTZ) {
          CalcNodeForceLorentz(*temp, elemNodeToCouplingNode_[forcesCount],
                               actCoupling, couplingNodes->GetSize());
          forcesCount++;
        }

        break;
          
      case ELEM:
        EXCEPTION( "No Element coupling output" );
      }
    }
  }

  void MagPDE::
  CalcNodeForceLorentz(Vector<Double> & force, 
                       StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                       UInt actCoupling, UInt numCouplingNodes) {

    ENTER_FCN( "MagPDE::CalcNodeForceLorentz" );

    NodeStoreSol<Double> *solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);

    //get the coupling regions
    StdVector<std::string> couplRegions;
    StdVector<RegionIdType> regionIds;
    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
    ptgrid_->RegionNameToId( regionIds, couplRegions );
    MagLorentzForceOp<Double> *ForceOp; 

    ForceOp  = new MagLorentzForceOp<Double>( ptgrid_, this,
                                              eqnMap_, *solhelp, isaxi_, true  );

    force.Init(0.0);

    Vector<Double> Jeddy;

    UInt offset = 0;
    for (UInt reg=0; reg<couplRegions.GetSize(); reg++) {

      //find subdomain index
      Integer SDidx=-1; for (UInt sd=0; sd<subdoms_.GetSize(); sd++) {
        if (regionIds[reg] == subdoms_[sd]) {
          SDidx = sd;
          break;
        }
      }

      Double conductivity;
      materials_[subdoms_[SDidx]]->GetScalar(conductivity,MAG_CONDUCTIVITY,REAL);      
      
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion(subdoms_[SDidx] );
      
      EntityIterator it = actSDList.GetIterator();
      UInt actEl = 0;
      for ( it.Begin(); !it.IsEnd(); it++, actEl++ ) {
        
        StdVector<UInt> const & connecth = it.GetElem()->connect;
        Matrix<Double> elemForce;
        GetDerivSolVecOfElement(Jeddy,it,results_[0]);
        Jeddy *= -conductivity;

        ForceOp->CalcElemMagLorentzForce(elemForce, Jeddy, it);

        // Add the element force to the according coupling node
        for (UInt ielemnode=0; ielemnode<connecth.GetSize(); ielemnode++) {
          for( UInt idim=0; idim<dim_; idim++) {
            force[elemNodeToCouplingNode[actEl+offset][ielemnode]*dim_+idim]
              += elemForce[ielemnode][idim];
          }
        }
      }

      //in the case, that we have more than one coupling region!
      offset = actSDList.GetSize();
    }
  }

  bool MagPDE::HasOutput( SolutionType output ) {
    ENTER_FCN( "MagPDE::HasOutput" );
    if (output == MAG_FORCE_LORENTZ) {
      return true;
    }
    return false;
  }

} // end of namespace

