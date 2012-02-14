// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <utility>

#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Driver/solveStepMagHyst.hh"
#include "Driver/stdSolveStep.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "Forms/bdbInt.hh"
#include "Forms/curlCurlNodeInt.hh"
#include "Forms/laplaceInt.hh"
#include "Forms/linGradBDBInt.hh"
#include "Forms/linearForm.hh"
#include "Forms/magforceop.hh"
#include "Forms/massInt.hh"
#include "Forms/nLinMagHystInt2D.hh"
#include "Forms/nLincurlCurlNodeInt.hh"
#include "Forms/nonConformingInt.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Materials/baseMaterial.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/basePDE.hh"
#include "PDE/eqnMap.hh"
#include "Utils/ApproxData.hh"
#include "Utils/Coil.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/biotSavart.hh"
#include "Utils/coordSystem.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/result.hh"
#include "Utils/tools.hh"
#include "magneticPDE.hh"
#include "trapezoidal.hh"

namespace CoupledField {

DECLARE_LOG(magpde)
DEFINE_LOG(magpde, "magpde")

  // **************
  //  Constructor
  // **************
  MagPDE::MagPDE( Grid * aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magnetic";
    pdematerialclass_ = ELECTROMAGNETIC;
    maxTimeDerivOrder_ = 1;

    // check if we have a 3d setup
    is3d_ = false;
    is3d_ = param->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
  }


  // *************
  //  Destructor
  // *************
  MagPDE::~MagPDE() {
  }


  //****************
  // Initialize PDE
  //****************
  void MagPDE::Init(UInt sequenceStep, PtrParamNode base) {
    SinglePDE::Init(sequenceStep);
    
    // store regions on which the Lorentz force should be calculated
    PtrParamNode nodeStoreRes = myParam_->Get("storeResults", ParamNode::PASS);
    if (nodeStoreRes) {
      ParamNodeList resForceL
        = nodeStoreRes->GetListByVal("nodeResult", "type", "magForceLorentz");
      for (UInt i=0, n=resForceL.GetSize(); i<n; ++i) {
        if (resForceL[i]->Has("allRegions")) {
          regionsForceL_.insert(subdoms_.Begin(), subdoms_.End());
        } else {
          PtrParamNode nodeRegionList = resForceL[i]->Get("regionList", ParamNode::PASS);
          if (nodeRegionList) {
            ParamNodeList resultRegions
              = nodeRegionList->GetList("region");
            for (UInt iReg=0, nReg=resultRegions.GetSize(); iReg<nReg; ++iReg)
            {
              RegionIdType curRegId = 
                  ptgrid_->GetRegion().Parse(resultRegions[iReg]->Get("name")->As<std::string>());
              if (subdoms_.Find(curRegId) >= -1)
                regionsForceL_.insert(curRegId);
            }
          }
        }
      }
    }
  }
  
  
  // *********************************
  //  Read special boundary conitions
  // *********************************

  void MagPDE::ReadSpecialBCs() {

    // -----------------------------------
    //  Check for Biot-Savart formulation
    // -----------------------------------
    PtrParamNode bsNode = myParam_->Get("biotSavart", ParamNode::PASS );
    if( bsNode ) {
      
      // check, if we have transient simulation
      if( analysistype_ == BasePDE::TRANSIENT ) {
        EXCEPTION( "Biot-Savart is currently just implemented for static /"
                   << "harmonic simulations!" );
      }
      biotSavart_ = 
          shared_ptr<BiotSavart>(new BiotSavart());
      biotSavart_->SetFormulation(BiotSavart::VEC_POT);
      isBiotSavart_ = true;
    }
    
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

    SinglePDE::InitNonLin();
    
    //now do PDE specifics
    std::map<std::string, NonLinType>::iterator it;
    for ( it=nonLinTypes_.begin() ; it != nonLinTypes_.end(); it++ ) {
      if ( (*it).second == HYSTERESIS ) {
        nonLin_ = false;
        isHysteresis_ = true;
      }
      else if ( (*it).second == PERMEABILITY ) {
        nonLin_ = true;
        isHysteresis_ = false;
      }
    }
    
  }


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagPDE::DefineIntegrators() {


    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // flag for updatedLagrange formulation
    bool upLagrangeForm = true;

   // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Set current region and material
       actRegion = it->first;
       actMat    = it->second;

      // Get current region node
      std::string regionName = ptgrid_->GetRegion().ToString( actRegion );

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      //get possible nonlinearities defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 
    
      if ( nonLinTypes.GetSize() > 0 ) {

        if ( nonLinTypes.Find(HYSTERESIS) != -1 ) {
          EXCEPTION("Magnetics with hysteresis currently not supported");
          
          // hysteresis modeling in this region
          StdVector<Elem*> elemssd;
          ptgrid_->GetElems(elemssd, actRegion);
          UInt numElSD =  elemssd.GetSize();

          //allocate for hystersis modeling; 
          std::cout << "Do vec Hyst!" << std::endl;
          actMat->InitVecHyst(numElSD, actSDList, dim_);

          BaseForm *curlcurlHyst; 
          curlcurlHyst = new nLinMagHystInt2D( actMat, isaxi_, upLagrangeForm );
          curlcurlHyst->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
          
          BiLinFormContext * stiffContext = 
            new BiLinFormContext( curlcurlHyst, STIFFNESS );
          stiffContext->SetPtPdes(this, this);   
          stiffContext->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );     
          assemble_->AddBiLinearForm( stiffContext);
          
          
          //save bilinearForm
          pdeBilinearForms_[actRegion][curlcurlHyst->GetName()] = curlcurlHyst;
        }
        else if (  nonLinTypes.Find(PERMEABILITY) != -1 ) {
          // informs material that approx./interpol. for permeability is needed
          actMat->NeedApproxMatCurve( MAG_PERMEABILITY );

          BaseForm *curlcurlNL; 
          if( is3d_ ) {
            curlcurlNL = new nLinCurlCurlNode3DInt( actMat, upLagrangeForm );
          } else {
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
          pdeBilinearForms_[actRegion][curlcurlNL->GetName()] = curlcurlNL;
          

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
        pdeBilinearForms_[actRegion][curlcurl->GetName()] = curlcurl;


        if ( nonLin_ == true && isHysteresis_ == false ) {
          // we need nonlinear RHS linearform for both linear and nonlinear
          // subdomains; just in case of material nonlinearity!
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
      materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);

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
        linGradBDBInt * elecBiLin = new linGradBDBInt( actMat, ELEC_PERMITTIVITY,
                                                 FULL, true );
        pdeBilinearForms_[actRegion][elecBiLin->GetName()] = elecBiLin;
      }


      // If this subdomain is a coil we have to do special things
      for ( UInt coil = 0; coil < coilDef_.GetSize(); coil++ ) {
        if ( actRegion == coilRegionId_[coil] ) {
          std::string factor = coilDef_[coil]->value_ + "/" +
            lexical_cast<std::string>(coilDef_[coil]->windingCrossSection_);

          if ( is3d_ ) {
          	VolForceInt *coilSource3d = 
              new VolForceInt( dim_, coilDef_[coil]->phase_,
              									isaxi_ );


            StdVector<std::string> currDensity(3);
            currDensity[0] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[0]);
            currDensity[1] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[1]);
            currDensity[2] = factor + "*" + lexical_cast<std::string>(coilDef_[coil]->locFlowDir_[2]);
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
        if ( actRegion == magnetsDomain_[perm] ) {

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
          actMat->GetScalar(reluctivity,MAG_RELUCTIVITY,Global::REAL);

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
    
    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    PtrParamNode ncIfaceListNode
        = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);
    
    // Get index of LAGRANGE_MULT result, just in case of coupled magnetics
    UInt lmResultIdx = 0;
    for(UInt i=0, n=results_.GetSize(); i<n; i++) {
      if(results_[i]->resultType == LAGRANGE_MULT) {
        lmResultIdx = i;
        break;
      }
    }
    LOG_DBG2(magpde) << "NonMatching: Index of LAGRANGE_MULT result: "
                     << lmResultIdx;
    
    for( UInt i = 0, n = ncIFaces_.GetSize(); i < n; i++ ) {
      // get regionId of Lagrangian surface
      StdVector<std::string> keyVec, attrVec, valVec;
      std::string slaveSide;
      std::string ncIfaceName = ptgrid_->GetRegion().ToString(ncIFaces_[i]);

      if (!ncIfaceListNode) {
        EXCEPTION("No ncInterfaces defined in domain section.");
      }
      PtrParamNode curNciNode = ncIfaceListNode->GetByVal("ncInterface", "name",
                                                   ncIfaceName);
      slaveSide = curNciNode->Get("slaveSide")->As<std::string>();

      // Part 1: Define integrator M(u, Lambda) on
      //         non-conforming interface (master/slave side)
      LOG_DBG2(magpde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
      actNcList->SetRegion( ncIFaces_[i] );

      UInt numDofs = 1;
      if (is3d_ ) {
        numDofs = 3;
      }
      
      NonConformingInt * ncInt =
        new NonConformingInt( numDofs, isaxi_ );
//      MassInt * ncInt = new MassInt( -1.0, dim_, isaxi_ );

      NcBiLinFormContext * stiffIntDescr =
        new NcBiLinFormContext( ncInt , STIFFNESS );

      // Force assembling of M(u, Lambda)^T
      stiffIntDescr->SetCounterPart( true );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[lmResultIdx],
                                 actNcList, actNcList );

      assemble_->AddBiLinearForm( stiffIntDescr );


      // Part 2: Define integrator D(u, Lambda) on
      //         Lagrangian surface (slave side)
      LOG_DBG2(magpde) << "NonMatching: Defining mass integrator"
                        << " for D on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );
      actSDList->SetRegion( ptgrid_->GetRegion().Parse( slaveSide ) );

      // D(u, Lambda) has the form of a standard mass
      // integrator with factor 1.0
      MassInt * dMatInt = new MassInt( 1.0, numDofs, isaxi_ );
      BiLinFormContext * dMatContext =
        new BiLinFormContext( dMatInt, STIFFNESS );

      // Force assembling of D(u, Lambda)^T
      dMatContext->SetCounterPart( true );
      dMatContext->SetPtPdes( this, this );
      dMatContext->SetResults( results_[0], results_[lmResultIdx],
                               actSDList, actSDList );

      assemble_->AddBiLinearForm( dMatContext );

      // Give result LAGRANGE_MULT to equation numbering class
      eqnMap_->AddResult( *results_[lmResultIdx], actSDList );
    }
    
  }

  void MagPDE::DefineSolveStep()
  {

    if ( isHysteresis_ ) 
      solveStep_ = new SolveStepMagHyst(*this);
    else 
      solveStep_ = new StdSolveStep(*this);
  }
  void MagPDE::PreparePDE4Computation()
  {
    if ( isBiotSavart_ ) {
      PtrParamNode bsNode = myParam_->Get("biotSavart", ParamNode::PASS );
      biotSavart_->Init( bsNode, ptgrid_, eqnMap_ );
      //biotSavart_->ComputeBiotSavartField();
    }
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagPDE::InitTimeStepping() {

    PtrParamNode systemNode = FindLinearSystem(pdename_);
    TS_alg_ = new Trapezoidal( algsys_, systemNode );
  }




  // ***************
  //   CalcResults
  // ***************
  void MagPDE::CalcResults ( shared_ptr<BaseResult> res ) {

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
    
    case MAG_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( res, results_[0] );
      } else {
        ExtractRhsResult<Double>( res, results_[0] );
      }
      break;

    case MAG_FLUX_DENSITY:
      if( isComplex_ ) {
        CalcFluxDensity<Complex>( res );
      } else {
        CalcFluxDensity<Double>( res );
      }
      break;

    case MAG_HFIELD:
      if( isComplex_ ) {
        EXCEPTION("H-field computation just for static and transient analysis");
      } else {
        CalcHfield( res );
      }
      break;

    case MAG_EDDY_CURRENT:
      if( isComplex_ ) {
        CalcEddyCurrent<Complex>( res );
      } else {
        CalcEddyCurrent<Double>( res );
      }
      break;
    
    case MAG_ELEM_PERMEABILITY:
      if( isComplex_ ) {
        CalcPermeability<Complex>( res );
      } else {
        CalcPermeability<Double>( res );
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

    case MAG_FORCE_LORENTZ:
      if (isComplex_) {
        CalcForceLorentz<Complex>( res );
      } else {
        CalcForceLorentz<Double>( res );
      }
      break;
      
    default:
      WARN( "Resulttype not computable by magnetic PDE" );
    }

    // In any case, we should trigger calculation of magnetic
    // quantities
    if( isComplex_ ) {
      WriteUI2File<Complex>();
    } else {
      WriteUI2File<Double>();
    }
  }

  

  void MagPDE::CalcForceVWP( shared_ptr<BaseResult> result ) {

    // get postprocessing result
     Result<Double> &  actRes = 
      dynamic_cast<Result<Double>&>(*result);

     EntityIterator regionIt = actRes.GetEntityList()->GetIterator();

    // resize vector
    Vector<Double> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() * dim_ );
    
    std::map<RegionIdType, MagForceOp*>::iterator itForce;
    itForce =  ForcePostVWP_.find(regionIt.GetRegion());

    if ( itForce == ForcePostVWP_.end() ) {
      //initialize the forve operator for the corresponding region

      //get all regions, the PDE is defined on
      StdVector<RegionIdType> allRegionIds(materials_.size());
      std::map<RegionIdType, BaseMaterial*>::iterator itMat;
      UInt i=0;
      for ( itMat = materials_.begin(); itMat != materials_.end(); itMat++ ) {
        allRegionIds[i] = itMat->first;
        i++;
      }

      UInt numNodes = ptgrid_->GetNumNodes(  regionIt.GetRegion() );
      StdVector<UInt> nodes;
      nodes.Reserve(numNodes);    
      ptgrid_->GetNodesByRegion(nodes, regionIt.GetRegion());
      
      NodeStoreSol<Double> * solhelp = 
        dynamic_cast<NodeStoreSol<Double> *>(sol_);
      MagForceOp* forceOp = new  MagForceOp(ptgrid_, this,  eqnMap_, 
                                            *solhelp, dim_, materials_, 
                                            isaxi_, true );
      
      //as neighbours we define all region Ids of the PDE
      forceOp->Setup(allRegionIds, nodes);
      ForcePostVWP_[regionIt.GetRegion()] = forceOp ;  
    }

    //get correct force operator         
    itForce =  ForcePostVWP_.find(regionIt.GetRegion());
    MagForceOp* forceOp = itForce->second;

    //compute force
    Vector<Double> totalForce, nodesForces;
    nodesForces.Resize(forceOp->GetNumCouplingNodes()*dim_);

    forceOp->CalcNodeForce(nodesForces, totalForce);
    for( UInt iDof = 0; iDof < dim_; iDof++ ) {
      actVal[ regionIt.GetPos()*dim_ + iDof] = totalForce[iDof];
    }

  }

  template<class TYPE>
  void MagPDE::CalcFluxDensity( shared_ptr<BaseResult> result ) {

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

  void MagPDE::CalcHfield( shared_ptr<BaseResult> result ) {

     Vector<Double> elemHfield(dim_);
    
    // fetch result object and convert data
    Result<Double> &  actSol = 
      dynamic_cast<Result<Double>&>(*result);
    Vector<Double> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );

    // loop over elements
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      CalcHfieldAtIP( it, 0, elemHfield );
      for( UInt iDof = 0; iDof < dim_; iDof++ ) {
        actVal[it.GetPos()*dim_ + iDof] = elemHfield[iDof];
      }
    }
  }
  
  template<class TYPE>
  void MagPDE::CalcEddyCurrent( shared_ptr<BaseResult> result ) {
    
    
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
        ->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
      
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
        ->GetScalar(reluctivity,MAG_RELUCTIVITY,Global::REAL);
      
      // Create stiffness integrator
      //get possible nonlinearities defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[regionIt.GetRegion()]; 

      if ( nonLinTypes.Find(PERMEABILITY) != -1 ) {
        
        if( is3d_ ) {
          bilinear_stiff = new nLinCurlCurlNode3DInt( materials_[regionIt.GetRegion()],
                                                      true );
        } else  {
          bilinear_stiff = new nLinCurlCurlNode2DInt( materials_[regionIt.GetRegion()],
                                                      isaxi_, true );
        }
        
        bilinear_stiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        
        // VERY IMPORTANT: Set nonlinear-method "fixpoint", as otherwise also
        // the Frechet part of the stiffness is calculated!
        bilinear_stiff->SetNonLinMethod( FIXEDPOINT );
      } else {
        if (is3d_ ) {
          bilinear_stiff = new CurlCurlNode3DInt( materials_[regionIt.GetRegion()] );
        } else {
          bilinear_stiff = new CurlCurlNode2DInt( materials_[regionIt.GetRegion()], isaxi_);
        }
        bilinear_stiff->SetAnsatzFct( results_[0]->fctType, results_[0]->fctType );
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
  void MagPDE::CalcFlux( shared_ptr<Coil> coil, 
                         TYPE& flux,
                         bool timeDeriv ) {


    Vector<Double> elemMidPoint, flowDir(3);

    // get math parser
    MathParser * parser = domain->GetMathParser();
    MathParser::HandleType h1, h2, h3;
    h1 = parser->GetNewHandle();
    h2 = parser->GetNewHandle();
    h3 = parser->GetNewHandle();

    // get current coil
    shared_ptr<Coil> actCoil = coil;
    flux = 0.0;
    
    // get related region
    ElemList actSDList(ptgrid_ );
    actSDList.SetRegion( actCoil->coilRegionId_ );
    
    // get coordinate system of coil
    CoordSystem * coilCosy;
    if( actCoil->flowCoordSys_ ) {
      coilCosy = actCoil->flowCoordSys_;
    } else {
      coilCosy = domain->GetCoordSystem("default");
    }
    
    // if we are in 3d, set flow coordinate directions at parser
    if( is3d_ ) {
      parser->SetValue( h1, coilCosy->GetDofName(1), actCoil->locFlowDir_[0] );
      parser->SetValue( h2, coilCosy->GetDofName(2), actCoil->locFlowDir_[1] );
      parser->SetValue( h3, coilCosy->GetDofName(3), actCoil->locFlowDir_[2] );
      parser->SetExpr(h1, coilCosy->GetDofName(1) );
      parser->SetExpr(h2, coilCosy->GetDofName(2) );
      parser->SetExpr(h3, coilCosy->GetDofName(3) );
    }
    
    // loop over elements of region
    EntityIterator it = actSDList.GetIterator();
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      BaseFE * ptEl = it.GetElem()->ptElem;
        
      const UInt nrIntPts= ptEl->GetNumIntPoints();
      const Vector<Double> & intWeights = ptEl->GetIntWeights();  
      Double jacDet;
      
      StdVector<UInt> connect;
      connect = it.GetElem()->connect;
      
      Matrix<Double> ptCoord;
      ptgrid_->GetElemNodesCoord( ptCoord, connect, true );
      
      Vector<TYPE> magVecDeriv1Elem, temp;
      if( timeDeriv ) {
        GetDerivSolVecOfElement(magVecDeriv1Elem,it,results_[0]);
      } else {
        GetSolVecOfElement(magVecDeriv1Elem,it,results_[0]);
      }
      
      TYPE uiElem=0;
      
      if( is3d_ ) {
        
        // calculate global midpoint and register coordinates with
        // mathParser
        ptgrid_->GetGlobalElemMidPoint( it.GetElem()->elemNum, 
                                        elemMidPoint );
        parser->SetCoordinates( h1, *coilCosy, elemMidPoint );
        parser->SetCoordinates( h2, *coilCosy, elemMidPoint );
        parser->SetCoordinates( h3, *coilCosy, elemMidPoint );
        
        // evaluate flux direction
        flowDir[0] = parser->Eval( h1 );
        flowDir[1] = parser->Eval( h2 );
        flowDir[2] = parser->Eval( h3 );
        // normalize flux
        flowDir = flowDir/flowDir.NormL2();
        
        // calculate scalar product of potential and unit 
        // vector in current direction A * e_j
        temp = magVecDeriv1Elem;
        magVecDeriv1Elem.Resize((UInt)magVecDeriv1Elem.GetSize()/3);
        magVecDeriv1Elem.Init();
        UInt j = 0;
        for( UInt i = 0; i < magVecDeriv1Elem.GetSize(); i=i+3, j++ ) {
          magVecDeriv1Elem[j] += temp[i] * flowDir[0];
          magVecDeriv1Elem[j] += temp[i] * flowDir[1];
          magVecDeriv1Elem[j] += temp[i] * flowDir[2];
        }
          
      } 
      
      for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {
        Vector<Double> shapeFnc;
        jacDet = ptEl->CalcJacobianDetAtIp(actIntPt, ptCoord, 
                                           it.GetElem());    
        ptEl -> GetShFncAtIp(shapeFnc, actIntPt, it.GetElem() );
        
        if (isaxi_) {
          Vector<Double> coordAtIP;
          coordAtIP = ptCoord * shapeFnc;
          uiElem += magVecDeriv1Elem * 
            shapeFnc  * (2.0 * PI * coordAtIP[0]
                         * intWeights[actIntPt-1] * jacDet);
        }
        else {
          uiElem += magVecDeriv1Elem * shapeFnc 
            * intWeights[actIntPt-1] * jacDet;
        }
      }
      
      flux += uiElem;
    } 
    parser->ReleaseHandle(h1);
    parser->ReleaseHandle(h2);
    parser->ReleaseHandle(h3);
  }
  template<class TYPE>
  void MagPDE::CalcPermeability( shared_ptr<BaseResult> result ) {
    
    TYPE elemPerm;
    Vector<TYPE> elemFlux;
    Double bAbs, reluct;
   
   // fetch result object and convert data
   Result<TYPE> &  actSol = 
     dynamic_cast<Result<TYPE>&>(*result);
   Vector<TYPE> & actVal = actSol.GetVector();
   actVal.Resize( actSol.GetEntityList()->GetSize());

   // loop over elements
   EntityIterator it = actSol.GetEntityList()->GetIterator();
   for ( it.Begin(); !it.IsEnd(); it++ ) {
     
     // Determine regionId of element
     const Elem & actEl = *(it.GetElem());
     RegionIdType actRegion = actEl.regionId; 

     //get possible nonlinearities defined in this region
     StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 
     
     // Check, if region is nonlinear
      if ( nonLinTypes.Find(PERMEABILITY) != -1 )  {
       
       // Obtain nonlinear approximation functional
       ApproxData * approx  = materials_[actRegion]->GetNonlinFnc(MAG_PERMEABILITY);
       
       // Calculate flux density in element midpoint
       CalcFluxDensityAtIP( it, 0, elemFlux );
       bAbs = elemFlux.NormL2();
       reluct = approx->EvaluateFuncNu(bAbs);
       elemPerm = 1.0 / reluct;
     } else {
       materials_[actRegion]->GetScalar( elemPerm, MAG_PERMEABILITY,Global::REAL); 
     }
     actVal[it.GetPos() ] = elemPerm;
   }
  }

  template<class TYPE>
  void MagPDE::CalcForceLorentz( shared_ptr<BaseResult> result ) {
    Result<TYPE> &actSol = dynamic_cast<Result<TYPE>&>(*result);      
    std::map<UInt, UInt> nodeNumPos;
    Vector<TYPE> &actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
    EntityIterator nodeIt = actSol.GetEntityList()->GetIterator();
    StdVector<RegionIdType> regionsVec;
    
    for (nodeIt.Begin(); !nodeIt.IsEnd(); nodeIt++) {
      nodeNumPos[nodeIt.GetNode()] = nodeIt.GetPos();
    }
    
    std::set<RegionIdType>::iterator regIt = regionsForceL_.begin(),
                                     regItEnd = regionsForceL_.end();
    for ( ; regIt != regItEnd; ++regIt) {
        regionsVec.Push_back(*regIt);
    }
    
    CalcNodeForceLorentz(actVal, regionsVec, nodeNumPos);
    
    Vector<TYPE> sum(dim_);
    sum.Init();
    for( UInt i = 0; i < nodeNumPos.size(); i++ ) {
      for( UInt j = 0; j < dim_; j++ ) {
        sum[j] += actVal[i*dim_+j];
      }
    }
    std::cerr << "sum for region '" << actSol.GetEntityList()->GetName() << " is " 
    << sum << std::endl;
  }
  
  
  // ***************************************************
  //   Store currents/voltages and inductivity in file
  // ***************************************************
  template<typename T>
  std::string TypeToString( const T& s ) {
    return boost::lexical_cast<std::string>(s);
  }
  
  template<>
  std::string TypeToString( const Complex& s) {
    return lexical_cast<std::string>(std::abs(s)) + std::string("\t") +
    lexical_cast<std::string>(
    (std::abs(s.imag()) > 1e-16) ?                   
                                  std::atan2(s.imag(),s.real() )*180/PI : 0.0);
  }
  
  template<class TYPE>
  void MagPDE::WriteUI2File( ) {


    // get current time / frequency step
    
    MathParser * parser = domain->GetMathParser();
    MathParser::HandleType mHandle =  parser->GetNewHandle();
    parser->SetExpr( mHandle, "step" );
    UInt actStep = (UInt) parser->Eval(mHandle);



    // iterate over all coils
    for( UInt iCoil = 0; iCoil < coilDef_.GetSize(); iCoil++ ) {

      shared_ptr<Coil> actCoil = coilDef_[iCoil];

      // check if results was already written out for this step
      if (actStep == (UInt) actCoil->lastSaveStep_ ) { 
          continue;
      } else {
        actCoil->lastSaveStep_ = (UInt) actStep;
      }


        // if coil needs to compuet U, trigger calculation
        if( analysistype_ == TRANSIENT ||
            analysistype_ == HARMONIC ) {
          TYPE voltage = 0.0;
          if( actCoil->fileU_ ) {
            CalcFlux<TYPE>( actCoil, voltage, true );
            std::ofstream * uOut = actCoil->fileU_;
            *uOut << solveStep_->GetActStep() << " \t";
            *uOut << TypeToString (voltage / actCoil->windingCrossSection_) << " \n";
          }
        }

      // if coil needs to comput L, trigger calculation
      TYPE induct = 0.0;
      if( actCoil->fileL_ ) {
        CalcFlux<TYPE>( actCoil, induct, false );
        parser->SetExpr( mHandle, actCoil->value_);
        induct /= (actCoil->windingCrossSection_ *  parser->Eval(mHandle));
        std::ofstream * lOut = actCoil->fileL_;
        *lOut << solveStep_->GetActStep() << " \t";
        *lOut << TypeToString( induct ) << " \n";

      }
    }
    parser->ReleaseHandle( mHandle );
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagPDE::ReadCoils() {


    // Check if the element "coils" is present at all.
    // Otherwise leave
    PtrParamNode coilNode = myParam_->Get( "coils", ParamNode::PASS );
    if ( !coilNode )
      return;
    
    // Get single coil nodes
    ParamNodeList coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    if( coilNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Using the following coils:\n" );
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {
        
        // get region name of actual coil
        std::string regionName = coilNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptgrid_->GetRegion().Parse( regionName );

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



    // Check if the element "magnets" is present at all.
    // Otherwise leave
    PtrParamNode magnetNode = myParam_->Get( "magnets", ParamNode::PASS );
    if ( !magnetNode )
      return;

    // Get single magnet nodes
    ParamNodeList magnetNodes = magnetNode->GetChildren();

    // trigger definition of magnets
    if( magnetNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_,
              "Found permanent magnets in the following regions:\n" );
      
      Double tmpDir;
      for( UInt i = 0; i < magnetNodes.GetSize(); i++ ) {
        
        // get region name of actual magnet
        std::string regionName = magnetNodes[i]->Get("name")->As<std::string>();
        RegionIdType regionId = ptgrid_->GetRegion().Parse( regionName );
        
        magnetsDomain_.Push_back( regionId );
        
        // read orientation
        magnetNodes[i]->GetValue( "orientX", tmpDir );
        magnetsOriX_.Push_back( tmpDir );
          
        magnetNodes[i]->GetValue( "orientY", tmpDir );
        magnetsOriY_.Push_back( tmpDir );
          
        magnetNodes[i]->GetValue( "orientZ", tmpDir );
        magnetsOriZ_.Push_back( tmpDir );

        // report name to logfile
        Info->PrintF( pdename_, " %s\n", regionName.c_str());
      }
      Info->PrintF( "", "\n" );
    }
  }


  void MagPDE::DefineAvailResults() {
    
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
    res1->resultType = MAG_POTENTIAL;
    
    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->As<std::string>();
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      res1->fctType = fct;
      res1->definedOn = ResultInfo::NODE;
    } else if (approxType == "legendre" ) {
      shared_ptr<LegendreFct> fct(new LegendreFct);
      UInt order =  myParam_->Get("order")->As<UInt>();
      fct->SetIsoOrder( order );
      res1->definedOn = ResultInfo::PFEM;
      res1->fctType = fct;
    } else {
      EXCEPTION("Approximation type' " << approxType << "' not known");
    }
    
    if ( is3d_ ) {
      res1->dofNames = vecComponents;
    }
    else {
      res1->dofNames = "";
    }
    res1->unit = "Vs/m";
    if ( is3d_ ) {
      res1->entryType = ResultInfo::VECTOR;
    }
    else {
      res1->entryType = ResultInfo::SCALAR;
    }
    results_.Push_back( res1 );
    availResults_.insert( res1 );

    if (is3d_ && analysistype_ != STATIC ) {
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

    // === MAGNETIC RHS LOAD ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD;
    rhs->dofNames = vecComponents;
    rhs->unit = "Am";
    rhs->definedOn = ResultInfo::NODE;
    rhs->entryType = ResultInfo::VECTOR;
    rhs->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( rhs ); 

    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX_DENSITY;
    flux->dofNames = vecComponents;
    flux->unit = "Vs/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    flux->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( flux ); 
    
    // === MAGNETIC Field Intensity ===
    shared_ptr<ResultInfo> hfield(new ResultInfo);
    hfield->resultType = MAG_HFIELD;
    hfield->dofNames = vecComponents;
    hfield->unit = "A/m";
    hfield->definedOn = ResultInfo::ELEMENT;
    hfield->entryType = ResultInfo::VECTOR;
    hfield->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( hfield ); 

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
    
    // === PERMEABILITY  ===
    shared_ptr<ResultInfo> perm(new ResultInfo);
    perm->resultType = MAG_ELEM_PERMEABILITY;
    perm->dofNames = "";
    perm->unit = "Vs/Am";
    perm->definedOn = ResultInfo::ELEMENT;
    perm->entryType = ResultInfo::SCALAR;
    perm->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( perm );

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
    
    // === LORENTZ FORCE ===
    shared_ptr<ResultInfo> forceLorentz(new ResultInfo);
    forceLorentz->resultType = MAG_FORCE_LORENTZ;
    forceLorentz->dofNames = vecComponents;
    forceLorentz->unit = "N";
    forceLorentz->definedOn = ResultInfo::NODE;
    forceLorentz->entryType = ResultInfo::VECTOR;
    forceLorentz->fctType = shared_ptr<ConstFct>(new ConstFct());
    availResults_.insert( forceLorentz );

    // === VWP FORCE ===
    shared_ptr<ResultInfo> forceVWP(new ResultInfo);
    forceVWP->resultType = MAG_FORCE_VWP;
    forceVWP->dofNames = vecComponents;
    forceVWP->unit = "N";
    forceVWP->definedOn = ResultInfo::REGION; //    forceVWP->definedOn = ResultInfo::NODE;
    forceVWP->entryType = ResultInfo::VECTOR;
    forceVWP->fctType = shared_ptr<ConstFct>(new ConstFct());
    availResults_.insert( forceVWP );

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;

    LOG_DBG2(magpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    if(domainNCIfaceListNode)
    {
      PtrParamNode ncInterfaceListNode = myParam_->Get("ncInterfaceList", ParamNode::PASS );
      ParamNodeList pdeNCIfaceNodes;

      if(ncInterfaceListNode)
      {
        pdeNCIfaceNodes = ncInterfaceListNode->GetList("ncInterface");

        for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
          std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->As<std::string>();
          std::string domainIfaceName;

          PtrParamNode domainIfaceNode = domainNCIfaceListNode->GetByVal("ncInterface",
              "name",
              pdeIfaceName,
              ParamNode::PASS);
          if(!domainIfaceNode)
          {
            LOG_DBG2(magpde) << "NonMatching: Nonconforming "
                              << "interface '" << ncIfaceNames[i]
                              << "' does not exist in domain.";

            EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
          }

          ncIfaceNamesForPDE.Push_back(pdeIfaceName);
        }
        ptgrid_->GetRegion().Parse(ncIfaceNamesForPDE, ncIfaceIds);

        for (UInt i = 0; i < ncIfaceIds.GetSize(); i++) {
          ncIFaces_.Push_back(ncIfaceIds[i]);
        }

        // In the case of the presence of non-conforming interfaces,
        // a second resultdof object has to be created, which describes the
        // Lagrange multiplier
        if( ncIFaces_.GetSize() > 0 ) {
          LOG_DBG2(magpde) << "NonMatching: Defining new ResultDof "
                            << "Lagrange Multiplier (LM).";

          LOG_DBG3(magpde) << "NonMatching: Lagrange Multiplier DOFs: ";
          StdVector<std::string> lmDofNames;
          for( UInt i=0, n=res1->dofNames.GetSize(); i<n; i++ ) {
            lmDofNames.Push_back("LM_" + res1->dofNames[i]);
              LOG_DBG3(magpde) << "NonMatching: " << lmDofNames[i];
          }

          shared_ptr<ResultInfo> lagr ( new ResultInfo );
          lagr->resultType = LAGRANGE_MULT;
          lagr->dofNames = lmDofNames;
          lagr->fctType = results_[0]->fctType;
          lagr->definedOn = results_[0]->definedOn;
          results_.Push_back( lagr );
        }
      }

    }

  }

  void MagPDE::ReadSpecialResults() {
    
    // There is a small problem with the magnetic force VWP

    // The force itself is primarily calculated on nodes,
    // which makes it dependent on the discretization.

    // Therefore we should choose between two approaches:
    // a) Either we calculate the force-density or
    // b) We calculate only the sum of the nodes. However
    //    since we have no possibility to associate a 
    //    result on "namedNodesName", we should consider
    //    using a surface mesh for the determination of the
    //    virtual work principle
    


//     StdVector<std::string> vecComponents;
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
      
//       ptgrid_->GetRegion().Parse( forceRegionIds, forceRegions );
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

    // get element solution
    Vector<TYPE> elSol;
    sol_->GetElemSolution(elSol, it);

    RegionIdType actRegionId = it.GetElem()->regionId;

    //get possible nonlinearities defined in this region
    StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegionId]; 
     
    if( is3d_ ) {
      field.Resize(3);
    
      CurlCurlNode3DInt* curlOp = NULL;
      if ( nonLinTypes.Find(PERMEABILITY) != -1 ) {
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

      Matrix<Double> CornerCoords, bMat;
      ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect,true);

//      if ( regionNonLinType_[actRegionId] == HYSTERESIS ) {
//         std::string bilinearName = "nLinMagHystInt2D";
//         nLinMagHystInt2D* curlOpHyst = 
//           dynamic_cast<nLinMagHystInt2D*>(pdeBilinearForms_[actRegionId][bilinearName]);
      
//         //set element info
//         curlOpHyst->ExtractElemInfo( it );
 
//         // case 1: element midpoint
//         if( ip == 0 ) {
//           Vector<Double> intPoint;
//           it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
//           curlOpHyst->SetIntPoint(intPoint);
//           curlOpHyst->calcBMat(bMat, 0, CornerCoords);
//           curlOpHyst->UnsetIntPoint();     
//         } else {
//           // case2: real integration point
//           curlOpHyst->calcBMat(bMat, ip, CornerCoords);
//         }      
//       }
//       else {


      if ( nonLinTypes.Find(HYSTERESIS) != -1 ) {  
        EXCEPTION("Magnetics with hysteresis currently not supported");
//         std::string bilinearName = "nLinMagHystInt2D";
//         nLinMagHystInt2D* curlHystOp;
//         curlHystOp = dynamic_cast<nLinMagHystInt2D*>(pdeBilinearForms_[actRegionId][bilinearName]);

//         //set element info      
//         curlHystOp->ExtractElemInfo( it );

//         // case 1: element midpoint
//         if( ip == 0 ) {
//           Vector<Double> intPoint;
//           it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
//           curlHystOp->SetIntPoint(intPoint);
//           curlHystOp->calcBMat(bMat, 0, CornerCoords);
//           curlHystOp->UnsetIntPoint();     
//         } else {
//           // case2: real integration point
//           curlHystOp->calcBMat(bMat, ip, CornerCoords);
//         }      
      }
      else {
        CurlCurlNode2DInt* curlOp = NULL;
        if ( nonLinTypes.Find(PERMEABILITY) != -1 ) {
          std::string bilinearName = "nLinCurlCurlNode2DInt";
          curlOp = dynamic_cast<nLinCurlCurlNode2DInt*>(pdeBilinearForms_[actRegionId][bilinearName]);
        } else {
          curlOp = dynamic_cast<CurlCurlNode2DInt*>
            (pdeBilinearForms_[actRegionId]["CurlCurlNode2DInt"]);
        }
        
        //set element info      
        curlOp->ExtractElemInfo( it );
        // case 1: element midpoint
        if( ip == 0 ) {
          Vector<Double> intPoint;
          it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
          curlOp->SetIntPoint(intPoint);
          curlOp->CalcBMat(bMat, 0, CornerCoords);
          curlOp->UnsetIntPoint();     
        } else {
          // case2: real integration point
          curlOp->CalcBMat(bMat, ip, CornerCoords);
        }      
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


  void MagPDE::CalcHfieldAtIP( EntityIterator it, UInt ip,
                               Vector<Double>& field ) {

    RegionIdType actRegion = it.GetElem()->regionId;
  
    //get possible nonlinearities defined in this region
    StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion]; 

    if (  nonLinTypes.Find(HYSTERESIS) != -1 ) {
      EXCEPTION("Magnetics with hysteresis currently not supported");
//       materials_[actRegion]->GetVectorXHystVal( it.GetElem()->elemNum,
//                                                field ); 
    }
    else if ( nonLinTypes.Find(PERMEABILITY) != -1 ) {
      EXCEPTION("CalcHfieldAtIP for nonlinear BH curve not implemented");
    }
    else {
      //first compute flux density
      CalcFluxDensityAtIP( it, ip, field );

      Double reluctivity;
      materials_[actRegion]->GetScalar( reluctivity, MAG_RELUCTIVITY, Global::REAL);

      field *= reluctivity;
    }
    
  }
  

  template<class TYPE>
  void MagPDE::CalcEddyCurrentAtIP( EntityIterator it,
                                    UInt ip,
                                    Vector<TYPE>& jEddy ) {


    Vector<Double> lCoord, shpFnc;
    Vector<TYPE> magVecDeriv1Elem, elecPotElem;
    Double conductivity = 0.0;
    
     // Get the right material parameter for current element
    RegionIdType actRegionId = it.GetElem()->regionId;
    materials_[actRegionId]
      ->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
    BaseFE * ptEl = it.GetElem()->ptElem;
    
    
    if (is3d_ ) {
      jEddy.Resize(3);
    
      // Get electric bilinear form for regions with non-zero
      // conductivity
      linGradBDBInt* elecBiLin = 
        dynamic_cast<linGradBDBInt*>
        (pdeBilinearForms_[actRegionId]["linGradBDBInt"]);
      
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
          elecBiLin->CalcBMat( bMat, ip, ptCoord );
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

  
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;

    // Enable update of geometry
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();  


    cplNodeNumPos_.Resize( numCouplings );

    for ( UInt actCoupling = 0; actCoupling < numCouplings; actCoupling++ ){

      if (ptCoupling_->GetOutputQuantity(actCoupling) == MAG_FORCE_LORENTZ) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

        //get the element-node to coupling node matching
        StdVector<std::string> couplRegions;
        StdVector<RegionIdType> regionIds;
        ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
        ptgrid_->GetRegion().Parse(couplRegions, regionIds);
        
        // Check, that every coupling region is part of
        // the magnetic pde itself
        for( UInt iReg = 0; iReg < couplRegions.GetSize(); iReg++ ) {
          if( subdoms_.Find(regionIds[iReg]) == -1 ) {
            EXCEPTION( "Coupling region '" << couplRegions[iReg]
                       << "' is not contained in regions for"
                       << " magnetic PDE" );
          }
        }
        
        StdVector<UInt> * couplingnodes = NULL;
        ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
        
        if ( couplingnodes == NULL ) {
          EXCEPTION( "Pointer couplingnodes = NULL!" );
        }
        
        for( UInt iNode = 0; iNode < couplingnodes->GetSize(); iNode++ ) {
          UInt actNode = (*couplingnodes)[iNode];
          cplNodeNumPos_[actCoupling][actNode] = iNode;
        }
      }
      else if (ptCoupling_->GetOutputQuantity(actCoupling) == MAG_FORCE_VWP ) {
      
        // Initialization of coupling helper arrays
        StdVector<UInt> * couplingnodes = NULL;
        StdVector<std::string> * nRegions, coupRegions;
        StdVector<RegionIdType> nRegionIds;

        ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
        if (couplingnodes == NULL)
          std::cerr << "Couplingnodes = 0!!!!" << std::endl;

        // get volume neighbours lying next to coupling nodes, because 
        // these volume elements have to be  moved 'virtually'
        NodeStoreSol<Double> * solhelp = 
          dynamic_cast<NodeStoreSol<Double> *>(sol_);
        ForceOpVWP_ = new  MagForceOp(ptgrid_, this,  eqnMap_, 
                                      *solhelp, dim_, materials_, 
                                      isaxi_, true );
        
        ptCoupling_->GetOutputNeighbourRegion(actCoupling, nRegions);
        ptgrid_->GetRegion().Parse(*nRegions, nRegionIds);
        ForceOpVWP_->Setup(nRegionIds, *couplingnodes);

        ptCoupling_->GetOutputRegions(actCoupling, coupRegions);
      
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

      }
    }
  }

  void MagPDE::CalcOutputCoupling() {


    SolutionType quantity;
    StdVector<UInt> * couplingNodes = NULL;
    SingleVector * values = NULL;
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
          //get the coupling regions
          StdVector<std::string> couplRegions;
          StdVector<RegionIdType> regionIds;
          ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
          ptgrid_->GetRegion().Parse( couplRegions, regionIds );

          CalcNodeForceLorentz(*temp, regionIds, cplNodeNumPos_[forcesCount]);
          Vector<Double> sum(dim_);
          sum.Init();
          for( UInt i = 0; i < cplNodeNumPos_[forcesCount].size(); i++ ) {
            for( UInt j = 0; j < dim_; j++ ) {
              sum[j] += (*temp)[i*dim_+j];
            }
          }
          Info->PrintF(pdename_, "Sum of magnetic force (Lorentz):\n");
          Info->PrintVec(sum);
          forcesCount++;
        }
        else  if (quantity ==  MAG_FORCE_VWP) {
          Vector<Double> totalForce;
          ForceOpVWP_->CalcNodeForce(*temp, totalForce);
          
          // write information in .info-file
          Info->PrintF(pdename_, "Sum of magnetic force (VWM):\n");
          Info->PrintVec(totalForce);
          forcesCount++;
        }

        break;
          
      case ELEM:
        EXCEPTION( "No Element coupling output" );
      }
    }
  }

  template<class TYPE> void MagPDE::
  CalcNodeForceLorentz( Vector<TYPE> & force,
                        const StdVector<RegionIdType> regionIds,
                        const std::map<UInt, UInt> & nodeNumPos ) {
       
    force.Init();
    Vector<TYPE> fluxAtIp(dim_), elemForce, fAtIp, jAtIp;
    Matrix<Double> ptCoord;
    std::map<UInt, UInt>::const_iterator itPosEnd = nodeNumPos.end();
    
    for (UInt reg=0, numRegs=regionIds.GetSize(); reg < numRegs; ++reg) {

      //find subdomain index
      Integer sdIndex = subdoms_.Find( regionIds[reg] );
      if( sdIndex == -1 ) {
        EXCEPTION( "The region coupling region '" <<
            ptgrid_->GetRegion().ToString( regionIds[reg] )
            << "' was not found in magneticPDE" );
      }

      // Check if this region is a coil
      RegionIdType actRegionId = subdoms_[sdIndex] ;
      Integer coilIndex = coilRegionId_.Find(actRegionId);
      Double conductivity;
      materials_[actRegionId]->GetScalar(conductivity,MAG_CONDUCTIVITY,Global::REAL);
      ElemList actSDList(ptgrid_);
      actSDList.SetRegion(actRegionId);

      EntityIterator it = actSDList.GetIterator();
      UInt actEl = 0;
      // iterate over all elements of regions
      for ( it.Begin(); !it.IsEnd(); it++, ++actEl ) {

        BaseFE * ptElem = it.GetElem()->ptElem;
        ptElem->SetAnsatzFct( results_[0]->fctType );
        UInt numFncs = ptElem->GetNumFncs(results_[0]->fctType);
        const UInt nrIntPts = ptElem->GetNumIntPoints();
        const Vector<Double> & intWeights = ptElem->GetIntWeights();  

        Vector<Double> shpFncAtIp;         
        ptgrid_->GetElemNodesCoord( ptCoord, it.GetElem()->connect, true );

        // iterate over all integration points
        fAtIp.Resize( dim_ * numFncs );
        elemForce.Resize( dim_ * numFncs );
        elemForce.Init();
        for( UInt ip = 1; ip <= nrIntPts; ++ip ) {
          ptElem->GetShFncAtIp(shpFncAtIp, ip, it.GetElem() );

          // CHECK: If this region is a current coil, we simply take the
          // prescribed current density value
          if( coilIndex != -1 ) {
            MathParser * mParser =  domain->GetMathParser();
            std::string factor = coilDef_[coilIndex]->value_ + "/" 
            + lexical_cast<std::string>(coilDef_[coilIndex]->windingCrossSection_ );
            mParser->SetExpr( mHandle_, factor );
            Double currDens = mParser->Eval(mHandle_);
            if( is3d_ ) {
              // take flow direction into account
            } else {
              jAtIp.Resize(1);
              jAtIp[0] = currDens;
            }
          } else {
            // calculate jEddy at integration point
            CalcEddyCurrentAtIP( it, ip, jAtIp );
          }

          // calculate flux density at ip
          CalcFluxDensityAtIP( it, ip, fluxAtIp );

          // calculate cross product 
          fAtIp.Init();
          if( is3d_ ) {
            Vector<TYPE> tempCross;
            jAtIp.CrossProduct( fluxAtIp, tempCross );
            for (UInt iFnc=0; iFnc<numFncs; ++iFnc) {
              fAtIp[iFnc*3+0] = -tempCross[0] * shpFncAtIp[iFnc];
              fAtIp[iFnc*3+1] = -tempCross[1] * shpFncAtIp[iFnc];
              fAtIp[iFnc*3+2] = -tempCross[2] * shpFncAtIp[iFnc];
            }
          } else if (isaxi_ ) {

            // calculate pseudo cross product for axi-symmetric case
            for (UInt iFnc=0; iFnc<numFncs; ++iFnc) {
              fAtIp[iFnc*2] -= jAtIp[0] * fluxAtIp[1] 
                                                   * shpFncAtIp[iFnc];
              fAtIp[iFnc*2+1] += jAtIp[0] * fluxAtIp[0]   
                                                     * shpFncAtIp[iFnc];
            } 
          } else{
            // calculate pseudo cross product
            for (UInt iFnc=0; iFnc<numFncs; ++iFnc) {
              fAtIp[iFnc*2] += jAtIp[0] * fluxAtIp[1] 
                                                   * shpFncAtIp[iFnc];
              fAtIp[iFnc*2+1] -= jAtIp[0] * fluxAtIp[0]   
                                                     * shpFncAtIp[iFnc];
            } 
          }
          Double jacDet = ptElem->CalcJacobianDetAtIp(ip, ptCoord, 
              it.GetElem() );
          if( isaxi_ ) {
            Vector<Double> CoordAtIP;
            CoordAtIP = ptCoord * shpFncAtIp;
            jacDet *=  2 * PI * CoordAtIP[0];
          }
          elemForce -= fAtIp * (jacDet * intWeights[ip-1] );

        }
        StdVector<UInt> const & connecth = it.GetElem()->connect;

        // Add the element force to the according coupling node
        for (UInt ielemnode=0, n=connecth.GetSize(); ielemnode<n; ++ielemnode) {
          std::map<UInt, UInt>::const_iterator itPos
              = nodeNumPos.find(connecth[ielemnode]);
          if (itPos != itPosEnd) {
            UInt pos = itPos->second;
            for( UInt idim=0; idim<dim_; ++idim) {
              force[pos*dim_+idim] += elemForce[ielemnode*dim_+idim];
            }
          }
        }
      }


    }
  }


  bool MagPDE::HasOutput( SolutionType output ) {
    if (output == MAG_FORCE_LORENTZ) {
      return true;
    }
    else if (output == MAG_FORCE_VWP) {
      return true;
    }
    return false;
  }

} // end of namespace

