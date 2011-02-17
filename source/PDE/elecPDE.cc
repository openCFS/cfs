// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <set>

#include "elecPDE.hh"

#include "Forms/linGradBDBInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Forms/nLinElecHystInt.hh"
#include "Forms/nLinElecMaterial.hh"
#include "Forms/elecchargeop.hh"
#include "Forms/laplaceInt.hh"
#include "Forms/FlatShellElecInt.hh"
#include "Forms/nonConformingInt.hh"
#include "Forms/singleEntryInt.hh"
#include "Forms/massInt.hh"
#include "Forms/gradfieldop.hh"
#include "Forms/piezoPolarizationMatrixRHSInt.hh"
#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Utils/StdVector.hh"
#include "Driver/solveStepElec.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/hysteresis.hh"
#include "Utils/preisach.hh"
#include "Optimization/Design/DesignSpace.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {

  DECLARE_LOG(elecpde)
  DEFINE_LOG(elecpde, "pde.electrostatic")

  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE( Grid* aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "electrostatic";
    pdematerialclass_ = ELECTROSTATIC;
    maxTimeDerivOrder_ = 0;
 
    nonLin_    = false;
    nonLinMaterial_ = false;
    isAlwaysStatic_ = true;
    isPiezoCoupled_ = false;
    isThermoCoupled_= false;

    needSolPrev_ = true;
 
    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);
  }
  

  void ElecPDE::InitNonLin() {

    // Check, if "nonLinList" is present
    PtrParamNode nonLinListNode = myParam_->Get("nonLinList", ParamNode::PASS );
    if( nonLinListNode ) { 

      // Get nonlinear types
      ParamNodeList nonLinNodes = nonLinListNode->GetChildren();
      for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

        std::string actTypeString = nonLinNodes[i]->GetName();
        std::string actId = nonLinNodes[i]->Get("id")->As<std::string>();

        NonLinType actType;
        String2Enum( actTypeString, actType );
        nonLinIdType_[actId] = actType;
      }
    }
    
    // Run over all region and set entry in "regionNonLinId"
    ParamNodeList regionNodes = 
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
      
      regionNodes[i]->GetValue( "name", actRegionName );
      regionNodes[i]->GetValue( "nonLinId", actNonLinId );
      
      if( actNonLinId == "" )
        continue;

      actRegionId = ptgrid_->GetRegion().Parse(actRegionName);

      // Check nonLinId was already registerd
      if( nonLinIdType_.find( actNonLinId) == nonLinIdType_.end() ) {
        EXCEPTION( "NonLinearity with id '" << actNonLinId 
                   << "' was not defined in 'nonLinList'" );
      }
      
      regionNonLinId_[actRegionId] = actNonLinId;

      // get related type of nonlinearity
      NonLinType actType = nonLinIdType_[actNonLinId];
      regionNonLinType_[actRegionId] = actType;

      // check type
      if( actType == HYSTERESIS ) {
        isHysteresis_ = true;
      }

      if( actType == MATERIAL ) {
        nonLin_ = true;
        nonLinMaterial_ = true;
      }

      // Log to info file
      std::string nonLinString;
      Enum2String( nonLinIdType_[actNonLinId], nonLinString );
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(), 
                    nonLinString.c_str() );
      
    }
    Info->PrintF( pdename_, "\n" );

  }

  void ElecPDE::DefineIntegrators() {

    
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // flag for updatedLagrange formulation
    bool upLagrangeForm = true;
 
    //transform the type
    SubTensorType tensorType;

    if ( dim_ == 3 ) {
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

    // if the pde is piezo-coupled, the electrostatic entries
    // have to multiplied with -1
    std::string factor = "1.0";
    if ( isPiezoCoupled_ == true || isThermoCoupled_== true )
      factor = "-1.0";  

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region node
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      // isPiezoHyst (isMicroPiezo) = true means, that the bilinear-forms 
      // will be defined in PiezoCupling.cc!!
      bool isPiezoHyst  = false;
      bool isMicroPiezo = false;
      if ( isPiezoCoupled_ == true ) {
        isPiezoHyst = IsRegionPiezoHyst( regionName );
        isMicroPiezo = IsRegionMicroPiezo( regionName );
      }

      if ( !isPiezoHyst && !isMicroPiezo ) {
        // check for nonlinearity
        if ( regionNonLinType_[actRegion] == HYSTERESIS) {
          StdVector<Elem*> elemssd;
          ptgrid_->GetElems(elemssd, actRegion);
          UInt numElSD =  elemssd.GetSize();
          
          
          //allocate for hystersis modeling
          actSDMat->InitHyst(numElSD, actSDList);
          
          nlinElecHystInt* nlForm;
          nlForm = new nlinElecHystInt( actSDMat, tensorType,
                                        upLagrangeForm );

          nlForm->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
          nlForm->Set4Hyst(ptgrid_, this, eqnMap_,  results_[0]);
          nlForm->SetFactor(factor);
          
          BiLinFormContext * stiffIntDescr = 
            new BiLinFormContext(nlForm, STIFFNESS );
          
          stiffIntDescr->SetPtPdes(this, this);
          stiffIntDescr->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
          
          assemble_->AddBiLinearForm( stiffIntDescr );
        }
        
        else if (  regionNonLinType_[actRegion] == MATERIAL ) {
          
          std::string nlfnc = materials_[actRegion]->GetNonlinFileName();
          materials_[actRegion]->GetScalar(nlfnc,NONLIN_DATA_NAME);           
          
          ApproxData *nlinFnc = new SmoothSpline(nlfnc);
          //  oder        ApproxData *nlinFnc = new LinInterpolate(nlfnc);
          nlinFnc->CalcBestParameter();
          nlinFnc->CalcApproximation();
          
          //         if (dim_ == 3)
          //           {   
          nLinElec3dInt_Material* nLinMaterial;
          nLinMaterial = new nLinElec3dInt_Material(nlinFnc, actSDMat, FULL);    
          //          }
          
          nLinMaterial->SetSolution(dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
          nLinMaterial->Set4NonLinMaterial(ptgrid_, this, eqnMap_,  results_[0]);
          
          BiLinFormContext * stiffNLMaterialDescr = 
            new BiLinFormContext(nLinMaterial, STIFFNESS );
          
          stiffNLMaterialDescr->SetPtPdes(this, this);
          stiffNLMaterialDescr->SetResults( results_[0], results_[0],
                                            actSDList, actSDList );
          assemble_->AddBiLinearForm(stiffNLMaterialDescr);
        }
        
        else {

          // --- standard real-valued stiffness integrator ---
          linGradBDBInt *  linElecForm = 
              new linGradBDBInt( actSDMat, ELEC_PERMITTIVITY, 
                                 tensorType, upLagrangeForm );
          linElecForm->SetFactor( factor );
          
          BiLinFormContext * stiffIntDescr = 
            new BiLinFormContext(linElecForm, STIFFNESS );
          
          stiffIntDescr->SetPtPdes(this, this);
          stiffIntDescr->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
          
          assemble_->AddBiLinearForm( stiffIntDescr );
          
          
          // --- check for complex valued material parameter ---
          // check for piezo-coupling and 
          if( complexMatData_[actRegion] == true ) {
            Global::ComplexPart matType = Global::IMAG; 

            linGradBDBInt * linElecFormC  = new 
                linGradBDBInt( materials_[actRegion], 
                               ELEC_PERMITTIVITY,
                               tensorType,
                            upLagrangeForm);
            linElecFormC->SetFactor( factor );
            linElecFormC->SetMatDataType(matType);

            BiLinFormContext * stiffIntDescrC = 
                new BiLinFormContext(linElecFormC, STIFFNESS);

            stiffIntDescrC->SetPtPdes(this, this);
            stiffIntDescrC->SetEntryType(matType);
            stiffIntDescrC->SetResults( results_[0], results_[0],
                                        actSDList, actSDList );
            assemble_->AddBiLinearForm( stiffIntDescrC );
          }
        }
      }  
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
      
    }

    // Define integrators for composite materials
    // (only for subType "flatShell")
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
         compIt++ ) {
      
      // Get current subdomain and composite material
      RegionIdType actRegion = compIt->first;
      Composite * composite = &compIt->second;
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      
      FlatShellElecInt * compElecInt = new FlatShellElecInt( composite, false );
      BiLinFormContext * stiffContext = 
        new BiLinFormContext( compElecInt, STIFFNESS);
      stiffContext->SetPtPdes( this, this );
      stiffContext->SetResults( results_[0], results_[0],
                                actSDList, actSDList );
      assemble_->AddBiLinearForm( stiffContext );
      eqnMap_->AddResult( *results_[0], actSDList );
     
      // --- check for complex valued material parameter ---
      if( complexMatData_[actRegion] == true ) {
        Global::ComplexPart matType = Global::IMAG;
        FlatShellElecInt * compElecIntC = new FlatShellElecInt( composite, false );
        compElecIntC->SetMatDataType(matType);
        BiLinFormContext * stiffContextC = 
          new BiLinFormContext( compElecIntC, STIFFNESS);
        stiffContextC->SetPtPdes( this, this );
        stiffContextC->SetEntryType(matType);
        stiffContextC->SetResults( results_[0], results_[0],
                                  actSDList, actSDList );
        assemble_->AddBiLinearForm( stiffContextC );
      }
      
    }

    // **********************************************************************
    //   inhom. Neumann boundary condition
    // **********************************************************************
    for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) 
    {
      InhomNeumannBc const & actBc = *inBcs_[iBc];

      
      LinearSurfForm *neumannBC = new LinNeumannInt(actBc.value, actBc.phase, 
                                                     NO_MATERIAL, isaxi_);
      neumannBC->SetVoluInfo(materials_);
      LinearFormContext* neumannContext = new LinearFormContext(neumannBC);
      neumannContext->SetPtPde(this);
      neumannContext->SetResult(actBc.result, actBc.entities);
      assemble_->AddLinearForm(neumannContext);
      
      // Give result to equation numbering class
      eqnMap_->AddResult(*actBc.result, actBc.entities);
    }
    
    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    
    // Get index of LAGRANGE_MULT result, just in case... who knows...
    UInt lmResultIdx = 0;
    for(UInt i=0, n=results_.GetSize(); i<n; i++) {
      if(results_[i]->resultType == LAGRANGE_MULT) {
        lmResultIdx = i;
        break;
      }
    }
    LOG_DBG2(elecpde) << "NonMatching: Index of LAGRANGE_MULT result: "
                     << lmResultIdx;
    
    for( UInt i = 0; i < ncIFaces_.GetSize(); i++ ) {
      
      // get regionId of Lagrangian surface
      StdVector<std::string> keyVec, attrVec, valVec;
      std::string slaveSide;
      std::string ncIfaceName = ptgrid_->GetRegion().ToString(ncIFaces_[i]);

      PtrParamNode ncIfaceListNode;
      ncIfaceListNode = param->Get("domain")->Get("ncInterfaceList");

      slaveSide = ncIfaceListNode->
          GetByVal("ncInterface", "name",
                   ncIfaceName)->Get("slaveSide")->As<std::string>();

      // Part 1: Define integrator M(Psi, Lambda) on
      //         non-conforming interface
      LOG_DBG2(elecpde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
      actNcList->SetRegion( ncIFaces_[i] );
      
      NonConformingInt * ncInt = 
        new NonConformingInt( 1, isaxi_ );

      NcBiLinFormContext * stiffIntDescr = 
     	  new NcBiLinFormContext( ncInt , STIFFNESS );

      // Force assembling of M(Psi, Lambda)^T
      stiffIntDescr->SetCounterPart( true );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[lmResultIdx],
                                 actNcList, actNcList );
      
      assemble_->AddBiLinearForm( stiffIntDescr );


      // Part 2: Define integrator D(Psi, Lambda) on
      //         Lagrangian surface
      LOG_DBG2(elecpde) << "NonMatching: Defining mass integrator"
                        << " for D on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );      
      actSDList->SetRegion( ptgrid_->GetRegion().Parse(slaveSide));

      // D(Psi, Lambda) has the form of a standard mass
      // integrator with factor 1.0
      MassInt * dMatInt = new MassInt( 1.0, 1, isaxi_ );
      BiLinFormContext * dMatContext = 
        new BiLinFormContext( dMatInt, STIFFNESS );

      // Force assembling of D(Psi, Lambda)^T
      dMatContext->SetCounterPart( true );
      dMatContext->SetPtPdes( this, this );
      dMatContext->SetResults( results_[0], results_[lmResultIdx],
                               actSDList, actSDList );
      
      assemble_->AddBiLinearForm( dMatContext );

      // Give result LAGRANGE_MULT to equation numbering class
      eqnMap_->AddResult( *results_[lmResultIdx], actSDList );
    }


    // define integrators for electric impedances
    DefineImpedanceIntegrators();
  }
  
  void ElecPDE::DefineImpedanceIntegrators() {
    
    // The definition of the impedances is taken from:

    // Jian S. Wang and Dale F. Ostergaard:
    // A Finite Element-Electric Circuit Coupled Simulation Method for
    // Piezoelectric Transducer
    // IEEE Ultrasoncics Symposium, 1999

    // =======================================================================
    // Integrators for electric impedances
    // =======================================================================
    
    for( UInt i = 0; i < impedances_.GetSize(); i++ ) {
      
      std::string real, imag, temp;
      real = "";
      imag = "";
      
      // *** Resistance ***
      if ( impedances_[i].resistance > EPS ) {
        imag = "( 1.0 / ((2*pi*f) * ";
        imag += lexical_cast<std::string>(impedances_[i].resistance)+ " ))";
      }
      
      // *** Inductance ***
      if ( impedances_[i].inductance > EPS ) {
        real = "(1.0 / ((2*pi*f) * (2*pi*f) * ";
        real += lexical_cast<std::string>(impedances_[i].inductance) + " ))";
      }
      
      // *** Capcitance *** 
      if ( impedances_[i].capacitance > EPS ) {
        real += "-" + lexical_cast<std::string>( impedances_[i].capacitance );
      } 
      
      // perform consistency check
      if( real == "" ) real = "0.0";
      if( imag == "" ) imag = "0.0";

      // generate integrators
      // a) diagonal entries
      SingleEntryInt * stiffInt1 = 
        new SingleEntryInt( real, imag,  1, 1 );
      BiLinFormContext * stiffIntContext1 = 
        new BiLinFormContext( stiffInt1, STIFFNESS );
      stiffIntContext1->SetPtPdes(this, this);
      stiffIntContext1->SetResults( results_[0], results_[0],
                                   impedances_[i].node1,
                                   impedances_[i].node1 );
      assemble_->AddBiLinearForm( stiffIntContext1 );

      SingleEntryInt * stiffInt2 = 
        new SingleEntryInt( real, imag,  1, 1 );
      BiLinFormContext * stiffIntContext2 = 
        new BiLinFormContext( stiffInt2, STIFFNESS );
      stiffIntContext2->SetPtPdes(this, this);
      stiffIntContext2->SetResults( results_[0], results_[0],
                                    impedances_[i].node2,
                                    impedances_[i].node2 );
      assemble_->AddBiLinearForm( stiffIntContext2 );
      
      // b) off-diagonal entries
      SingleEntryInt * stiffInt3 = 
        new SingleEntryInt( "-"+real, "-"+imag,  1, 1 );
      BiLinFormContext * stiffIntContext3 = 
        new BiLinFormContext( stiffInt3, STIFFNESS );
      stiffIntContext3->SetPtPdes(this, this);
      stiffIntContext3->SetResults( results_[0], results_[0],
                                    impedances_[i].node1,
                                    impedances_[i].node2 );
      stiffIntContext3->SetCounterPart( true );
      assemble_->AddBiLinearForm( stiffIntContext3 );
    }
  }

  void ElecPDE::DefineSolveStep() {
    solveStep_ = new SolveStepElec(*this);
  }

  void ElecPDE::ReadSpecialBCs( ) 
  {
     ReadImpedances();
  }
  
  

  void ElecPDE::ReadImpedances( ) 
  {
  
    // Get values from parameter file
    ParamNodeList impedNodes = 
      myParam_->Get("bcsAndLoads")->GetList("impedance");
   
    // make sure we are in harmonic mode
    if( impedNodes.GetSize() > 0 && analysistype_ != HARMONIC) {
      EXCEPTION( "Impedances are only available in harmonic simulation "
                 << "at the moment!" );
    }

    std::string node1, node2;
    Double resistance, capacitance, inductance;
    StdVector<UInt> nodeVec, oneNodeVec(1);
    
    // resize impedances and fill in data
    PtrParamNode list = infoNode_->Get(ParamNode::HEADER)->Get("impedances");
    std::ostringstream out;
    impedances_.Resize( impedNodes.GetSize() );
    for( UInt i = 0; i < impedNodes.GetSize(); i++ )
    {
      PtrParamNode node = impedNodes[i];
      
      // get data from node
      node1       = node->Get("node1")->As<std::string>();
      node2       = node->Get("node2")->As<std::string>();
      resistance  = node->Get("resistance")->As<Double>();
      inductance  = node->Get("inductance")->As<Double>();
      capacitance = node->Get("capacitance")->As<Double>();

      Impedance& imp = impedances_[i];
      
      // fill data for impedance node
      imp.node1 = shared_ptr<NodeList>( new NodeList(ptgrid_) );
      ptgrid_->GetNodesByName( nodeVec, node1);
      oneNodeVec[0] = nodeVec[0];
      imp.node1->SetNodes( oneNodeVec );
      imp.node2 = shared_ptr<NodeList>( new NodeList(ptgrid_) );
      ptgrid_->GetNodesByName( nodeVec, node2);
      oneNodeVec[0] = nodeVec[0];
      imp.node2->SetNodes( oneNodeVec );
      imp.resistance = resistance;
      imp.capacitance = capacitance;
      imp.inductance = inductance;

      PtrParamNode in = list->Get("impedance", ParamNode::APPEND);
      in->Get("node1")->SetValue(node1);
      in->Get("node2")->SetValue(node2);
      in->Get("resistance")->SetValue(imp.resistance);
      in->Get("inductance")->SetValue(imp.inductance);
      in->Get("capacitance")->SetValue(imp.capacitance);
    }
  }
  
  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  // ***************
  //   CalcResults
  // ***************
   void ElecPDE::CalcResults( shared_ptr<BaseResult> res ) {
     
     switch (res->GetResultInfo()->resultType ) {
     case ELEC_POTENTIAL:
       if( isComplex_ ) {
         ExtractResult<Complex>( res, sol_ );
       } else {
         ExtractResult<Double>( res, sol_ );
       }
       break;

     case ELEC_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( res, results_[0] );
      } else {
        ExtractRhsResult<Double>( res, results_[0] );
      }
      break;
       
     case ELEC_FIELD_INTENSITY:
       if( isComplex_ ) {
         CalcElectricField<Complex>( res );
       } else {
         CalcElectricField<Double>( res );
       }
       break;
       
     case ELEC_FLUX_DENSITY:
       if( isComplex_ ) {
         CalcElectricFluxDensity<Complex>( res );
       } else {
         CalcElectricFluxDensity<Double>( res );
       }
       break;

     case ELEC_ENERGY:
       if( isComplex_ ) {
         CalcEnergy<Complex>( res );
       } else {
         CalcEnergy<Double>( res );
       }
       break;

     case ELEC_CHARGE:
       if( isComplex_ ) {
         CalcCharges<Complex>( res );
       } else {
         CalcCharges<Double>( res );
       }
       break;

     case ELEC_POLARIZATION:
       CalcPolarizationField( res );
       break;
       
     case ELEC_PSEUDO_POLARIZATION:
       if(domain->GetErsatzMaterial(false) == NULL) // no excpetion
         res->Init();
       else     
         domain->GetErsatzMaterial()->ExtractResults(res, isComplex_);
       break;

     case GRAD_ELEC_POTENTIAL:
       if(isComplex_)
         CalcGradSolution<Complex>(res);
       else
         CalcGradSolution<Double>(res);
       break;

     default:
       WARN( "Result type not computable by electric PDE" );
     }
     
   }


  template <class TYPE>
  void ElecPDE::CalcElectricField( shared_ptr<BaseResult> sol )
  {
    Vector<Double> lCoord;
    Vector<TYPE> tempE;
    NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
    GradientFieldOp<TYPE> * FieldOp = 
      new GradientFieldOp<TYPE>(ptgrid_, this, 
                                eqnMap_, solhelp, 
                                results_[0]->fctType,isaxi_);
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*sol);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
    
    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      it.GetElem()->ptElem->GetCoordMidPoint(lCoord);
      FieldOp->CalcElemGradField( tempE, it, lCoord, 1); 
      
      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
        actVal[it.GetPos()*dim_ + iDim] = tempE[iDim];
      }
    }
    delete FieldOp;
    
  }
  
  template <class TYPE>
  void ElecPDE::CalcElectricFluxDensity( shared_ptr<BaseResult> sol )
  {
     
    // first calculate electric field intensity
    CalcElectricField<TYPE>( sol );
    
    // then, run over all elements and multiply by element permittivity
    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*sol);
      EntityIterator it = actSol.GetEntityList()->GetIterator();

      Vector<TYPE> & actVal = actSol.GetVector();

      Matrix<Double> permit;
      BaseMaterial * ptMat = NULL;
      Vector<TYPE> tempE (dim_), tempD (dim_);
      SubTensorType tensorType;
      if( dim_ == 2 ) {
        tensorType = PLANE;
      } else {
        tensorType = FULL;
      }
      
      // loop over elements
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        ptMat = materials_[it.GetElem()->regionId];
        ptMat->GetTensor( permit, ELEC_PERMITTIVITY, Global::REAL, tensorType );
        
        // build temporary element vector
        tempE.Init();
        for(UInt iDim = 0; iDim < dim_; iDim++ ) {
          tempE[iDim] = actVal[it.GetPos()*dim_ + iDim];
        }

        tempD = permit * tempE;

        // loop over dofs
        for(UInt iDim = 0; iDim < dim_; iDim++ ) {
          actVal[it.GetPos()*dim_ + iDim] = tempD[iDim];
        }
      }

  }
  
  


  void ElecPDE::CalcPolarizationField( shared_ptr<BaseResult> res ) {

    //we assume, that the actual solution is stored in sol_!
    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);

    GradientFieldOp<Double> * FieldOp = 
      new GradientFieldOp<Double>(ptgrid_,this, eqnMap_,
 				  *solhelp, results_[0]->fctType, isaxi_);
  
    Vector<Double> LCoord, Efield, vecE, vecP;
    Double Ecomp, Pval;
    UInt ielGlob;
    vecP.Resize(2);

    // Fetch current result object and resize it correctly
    Result<Double> &  actSol = 
      dynamic_cast<Result<Double>&>(*res);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    
    Vector<Double> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
    actVal.Init();

    // determine material from first element
    it.Begin();
    BaseMaterial * actMat = materials_[it.GetElem()->regionId];
    
    Hysteresis* hyst = actMat->getHysteresis();
    if ( hyst!= NULL ) {
      // get direction of polarization
      std::string str;
      actMat->GetScalar(str, P_DIRECTION);
      Directions dirP;
      String2Enum(str,dirP);

      // iterate over all element in list
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        
        //compute the electric field intensity
        it.GetElem()->ptElem->GetCoordMidPoint(LCoord);
        FieldOp->CalcElemGradField( Efield, it, LCoord, 1);
        
        vecE = Efield / Efield.NormL2();
        //get correct component of electric field for scalar Preisach model
        //materials_[ subdoms_[isd] ]->
        Ecomp = Efield[dirP]; //.NormL2(); //[comp]; 
        ielGlob = it.GetElem()->elemNum ;
        Pval = actMat->ComputeScalarHystVal( ielGlob, Ecomp );
        vecP.Init();
//          for ( UInt i=0; i<vecE.GetSize(); i++)
//            vecE[i] = abs( vecE[i] );

//         vecP = vecE * Pval;
        vecP[dirP] = Pval;
        for( UInt iDof = 0; iDof < dim_; iDof++ ) {
          actVal[it.GetPos()*dim_ + iDof] = vecP[iDof];
        }
      }  
    } else {
      WARN( "No hysteresis defined for list '"
            << actSol.GetEntityList()->GetName() );
    }
  }



  template <class TYPE>
  void ElecPDE::CalcCharges( shared_ptr<BaseResult> res ) {


    NodeStoreSol<TYPE> * solhelp = dynamic_cast<NodeStoreSol<TYPE>*>(sol_);
    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<TYPE> elemDField;
    Elem * ptVolElem = 0;
    BaseFE * ptSurfElemFE = 0, * ptVolElemFE = 0;
    Double permittivity = 0.0;
    TYPE elemNormalD = 0.0;
    TYPE charge = 0.0;
    Double normSign = 0;
    UInt pdeElemNum = 0;

    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);
    
    // Create operator for electric flux density and charge calculation
    GradientFieldOp<TYPE> * dFieldOp = 
      new GradientFieldOp<TYPE>(ptgrid_, this, eqnMap_, *solhelp,
                                results_[0]->fctType, isaxi_);
    ElecChargeOp<TYPE> * chargeOp = 
      new ElecChargeOp<TYPE>(ptgrid_, this, eqnMap_, results_[0], isaxi_ );
    
    
    // do cast of restult and resize solution vector
    Result<TYPE> &  actRes = 
      dynamic_cast<Result<TYPE>&>(*res);
    EntityIterator it = actRes.GetEntityList()->GetIterator();
    
    Vector<TYPE> & actVal = actRes.GetVector();
    actVal.Resize( actRes.GetEntityList()->GetSize() );
    
    // loop over all surface elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      // Determine, which volume element is the right neighbour for the 
      // calculation:
      // We assume the normal to point out of the surface element into the
      // direction of the neighboring electric volume element.
      if ( surfNeighborRegions_[res] ==
           it.GetSurfElem()->ptVolElem1->regionId) {
        ptVolElem = it.GetSurfElem()->ptVolElem1;
        normSign = -1.0;
      } else {
        ptVolElem = it.GetSurfElem()->ptVolElem2;
        normSign = 1.0;
      }
      
      // Check, that a volume element at all is present
      if( ptVolElem == NULL ) {
        EXCEPTION( "The surface element #" << it.GetSurfElem()->elemNum
                   << " does not belong to a region defined in the "
                   << "electrostatic Pde" );
      }


      normSign *= (Double) it.GetSurfElem()->normalSign;
      
      ptSurfElemFE = it.GetSurfElem()->ptElem; 
      ptVolElemFE = ptVolElem->ptElem;
      
      const StdVector<UInt> & surfConnect = it.GetSurfElem()->connect;
      const StdVector<UInt> & volConnect = ptVolElem->connect;
      
      // calculate volume integration coordinates from
      // surfe integration coordinat for evalauting the 
      // electric flux density on the surface of the volume
      // element
      ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
      ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                             lCoordSurf, lCoordVol);

      // get local element number for volume element (for hysteresis value)
      pdeElemNum = eqnMap_->Mesh2PdeElem(ptVolElem->elemNum);
      
      // Get the right material parameter for actual volume element
      BaseMaterial * myMat = materials_[ptVolElem->regionId];
      myMat->GetScalar(permittivity,ELEC_PERMITTIVITY,Global::REAL);

      // check for hysteresis
      Hysteresis* hyst = NULL;
      hyst = myMat->getHysteresis();      

      // Calc electric flux density
      ElemList tempList(ptgrid_);
      tempList.SetElement( ptVolElem );
      EntityIterator tempIt = tempList.GetIterator();

      // D = eps_0 E * P
      if ( hyst ) {
        //compute electric field
        dFieldOp->CalcElemGradField(elemDField, tempIt, 
                                    lCoordVol,1.0);
        // get direction of polarization
        std::string str;
        myMat->GetScalar(str, P_DIRECTION);
        Directions dirP;
        String2Enum(str,dirP);
        
        Vector<Double> vecP(dim_);
        vecP.Init();
        vecP[dirP] = myMat->GetScalarHystVal( pdeElemNum );
        
        // D = eps_0 E * P
        elemDField *= 8.854e-12;
        elemDField =  elemDField + vecP;
      } else {
        dFieldOp->CalcElemGradField(elemDField, tempIt, 
                                    lCoordVol,permittivity);
      }
      
      // Calc global normal
      ptgrid_->CalcSurfNormal(normal, *(it.GetSurfElem()));
      
      normal *= normSign;
      
      // Since the routine CalcLineNormal always computes a normal
      // which points in the OPPOSITE direction of the volume element,
      // we have to multiply the normal with -1 to get the correct sign for
      // the charges
      elemNormalD = elemDField * normal;
      
      chargeOp->CalcElemCharge(charge, it, 
                               lCoordSurf, elemNormalD);
      actVal[it.GetPos()] = charge;
    }
    delete chargeOp;
    delete dFieldOp;
  }

  template <class TYPE>
  void ElecPDE::CalcEnergy( shared_ptr<BaseResult> vals )
  {

    Matrix<Double> elemmat;  
    SubTensorType tensorType;

    if ( dim_ == 3 ) {
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
   
    Vector<TYPE> help;
    std::string factor = "1.0";
    TYPE energy;

    Result<TYPE> &  actSol = 
      dynamic_cast<Result<TYPE>&>(*vals);      
    EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
    
    // resize vector
    Vector<TYPE> & actVal = actSol.GetVector();
    actVal.Resize( actSol.GetEntityList()->GetSize() );
    
    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();
      
      linGradBDBInt * bilinear_stiff = 
        new linGradBDBInt(materials_[regionIt.GetRegion()],ELEC_PERMITTIVITY,
                          tensorType);

      // Loop over elements
      energy = 0;
      Vector<TYPE> elpot;
      for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        
        bilinear_stiff->SetFactor(factor);
        bilinear_stiff->CalcElementMatrix( elemmat, elemIt, elemIt );
        
        sol_->GetElemSolution(elpot, elemIt );
        help =  elemmat * elpot;
        energy += 0.5 * (help * elpot);
        
        LOG_DBG3(elecpde) << "CalcEnergy: elem=" << elemIt.GetElem()->elemNum << " mat=" << elemmat.ToString();
        LOG_DBG3(elecpde) << "CalcEnergy: elemIt=" << elemIt.GetElem()->elemNum << " elpot=" << elpot.ToString();
        LOG_DBG3(elecpde) << "CalcEnergy: elemIt=" << elemIt.GetElem()->elemNum << " sum energy -> " << energy;
      }  
      actVal[regionIt.GetPos()] = energy;
      delete bilinear_stiff;
    }
  }


  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void ElecPDE::InitCoupling(PDECoupling * Coupling)
  {
  
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    
    StdVector<std::string> * nRegions;
    StdVector<RegionIdType> nRegionIds;
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();
    
    nonLin_ = false;

    // Initialization of coupling helper arrays
    StdVector<UInt> * couplingnodes = NULL;

    for (UInt actCoupling=0; actCoupling<numCouplings; actCoupling++) {
      // Initialize arrays for coupling surface elements
      if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP
          || ptCoupling_->GetOutputQuantity(actCoupling) == 
          ELEC_INTERFACE_FORCE) {
      
        ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
        if (couplingnodes == NULL)
          std::cerr << "Couplingnodes = 0!!!!" << std::endl;
      
        // if quantity is elecFocrceVWP, get volume neighbours lying next to
        // coupling nodes, because these volume elements have to be 
        // moved 'virtually'
        if (ptCoupling_->GetOutputQuantity(actCoupling) == ELEC_FORCE_VWP) {
          NodeStoreSol<Double> * solhelp = 
            dynamic_cast<NodeStoreSol<Double> *>(sol_);
          ForceOp_ = new  ElecForceOp(ptgrid_, this,  eqnMap_, 
                                      *solhelp, dim_, materials_, 
                                      results_[0], isaxi_, true );

          ptCoupling_->GetOutputNeighbourRegion(actCoupling, nRegions);
          ptgrid_->GetRegion().Parse(*nRegions, nRegionIds);
          ForceOp_->Setup(nRegionIds, *couplingnodes);
        }
      
        else if ( ptCoupling_->GetOutputQuantity(actCoupling)
                  == ELEC_INTERFACE_FORCE) {
          EXCEPTION( "Currently ELEC_INTERFACE_FORCE not supported" );
        }

        else {
          EXCEPTION( "Coupling " << SolutionTypeEnum.ToString(ptCoupling_->GetOutputQuantity(actCoupling)) <<  " not known! ");
        }
      
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);
      
      
      } // end for (actNode)
    } 
  }
  


  void ElecPDE::CalcOutputCoupling()
  {

    SolutionType quantity;
    StdVector<UInt> * couplingNodes     = NULL;
    SingleVector * values = 0;
    UInt forcesCount = 0;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for (UInt actCoupling=0; 
         actCoupling<ptCoupling_->GetNumOutputCouplings(); 
         actCoupling++)
      {
        quantity = ptCoupling_->GetOutputQuantity(actCoupling);
        ptCoupling_->GetOutputValues(actCoupling, values);

        Vector<Double> * temp = dynamic_cast<Vector<Double> *>(values);
      
        switch(ptCoupling_->GetOutputType(actCoupling))
          {
          
          case NODE:        
            ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
          
            if (quantity == ELEC_POTENTIAL)
              sol_->NodeSolutionToCoupling(*values, *couplingNodes);
            
            if (quantity == ELEC_FORCE_VWP)
              {
                Vector<Double> totalForce;
                ForceOp_->CalcNodeForce(*temp, totalForce);

                // write information in .info-file
                Info->PrintF(pdename_, "Sum of electrostatic force (VWM):\n");
                Info->PrintVec(totalForce);
                forcesCount++;
              }
                  
            break;
          
          case ELEM:
            EXCEPTION("No Element coupling output");
          }
      }
  }


  bool ElecPDE::HasOutput(SolutionType output)
  {
  
    switch (output)
      {
      case ELEC_FORCE_VWP:
        return true;
        break;
      case ELEC_POTENTIAL:
        return true;
        break;
      case ELEC_FIELD_INTENSITY:
        return true;
        break;
      case ELEC_INTERFACE_FORCE:
        return true;
        break;
      default:
        return false;
        break;
      }
    return false;
  }


  void ElecPDE::SetPiezoCoupling()
  {
  
    isPiezoCoupled_ = true;

  }
  
  void ElecPDE::SetThermoCoupling()
  {
  
    isThermoCoupled_ = true;

  }
  
  void ElecPDE::DefinePolarizationMatrixIntegrators(const Vector<Double> &vals,
      StdVector<LinearFormContext*> *linForms, const int num)
  {
    LinearForm * polRHSElec = new PiezoPolarizationMatrixElecRHSInt(vals, num);
    
    StdVector<std::string> region_names;
    domain->GetGrid()->GetRegionNames(region_names);
    LOG_DBG(elecpde) << "region names = " << region_names.ToString();
    // FIXME: hardcoded region name!!
    shared_ptr<EntityList> entlist = domain->GetGrid()->GetEntityList(EntityList::SURF_ELEM_LIST, 
        region_names[1], EntityList::NO_TYPE);
    
    LinearFormContext *linRhs = new LinearFormContext(polRHSElec);
    linRhs->SetPtPde(this);
    linRhs->SetResult(results_[0], entlist);
    
    if(linForms != NULL)
      linForms->Push_back(linRhs);
    else
      assemble_->AddLinearForm(linRhs);
  }

  void ElecPDE::DefineAvailResults() {
    
    // Electric Potential

    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->As<std::string>();

    if ( approxType == "lagrange" ) {
      // Create new resultDof object
      shared_ptr<ResultInfo> res1( new ResultInfo);
      res1->resultType = ELEC_POTENTIAL;

      // check for special subtype 
      if( subType_  != "flatShell" ) {
        shared_ptr<AnsatzFct> fct(new LagrangeFct);
        res1->dofNames = "";
        res1->unit = "V";
        res1->definedOn = ResultInfo::NODE;
        res1->entryType = ResultInfo::SCALAR;
        res1->fctType = fct;
      } else {
        
        // check number of composite materials:
        // Currently, we can handle just one electrostatic composite
        // material, as we would have a different number of electric
        // unknowns for different regions
        if( compositeMaterials_.size() > 1 ) {
          WARN("Currently the electrostatic flatShell PDE can only "
              "handle ONE type of composite material!");
        }
          
        Composite & actComp = compositeMaterials_.begin()->second;
        UInt numLaminas = actComp.thickness.GetSize();
        for( UInt i=0; i<numLaminas; i++ ) {
          std::string dofName = "ep";
          dofName += lexical_cast<std::string>(i+1);
          res1->dofNames.Push_back( dofName );
        }
        shared_ptr<AnsatzFct> fct(new LagrangeFct);
        res1->unit = "V";
        res1->fctType = fct;
        res1->definedOn = ResultInfo::ELEMENT;
        res1->entryType = ResultInfo::SCALAR;
      }
        
      results_.Push_back( res1 );
      availResults_.insert( res1 );
      
    }else if(  approxType == "spectral" ) {
      shared_ptr<ResultInfo> res1( new ResultInfo);
      UInt order = myParam_->Get("order")->As<UInt>();
      shared_ptr<SpectralFct> fct(new SpectralFct);
      res1->definedOn = ResultInfo::PFEM;
      fct->SetOrder(order);
      res1->dofNames = "";
      res1->unit = "V";
      res1->fctType = fct;

      res1->resultType = ELEC_POTENTIAL;
    	res1->entryType = ResultInfo::SCALAR;
      results_.Push_back( res1 );
      availResults_.insert( res1 );

    } else {
      // Create new resultDof object
      shared_ptr<ResultInfo> res1( new ResultInfo);
      shared_ptr<LegendreFct> fct(new LegendreFct);
      if( myParam_->Get("isIsotropic")->As<bool>() ) {
        UInt order =  myParam_->Get("order")->As<UInt>();
        fct->SetIsoOrder( order );
      } else {
        Matrix<UInt> orderMat;
        ParamTools::AsTensor<unsigned int>(myParam_->Get("anisotropic"), dim_, 1, orderMat); 
        fct->SetAnisoOrder( orderMat );
      }
      
      if( subType_ == "flatShell" ) {
        EXCEPTION( "Subtype 'flatShell' not working with Legendre functions." );
      }

      res1->resultType = ELEC_POTENTIAL;
      res1->entryType = ResultInfo::SCALAR;
      res1->dofNames = "";
      res1->unit = "V";
      res1->definedOn = ResultInfo::PFEM;
      res1->fctType = fct;
      results_.Push_back( res1 );
      availResults_.insert( res1 );
    }
    
    // Electric RHS Load
    // create new resultDof object
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ELEC_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "C";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    availResults_.insert( rhs );

    // Electric Field Intensity
    // create new resultDof object
    shared_ptr<ResultInfo> res ( new ResultInfo );
    res->resultType = ELEC_FIELD_INTENSITY;
    res->SetVectorDOFs(dim_, isaxi_);
    res->unit = "V/m";
    res->definedOn = ResultInfo::ELEMENT;
    res->entryType = ResultInfo::VECTOR;
    availResults_.insert( res );
    
    // Electric Flux Density
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = "C/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    availResults_.insert( flux );
    
    // Electric charge
    shared_ptr<ResultInfo> charge( new ResultInfo );
    charge->resultType = ELEC_CHARGE;
    charge->definedOn = ResultInfo::SURF_ELEM;
    charge->entryType = ResultInfo::SCALAR;
    charge->dofNames = "";
    charge->unit = "C";
    availResults_.insert( charge );

    // Electric energy
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_ENERGY;
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = "Ws";
    availResults_.insert ( energy );

    // Electric polarization
    shared_ptr<ResultInfo> pol( new ResultInfo );
    pol->resultType = ELEC_POLARIZATION;
    pol->definedOn = ResultInfo::ELEMENT;
    pol->entryType = ResultInfo::VECTOR;
    pol->SetVectorDOFs(dim_, isaxi_);
    pol->unit = "C/m^2";
    availResults_.insert( pol );

    // pesudo electric polarization for piezo simp
    shared_ptr<ResultInfo> pseudoPol( new ResultInfo );
    pseudoPol->resultType = ELEC_PSEUDO_POLARIZATION;
    pseudoPol->definedOn = ResultInfo::ELEMENT;
    pseudoPol->entryType = ResultInfo::SCALAR;
    pseudoPol->dofNames = "";
    pseudoPol->unit = "";
    availResults_.insert( pseudoPol );

    // grad elec potential
    shared_ptr<ResultInfo> grad_sol(new ResultInfo);
    grad_sol->resultType = GRAD_ELEC_POTENTIAL;
    grad_sol->SetVectorDOFs(dim_, isaxi_);
    grad_sol->unit = "V/m";
    grad_sol->entryType = ResultInfo::VECTOR;
    grad_sol->definedOn = ResultInfo::NODE;
    grad_sol->fctType = shared_ptr<ConstFct>(new ConstFct() ); // ??? - copy and paste! !?
    availResults_.insert( grad_sol );

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;
    
    LOG_DBG2(elecpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode elecPDENCIfaceListNode;
    elecPDENCIfaceListNode = param->GetByVal("sequenceStep", std::string("index"), sequenceStep_)
    ->Get("pdeList/electrostatic/ncInterfaceList", ParamNode::PASS);
    
    if(!elecPDENCIfaceListNode)
      return;

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    if(!domainNCIfaceListNode)
    {
      EXCEPTION("No nonmatching interfaces have been specified in domain!");
    }

    ParamNodeList pdeNCIfaceNodes;
    pdeNCIfaceNodes = elecPDENCIfaceListNode->GetList("ncInterface");

    for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
      std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->As<std::string>();

      PtrParamNode domainIfaceNode = domainNCIfaceListNode
          ->GetByVal("ncInterface", "name", pdeIfaceName, ParamNode::PASS);
      if(!domainIfaceNode)
      {
        LOG_DBG2(elecpde) << "NonMatching: Nonconforming "
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
      LOG_DBG2(elecpde) << "NonMatching: Defining new ResultDof Lagrange.";
      shared_ptr<ResultInfo> lagr ( new ResultInfo );
      lagr->resultType = LAGRANGE_MULT;
      lagr->dofNames = "l";
      lagr->fctType = results_[0]->fctType;
      lagr->definedOn = results_[0]->definedOn;
      results_.Push_back( lagr );
    } 
  }
}
