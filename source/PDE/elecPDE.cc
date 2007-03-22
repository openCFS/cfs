// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>

#include "elecPDE.hh"

#include "Forms/linElecInt.hh"
#include "Forms/nLinElecHystInt.hh"
#include "Forms/nLinElecMaterial.hh"
#include "Forms/elecchargeop.hh"
#include "Forms/laplaceInt.hh"
#include "Forms/FlatShellElecInt.hh"
#include "Forms/nonConformingInt.hh"
#include "Forms/massInt.hh"
#include "Forms/gradfieldop.hh"
#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
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

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField {

  DECLARE_LOG(elecpde)
  DEFINE_LOG(elecpde, "pde.electrostatic")

  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE( Grid* aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    ENTER_FCN( "ElecPDE::ElecPDE" );

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

    // Check the subtype of the problem
    paramNode->Get("subType", subType_);
  }
  

  void ElecPDE::InitNonLin() {
    ENTER_FCN( "ElecPDE::InitNonLin" );

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
      nonLinIdType_[actId] = actType;

      // check type
      if( actType == HYSTERESIS ) {
        isHysteresis_ = true;
      }

      if( actType == MATERIAL ) {
        nonLin_ = true;
        nonLinMaterial_ = true;
      }
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
    Info->PrintF( pdename_, "\n" );

  }

  void ElecPDE::DefineIntegrators() {

    ENTER_FCN( "ElecPDE::DefineIntegerators" );
    
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
    if ( isPiezoCoupled_ == true )
      factor = "-1.0";  

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region node
      std::string regionName = ptgrid_->RegionIdToName( actRegion );

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

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
        nlForm->SetFactor( GenStr(factor) );

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
        linElecInt *  linElecForm = new linElecInt( actSDMat, tensorType,
                                                    upLagrangeForm );
        linElecForm->SetFactor( factor );

        BiLinFormContext * stiffIntDescr = 
          new BiLinFormContext(linElecForm, STIFFNESS );

        stiffIntDescr->SetPtPdes(this, this);
        stiffIntDescr->SetResults( results_[0], results_[0],
                                   actSDList, actSDList );
	
        assemble_->AddBiLinearForm( stiffIntDescr );
      
      
        // --- check for complex valued material parameter ---
        // check for piezo-coupling and 
        if( isPiezoCoupled_ ) {
          ParamNode * dataTypeNode = 
            param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
            ->Get("couplingList")->Get("direct")->Get("piezoDirect")
            ->Get("materialDataType", false );
          bool isImag = false;

          if( dataTypeNode ) 
            isImag = dataTypeNode->Get("type")->AsString() 
              == "imagMaterialParameter";

          if( isImag ) {
            DataType matType = IMAG; 
            
            linElecInt * linElecFormC  = new 
              linElecInt( materials_[actRegion], tensorType,
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
      
      FlatShellElecInt * compElecInt = new FlatShellElecInt( composite );
      BiLinFormContext * stiffContext = 
        new BiLinFormContext( compElecInt, STIFFNESS);
      stiffContext->SetPtPdes( this, this );
      stiffContext->SetResults( results_[0], results_[0],
                                actSDList, actSDList );
      assemble_->AddBiLinearForm( stiffContext );
      eqnMap_->AddResult( *results_[0], actSDList );
      
    }

    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    for( UInt i = 0; i < ncIFaces_.GetSize(); i++ ) {
      
      // get regionId of Lagrangian surface
      StdVector<std::string> keyVec, attrVec, valVec;
      std::string slaveSide;
      std::string ncIfaceName = ptgrid_->RegionIdToName(ncIFaces_[i]);

      ParamNode* ncIfaceListNode;
      ncIfaceListNode = param->Get("domain")->Get("ncInterfaceList");

      slaveSide = ncIfaceListNode->Get("ncInterface",
                                       "name",
                                       ncIfaceName)->Get("slaveSide")->AsString();

      // Part 1: Define integrator M(Psi, Lambda) on
      //         non-conforming interface
      LOG_DBG2(elecpde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->RegionIdToName(ncIFaces_[i]) << "'.";
      shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
      actNcList->SetRegion( ncIFaces_[i] );
      
      NonConformingInt * ncInt = 
        new NonConformingInt( 1, isaxi_ );

      NcBiLinFormContext * stiffIntDescr = 
	new NcBiLinFormContext( ncInt , STIFFNESS );

      // Force assembling of M(Psi, Lambda)^T
      stiffIntDescr->SetCounterPart( true );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[1],
                                 actNcList, actNcList );
      
      assemble_->AddBiLinearForm( stiffIntDescr );


      // Part 2: Define integrator D(Psi, Lambda) on
      //         Lagrangian surface
      LOG_DBG2(elecpde) << "NonMatching: Defining mass integrator"
                        << " for D on interface '"
                        << ptgrid_->RegionIdToName(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );      
      actSDList->SetRegion( ptgrid_->RegionNameToId( slaveSide ) );

      // D(Psi, Lambda) has the form of a standard mass
      // integrator with factor 1.0
      MassInt * dMatInt = new MassInt( 1.0, 1, isaxi_ );
      BiLinFormContext * dMatContext = 
        new BiLinFormContext( dMatInt, STIFFNESS );

      // Force assembling of D(Psi, Lambda)^T
      dMatContext->SetCounterPart( true );
      dMatContext->SetPtPdes( this, this );
      dMatContext->SetResults( results_[0], results_[1],
                               actSDList, actSDList );
      
      assemble_->AddBiLinearForm( dMatContext );

      // Give result LAGRANGE_MULT to equation numbering class
      eqnMap_->AddResult( *results_[1], actSDList );
    }
  }

  void ElecPDE::DefineSolveStep() {
    ENTER_FCN( "ElecPDE::DefineSolveStep" );
    solveStep_ = new SolveStepElec(*this);
  }

  void ElecPDE::ReadSpecialBCs( ) {
    ENTER_FCN( "ElecPDE::ReadSpecialBCs" );
  
    // Get values from parameter file
    StdVector<ParamNode*> impedNodes = 
      myParam_->Get("bcsAndLoads")->GetList("impedance");
   
    // make sure we are in harmonic mode
    if( impedNodes.GetSize() > 0 && analysistype_ != HARMONIC) {
      EXCEPTION( "Impedances are only available in harmonic simulation "
                 << "at the moment!" );
    }

    std::string node1, node2;
    Double resistance, capacitance, inductance;
    StdVector<UInt> nodes;
    
    // resize impedances and fill in data
    std::ostringstream out;
    impedances_.Resize( impedNodes.GetSize() );
    for( UInt i = 0; i < impedNodes.GetSize(); i++ ) {
      
      // get data from node
      node1 = impedNodes[i]->Get("node1")->AsString();
      node2 = impedNodes[i]->Get("node2")->AsString();
      resistance = impedNodes[i]->Get("resistance")->AsDouble();
      inductance = impedNodes[i]->Get("inductance")->AsDouble();
      capacitance = impedNodes[i]->Get("capacitance")->AsDouble();
      
      // fill data for impedance node
      ptgrid_->GetNodesByName( nodes, node1 );
      impedances_[i].node1 = nodes[0];
       ptgrid_->GetNodesByName( nodes, node2 );
      impedances_[i].node2 = nodes[0];
      impedances_[i].resistance = resistance;
      impedances_[i].capacitance = capacitance;
      impedances_[i].inductance = inductance;

      out.clear();
      out.str("");
      out << "Impedance defined between node " << impedances_[i]. node1
          << " and node " << impedances_[i].node2  << " with:\n";
      Info->PrintF( pdename_, out.str().c_str() );
      out.clear();
      out.str("");
      out << "\tR = " << impedances_[i].resistance 
          << "\tL = " << impedances_[i].inductance         
          << "\tC = " << impedances_[i].capacitance << std::endl;
      Info->PrintF( pdename_, out.str().c_str() );
      Info->PrintF( pdename_, "\n" );
    }
  }
  
  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  // ***************
  //   CalcResults
  // ***************
   void ElecPDE::CalcResults( shared_ptr<BaseResult> res ) {
     ENTER_FCN( "ElecPDE::CalcResults" );
     
     switch (res->GetResultInfo()->resultType ) {
     case ELEC_POTENTIAL:
       if( isComplex_ ) {
         ExtractResult<Complex>( res, sol_ );
       } else {
         ExtractResult<Double>( res, sol_ );
       }
       break;
       
     case ELEC_FIELD_INTENSITY:
       if( isComplex_ ) {
         CalcElectricField<Complex>( res );
       } else {
         CalcElectricField<Double>( res );
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
       
     default:
       Warning( "Resulttype not computable by electric PDE",
              __FILE__, __LINE__ );
     }
     
   }


  template <class TYPE>
  void ElecPDE::CalcElectricField( shared_ptr<BaseResult> sol )
  {
    ENTER_FCN( "ElecPDE::PostProcess" );
    Vector<Double> lCoord;
    Vector<TYPE> tempE;
    NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
    GradientFieldOp<TYPE> * FieldOp = 
      new GradientFieldOp<TYPE>(ptgrid_, this, 
                                eqnMap_, solhelp, 
                                ELEC_POTENTIAL, results_[0],isaxi_);
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


  void ElecPDE::CalcPolarizationField( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "ElecPDE::CalcPolarizationField" );

    //we assume, that the actual solution is stored in sol_!
    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);

    GradientFieldOp<Double> * FieldOp = 
      new GradientFieldOp<Double>(ptgrid_,this, eqnMap_,
 				  *solhelp, ELEC_POTENTIAL, 
 				  results_[0], isaxi_);
  
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
        vecP[dirP] = Pval;
        for( UInt iDof = 0; iDof < dim_; iDof++ ) {
          actVal[it.GetPos()*dim_ + iDof] = vecP[iDof];
        }
      }  
    } else {
      *warning << "No hysteresis defined for list '" <<
        actSol.GetEntityList()->GetName() << std::endl; 
      Warning( __FILE__, __LINE__ );}
  }



  template <class TYPE>
  void ElecPDE::CalcCharges( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "ElecPDE::CalcCharges" );


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
                                ELEC_POTENTIAL, results_[0], isaxi_);
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
      // calculation
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
      myMat->GetScalar(permittivity,ELEC_PERMITTIVITY,REAL);

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
    ENTER_FCN( "ElecPDE::CalcEnergy" );

    Matrix<Double> elemmat;  
    SubTensorType tensorType;

    if ( dim_ == 3 ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == TRUE ) {
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
      
      linElecInt * bilinear_stiff = 
        new linElecInt(materials_[regionIt.GetRegion()],tensorType);

      // Loop over elements
      energy = 0;
      for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        
        bilinear_stiff->SetFactor(factor);
        bilinear_stiff->CalcElementMatrix( elemmat, elemIt, elemIt );
        
        Vector<TYPE> elpot;
        sol_->GetElemSolution(elpot, elemIt );
        help =  elemmat * elpot;
        energy += 0.5 * (help * elpot);
        
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
    ENTER_FCN( "ElecPDE::InitCoupling" );
  
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    
    StdVector<std::string> * nRegions;
    StdVector<RegionIdType> nRegionIds;
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();
    
    nonLin_ = false;

    // Initialization of coupling helper arrays
    std::string quantity;
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
          ptgrid_->RegionNameToId(nRegionIds,*nRegions);
          ForceOp_->Setup(nRegionIds, *couplingnodes);
        }
      
        else if ( ptCoupling_->GetOutputQuantity(actCoupling)
                  == ELEC_INTERFACE_FORCE) {
          EXCEPTION( "Currently ELEC_INTERFACE_FORCE not supported" );
        }

        else {
          Enum2String(ptCoupling_->GetOutputQuantity(actCoupling), quantity);
          EXCEPTION( "Coupling " << quantity <<  " not known! ");
        }
      
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);
      
      
      } // end for (actNode)
    } 
  }
  


  void ElecPDE::CalcOutputCoupling()
  {
    ENTER_FCN( "ElecPDE::CalcOutputCoupling" );

    SolutionType quantity;
    StdVector<UInt> * couplingNodes     = NULL;
    CFSVector * values = 0;
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
    ENTER_FCN( "ElecPDE::HasOutput" );
  
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
    ENTER_FCN( "ElecPDE::SetPiezoCoupling" );
  
    isPiezoCoupled_ = true;

  }

  void ElecPDE::DefineAvailResults() {
    ENTER_FCN( "ElecPDE::DefineAvailResults" );
    
    // Determine vectorial dofNames
    StdVector<std::string> vecDofNames;
    if( isaxi_) {
      vecDofNames = "r", "z";
    } else if (dim_ == 2) {
      vecDofNames = "x", "y";
    } else {
      vecDofNames = "x", "y", "z";
    }
    
    // Electric Potential

    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->AsString();

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
        // Determine number of laminas for setting number of dofs
        StdVector<ParamNode*>  laminas = param->Get("domain")
          ->GetList("composite");
        
        if (laminas.GetSize() == 0) {
          res1->dofNames.Push_back("ep");
        } else {
          for( UInt i=0; i<laminas[0]->Count("lamina"); i++ ) {
            std::string dofName = "ep";
            dofName += GenStr(i+1);
            res1->dofNames.Push_back( dofName );
          }
        }
        shared_ptr<AnsatzFct> fct(new LagrangeFct);
        res1->unit = "V";
        res1->fctType = fct;
        res1->definedOn = ResultInfo::ELEMENT;
        res1->entryType = ResultInfo::SCALAR;
      }
        
      results_.Push_back( res1 );
      availResults_.insert( res1 );
      
    } else {
      // Create new resultDof object
      shared_ptr<ResultInfo> res1( new ResultInfo);
      shared_ptr<LegendreFct> fct(new LegendreFct);
      if( myParam_->Get("isIsotropic")->AsBool() ) {
        UInt order =  myParam_->Get("order")->AsUInt();
        fct->SetIsoOrder( order );
      } else {
        Matrix<UInt> orderMat;
        myParam_->GetDim1xDim2Tensor("anisotropic", dim_, 1, orderMat);
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

    // Electric Field Intensity
    // create new resultDof object
    shared_ptr<ResultInfo> res ( new ResultInfo );
    res->resultType = ELEC_FIELD_INTENSITY;
    res->dofNames = vecDofNames;
    res->unit = "V/m";
    res->definedOn = ResultInfo::ELEMENT;
    res->entryType = ResultInfo::VECTOR;
    availResults_.insert( res );
    
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
    pol->dofNames = vecDofNames;
    pol->unit = "C/m^2";
    availResults_.insert( pol );
   
    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;
    
    LOG_DBG2(elecpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    ParamNode* domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", false);

    if(domainNCIfaceListNode)
    {

        StdVector<ParamNode*> pdeNCIfaceNodes = 
            param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
            ->Get("pdeList")->Get("electrostatic")->Get("ncInterfaceList")
            ->GetList("ncInterface");
    
        for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
            std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->AsString();
        
            ParamNode* domainIfaceNode = domainNCIfaceListNode->Get("ncInterface",
                                                                    "name",
                                                                    pdeIfaceName,
                                                                    false);
            if(!domainIfaceNode)
            {
                LOG_DBG2(elecpde) << "NonMatching: Nonconforming "
                                  << "interface '" << ncIfaceNames[i]
                                  << "' does not exist in domain.";
            
                EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
            }

            ncIfaceNamesForPDE.Push_back(pdeIfaceName);
        }
        ptgrid_->RegionNameToId( ncIfaceIds, ncIfaceNamesForPDE );
    
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
}
