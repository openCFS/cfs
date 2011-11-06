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

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Utils/StdVector.hh"
#include "Driver/solveStepElec.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Driver/assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/hysteresis.hh"
#include "Utils/preisach.hh"
#include "Elements/fespaceH1.hh"
#include "Elements/fespaceH1Lagrange.hh"
#include "Elements/fefunction.hh"
#include "DataInOut/WriteInfo.hh"

//new integrator concept
#include "Forms/bdbInt.hh"
#include "Forms/gradientOp.hh"
#include "Forms/coefFunction.hh"
#include "Forms/buInt.hh"

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

      // Check nonLinId was already registered
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
    upLagrangeForm = true;
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
    shared_ptr<FeSpace> mySpace = feFunctions_[ELEC_POTENTIAL]->GetFeSpace();
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region node
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // ==========================================================
      // New implementation
      // ==========================================================
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);

      // --- standard real-valued stiffness integrator ---
      shared_ptr<CoefFunction > curCoef = actSDMat->GetCoefFunction(ELEC_PERMITTIVITY,tensorType,Global::REAL);


      BDBInt< GradientOperator ,FeH1,Double,Double >* stiffInt;
      stiffInt = new BDBInt<GradientOperator,FeH1,Double,Double >(curCoef,1.0 );
      //linElecInt *  linElecForm = new linElecInt( actSDMat, tensorType,
       //                                           upLagrangeForm );
      //linElecForm->SetFactor( factor );
      BiLinFormContext * stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );

      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
      stiffInt->SetFeSpace( feFunctions_[ELEC_POTENTIAL]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      biLinForms_[actRegion] = stiffInt;
    }

    // Define integrators for composite materials
    // (only for subType "flatShell")
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
         compIt++ ) {
      
      REFACTOR;
    }

    // **********************************************************************
    //   inhom. Neumann boundary condition
    // **********************************************************************
    for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) 
    {
      REFACTOR;
    }
    
    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    ///OUT FOR REFACTOR
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
      
      REFACTOR;
    }
  }

  void ElecPDE::DefineSolveStep() {
    solveStep_ = new SolveStepElec(*this);
  }

  void ElecPDE::ReadSpecialBCs( ) 
  {
     //ReadImpedances();
  }
  
  void ElecPDE::ReadImpedances( ) 
  {
    REFACTOR;
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
       feFunctions_[ELEC_POTENTIAL]->ExtractResult( res );
       break;

     case ELEC_RHS_LOAD:
       rhsFeFunctions_[ELEC_POTENTIAL]->ExtractResult( res );
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

     default:
       WARN( "Result '" << 
             SolutionTypeEnum.ToString(res->GetResultInfo()->resultType)
             << "' type not computable by electric PDE" );
     }
     
   }
   
   template <class TYPE>
   void ElecPDE::CalcElecPot( StdVector<const Elem*>& elems,
                              StdVector<LocPoint>& points,
                              Vector<TYPE>& values )  {
     Vector<TYPE> elemSol;
     Vector<Double> shape;
     values.Resize(elems.GetSize());
     values.Init();
     FeFunction<TYPE>& fct = 
           dynamic_cast<FeFunction<TYPE>& >(*(GetFeFunction(ELEC_POTENTIAL)));
     FeSpaceH1& space = dynamic_cast<FeSpaceH1&>(*(fct.GetFeSpace()));
     LocPointMapped lpm;
     for ( UInt iElem = 0; iElem < elems.GetSize(); ++iElem ) {
       
       const Elem * el = elems[iElem];
       FeH1*  ptFe = dynamic_cast<FeH1*>(space.GetFe( el->elemNum ));
       fct.GetElemSolution(elemSol, el);
       ptFe->GetShFnc(shape, points[iElem], el);
       values[iElem] = shape * elemSol;
     }
   }
   
   template <class TYPE>
   void ElecPDE::CalcElecField( StdVector<const Elem*>& elems,
                                StdVector<LocPoint>& points,
                                Vector<TYPE>& values )  {
     REFACTOR;
//<<<<<<< .mine
//     Vector<TYPE> tempE, elemSol;
//     NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
//     FeFunction<TYPE>& fct = 
//         dynamic_cast<FeFunction<TYPE>& >(*(GetFeFunction(ELEC_POTENTIAL)));
//     values.Resize( elems.GetSize() * dim_ );
//
//     // loop over elements
//     LocPointMapped lpm;
//     for ( UInt iElem = 0; iElem < elems.GetSize(); ++iElem ) {
//       const Elem * el = elems[iElem];
//       linElecInt * li = biLinForms_[el->regionId];
//       shared_ptr<ElemShapeMap> esm = ptgrid_->GetElemShapeMap( el, true );
//       BaseFE*  ptFe = fct.GetFeSpace()->GetFe( el->elemNum );
//       lpm.Set( points[iElem], esm );
//       solhelp.GetElemSolution( elemSol, el);
//       li->ApplyBMat( tempE, lpm, ptFe, elemSol);
//       // loop over dofs
//       for(UInt iDim = 0; iDim < dim_; iDim++ ) {
//         values[iElem*dim_ + iDim] = -tempE[iDim];
//       }
//     }
//=======
//     Vector<TYPE> tempE, elemSol;
//     NodeStoreSol<TYPE> & solhelp = dynamic_cast<NodeStoreSol<TYPE>&>(*sol_);
//     FeFunction<TYPE>& fct = 
//         dynamic_cast<FeFunction<TYPE>& >(*(GetFeFunction(ELEC_POTENTIAL)));
//     values.Resize( elems.GetSize() * dim_ );
//
//     // loop over elements
//     LocPointMapped lpm;
//     for ( UInt iElem = 0; iElem < elems.GetSize(); ++iElem ) {
//       const Elem * el = elems[iElem];
//       GradientOperator<FeH1,TYPE> li;
//       shared_ptr<ElemShapeMap> esm = ptgrid_->GetElemShapeMap( el, true );
//       BaseFE*  ptFe = fct.GetFeSpace()->GetFe( el->elemNum );
//       lpm.Set( points[iElem], esm );
//       solhelp.GetElemSolution( elemSol, el);
//       li.ApplyOp( tempE, lpm, ptFe, elemSol);
//       // loop over dofs
//       for(UInt iDim = 0; iDim < dim_; iDim++ ) {
//         values[iElem*dim_ + iDim] = -tempE[iDim];
//       }
//     }
//>>>>>>> .r11288
   }

   void ElecPDE::CalcField( SolutionType solType, StdVector<const Elem*>& elems,
                            StdVector<LocPoint>& points, SingleVector& values ) {

     // switch, depending on solutionType
     switch(solType) {
       case ELEC_POTENTIAL:
         if(isComplex_) {
           Vector<Complex>& temp = dynamic_cast<Vector<Complex>& >(values);
           CalcElecPot(elems, points, temp );
         } else {
           Vector<Double>& temp = dynamic_cast<Vector<Double>& >(values);
           CalcElecPot<Double>(elems, points, temp );
         }
         break;
       case ELEC_FIELD_INTENSITY:
         if(isComplex_) {
           Vector<Complex>& temp = dynamic_cast<Vector<Complex>& >(values);
           CalcElecField(elems, points, temp );
         } else {
           Vector<Double>& temp = dynamic_cast<Vector<Double>& >(values);
           CalcElecField<Double>(elems, points, temp );
         }
         break;
       default:
         EXCEPTION("Result not computable by this PDE");
     }
   }


   template <class TYPE>
   void ElecPDE::CalcElectricField( shared_ptr<BaseResult> sol ) {
     Vector<TYPE> tempE, elemSol;
     Result<TYPE> &  actSol = 
         dynamic_cast<Result<TYPE>&>(*sol);
     EntityIterator it = actSol.GetEntityList()->GetIterator();
     Vector<TYPE> & actVal = actSol.GetVector();
     actVal.Resize( actSol.GetEntityList()->GetSize() * dim_ );
     FeFunction<TYPE>& fct = 
         dynamic_cast<FeFunction<TYPE>& >(*(GetFeFunction(actSol.GetResultInfo()->resultType)));


     // loop over elements
     LocPointMapped lpm;
     for ( it.Begin(); !it.IsEnd(); it++ ) {
       const Elem * el = it.GetElem();
       GradientOperator<FeH1,TYPE> li;
       shared_ptr<ElemShapeMap> esm = ptgrid_->GetElemShapeMap( el, true );
       BaseFE*  ptFe = fct.GetFeSpace()->GetFe( it );
       lpm.Set( Elem::shapes[el->type].midPointCoord, esm );
       fct.GetElemSolution( elemSol, it.GetElem() );
       li.ApplyOp( tempE, lpm, ptFe, elemSol);
       // loop over dofs
       for(UInt iDim = 0; iDim < dim_; iDim++ ) {
         actVal[it.GetPos()*dim_ + iDim] = -tempE[iDim];
       }
     }
   }
   
  
  template <class TYPE>
  void ElecPDE::CalcElectricFluxDensity( shared_ptr<BaseResult> sol )
  {
    REFACTOR;
  }
  
  


  void ElecPDE::CalcPolarizationField( shared_ptr<BaseResult> res )
  {
    REFACTOR;
  }



  template <class TYPE>
  void ElecPDE::CalcCharges( shared_ptr<BaseResult> res ) {
    REFACTOR;
  }

  template <class TYPE>
  void ElecPDE::CalcEnergy( shared_ptr<BaseResult> vals )
  {
    Matrix<Double> elemmat;  
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
      BiLinearForm * bilinear_stiff = biLinForms_[regionIt.GetRegion()];

      EntityIterator elemIt = actSDList.GetIterator();

      // Loop over elements
      energy = 0;
      Vector<TYPE> elpot;
      for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {

        //COMMENTED OUT DUE TO REFACTORING
        //bilinear_stiff->SetFactor(factor);
        bilinear_stiff->CalcElementMatrix( elemmat, elemIt, elemIt );

        feFunctions_[ELEC_POTENTIAL]->GetEntitySolution(elpot, elemIt);
        help =  elemmat * elpot;
        energy += 0.5 * (help * elpot);

        LOG_DBG3(elecpde) << "CalcEnergy: elem=" << elemIt.GetElem()->elemNum << " mat=" << elemmat.ToString();
        LOG_DBG3(elecpde) << "CalcEnergy: elemIt=" << elemIt.GetElem()->elemNum << " elpot=" << elpot.ToString();
        LOG_DBG3(elecpde) << "CalcEnergy: elemIt=" << elemIt.GetElem()->elemNum << " sum energy -> " << energy;
      }  
      actVal[regionIt.GetPos()] = energy;
    }
  }


  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void ElecPDE::InitCoupling(PDECoupling * Coupling)
  {
  REFACTOR;
  }
  


  void ElecPDE::CalcOutputCoupling()
  {
    REFACTOR;
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
    REFACTOR;
  }

  void ElecPDE::DefineAvailResults() {
    
    // Electric Potential
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = ELEC_POTENTIAL;
    
    // check for special subtype 
    if( subType_  != "flatShell" ) {
      res1->dofNames = "";
      res1->unit = "V";
      res1->definedOn = ResultInfo::NODE;
      res1->entryType = ResultInfo::SCALAR;
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
      res1->unit = "V";
      res1->definedOn = ResultInfo::ELEMENT;
      res1->entryType = ResultInfo::SCALAR;
    }
    feFunctions_[ELEC_POTENTIAL]->SetResultInfo(res1);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    

    // Electric RHS Load
    // create new resultDof object
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ELEC_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "C";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    availResults_.insert( rhs );
    postProcResults_[ELEC_RHS_LOAD] = ELEC_POTENTIAL;
    
    // Electric Field Intensity
    // create new resultDof object
    shared_ptr<ResultInfo> res ( new ResultInfo );
    res->resultType = ELEC_FIELD_INTENSITY;
    res->SetVectorDOFs(dim_, isaxi_);
    res->unit = "V/m";
    res->definedOn = ResultInfo::ELEMENT;
    res->entryType = ResultInfo::VECTOR;
    availResults_.insert( res );
    postProcResults_[ELEC_FIELD_INTENSITY] = ELEC_POTENTIAL;
    
    // Electric Flux Density
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = "C/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    availResults_.insert( flux );
    postProcResults_[ELEC_FLUX_DENSITY] = ELEC_POTENTIAL;
    
    // Electric charge
    shared_ptr<ResultInfo> charge( new ResultInfo );
    charge->resultType = ELEC_CHARGE;
    charge->definedOn = ResultInfo::SURF_ELEM;
    charge->entryType = ResultInfo::SCALAR;
    charge->dofNames = "";
    charge->unit = "C";
    availResults_.insert( charge );
    postProcResults_[ELEC_CHARGE] = ELEC_POTENTIAL;

    // Electric energy
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_ENERGY;
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = "Ws";
    availResults_.insert ( energy );
    postProcResults_[ELEC_ENERGY] = ELEC_POTENTIAL;
    
    // Electric polarization
    shared_ptr<ResultInfo> pol( new ResultInfo );
    pol->resultType = ELEC_POLARIZATION;
    pol->definedOn = ResultInfo::ELEMENT;
    pol->entryType = ResultInfo::VECTOR;
    pol->SetVectorDOFs(dim_, isaxi_);
    pol->unit = "C/m^2";
    availResults_.insert( pol );
    postProcResults_[ELEC_POLARIZATION] = ELEC_POTENTIAL;
    
    // pesudo electric polarization for piezo simp
    shared_ptr<ResultInfo> pseudoPol( new ResultInfo );
    pseudoPol->resultType = ELEC_PSEUDO_POLARIZATION;
    pseudoPol->definedOn = ResultInfo::ELEMENT;
    pseudoPol->entryType = ResultInfo::SCALAR;
    pseudoPol->dofNames = "";
    pseudoPol->unit = "";
    availResults_.insert( pseudoPol );
    postProcResults_[ELEC_PSEUDO_POLARIZATION] = ELEC_POTENTIAL;

    // grad elec potential
    shared_ptr<ResultInfo> grad_sol(new ResultInfo);
    grad_sol->resultType = GRAD_ELEC_POTENTIAL;
    grad_sol->SetVectorDOFs(dim_, isaxi_);
    grad_sol->unit = "V/m";
    grad_sol->entryType = ResultInfo::VECTOR;
    grad_sol->definedOn = ResultInfo::NODE;
    availResults_.insert( grad_sol );
    postProcResults_[GRAD_ELEC_POTENTIAL] = ELEC_POTENTIAL;

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
      lagr->definedOn = results_[0]->definedOn;
      results_.Push_back( lagr );
    } 
  }
  
  std::map<SolutionType, shared_ptr<FeSpace> >
 ElecPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("elecPotential");
      crSpaces[ELEC_POTENTIAL] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[ELEC_POTENTIAL]->Init();
    }else{
      EXCEPTION("The formulation " << formulation << "of electric PDE is not known!");
    }
    return crSpaces;
  }
}
