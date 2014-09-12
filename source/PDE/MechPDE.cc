// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MechPDE.hh"

#include <sstream>
#include <iomanip>
#include <set>

#include "General/defs.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Driver/Assemble.hh"


// include elements
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ICModesInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormalTrans.hh"
#include "Forms/Operators/StrainOperator.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"
#include "Forms/BiLinForms/SingleEntryBiLinInt.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/EigenFrequencyDriver.hh"

#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"

namespace CoupledField {

DECLARE_LOG(mechpde)
DEFINE_LOG(mechpde, "mechpde")

MechPDE::MechPDE(Grid * aptgrid, PtrParamNode paramNode,PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    pdename_          = "mechanic";
    pdematerialclass_ = MECHANIC;

    nonLin_        = false;
    nonLinMaterial_= false;

    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;

    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    myParam_->GetValue("subType", subType_ );

    std::string probGeo;
    domain_->GetParamRoot()->Get("domain")->GetValue("geometryType", probGeo );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      stressDim_ = 6;
      tensorType_ = FULL;
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      stressDim_ = 4;
      tensorType_ = AXI;
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRAIN;
    }

    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRESS;
    }

    else if ( subType_ == "flatShell" ) {
      stressDim_ = 3;
    }
    else {
      EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
                 <<  pdename_ <<  "' does not fit to problem  geometry '"
                 << probGeo << "'"; );
    }
    
    // Sanity check: 3D can only be computed if 3D elements are present/
    if( subType_ == "3d" && ptGrid_->GetNumElemOfDim(3) == 0 ) {
      EXCEPTION("Can not calculate 3D mechanics without 3D elements in the grid!");
    }

    // thermal stress coefFunction
    thermalStress_.reset(new CoefFunctionMulti(CoefFunction::VECTOR,
                         stressDim_, 1, isComplex_, true));

    // thermal strain coefFunction
    thermalStrain_.reset(new CoefFunctionMulti(CoefFunction::VECTOR,
                         stressDim_, 1, isComplex_, true));
}


  MechPDE::~MechPDE()
  {

  }

  void MechPDE::ReadDampingInformation( ) {

//    bool identical = true; // i.e. same type of damping for all regions

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
          dampNodes[i]->GetValue( "ratioDeltaF", actRaylDamp->ratioDeltaF, ParamNode::PASS );
          dampNodes[i]->GetValue( "adjustDamping", actRaylDamp->adjustDamping, ParamNode::PASS );
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

//    if( regionNodes.GetSize() > 0 ) {
//      Info->PrintF( pdename_, "Damping in following region(s)\n" );
//    }

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
        
        // Check if analysis is harmonic
        bool isHarmonic = analysistype_ == BasePDE::HARMONIC;

        // Compute Rayleigh damping parameters
        materials_[actRegionId]->
         ComputeRayleighDamping( actRayl.alpha, actRayl.beta,
                                 dampFreq, actRayl.ratioDeltaF, 
                                 actRayl.adjustDamping, isHarmonic );
        regionRaylDamping_[actRegionId] = actRayl;

        PtrParamNode in = infoNode_->Get(ParamNode::HEADER)
            ->GetByVal("region", "name", ptGrid_->GetRegion().ToString(actRegionId));
        in->Get("alpha_M")->SetValue(actRayl.alpha);
        in->Get("alpha_K")->SetValue(actRayl.beta);
      }
    }

    // Check, if all entries are identical
    for ( UInt i = 1; i < dampingList_.size(); i++ ) {
      if ( dampingList_[regions_[i-1]] != dampingList_[regions_[i]] ) {
//        identical = false;
        break;
      }
    }
  }

  void MechPDE::ReadSoftening() {

    // Check if softeningList node is present
    std::string type, id;
    std::map<std::string, std::string> idSoftTypeMap;
    PtrParamNode softListNode = myParam_->Get("softeningList", ParamNode::PASS );
    PtrParamNode softInfo;
    if( softListNode ) {

      // Get child elements and read data
      ParamNodeList softNodes = softListNode->GetChildren();
      for( UInt i = 0; i < softNodes.GetSize(); i++ ) {
        type = softNodes[i]->GetName();
        softNodes[i]->GetValue( "id", id );
        idSoftTypeMap[id] = type;
      }

      if( softNodes.GetSize() ) {
        softInfo = infoNode_->Get("softeningList");
      }
    }

    // Now iterate over all regions and check for softening type
    ParamNodeList regionNodes =
        myParam_->Get("regionList")->GetList("region");
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get region Name
      std::string regionName = regionNodes[i]->Get("name")->As<std::string>();
      RegionIdType regionId = ptGrid_->GetRegion().Parse( regionName );
      
      // get softeningId of current region
      id = "";
      regionNodes[i]->GetValue( "softeningId", id, ParamNode::PASS );
      if (id == "") continue;

      // try to find related softening definition
      if( idSoftTypeMap.find( id) == idSoftTypeMap.end() ) {
        EXCEPTION( "Softening with id '" << id << "' for region '"
                   << regionName << "' could not be found in "
                   << "softeningList " );
      }
      // assign correct softening type to current region and map to log
      regionSoftening_[regionId] = idSoftTypeMap[id];
      PtrParamNode regionNode = softInfo->Get("region",ParamNode::APPEND);
      regionNode->Get("name")->SetValue(regionName);
      regionNode->Get("type")->SetValue(idSoftTypeMap[id]);
    }
  }

  void MechPDE::InitNonLin()
  {

    nonLin_ = false;
  }


  void MechPDE::DefineIntegrators() {

    // =======================================
    //  Get information about softening types
    // =======================================
    ReadSoftening();
    
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    
    // determine if we do bloch eigenfrequency analysis
    bool do_bloch = domain->GetDriver()->DoBlochModeEigenfrequency();

    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ )
    {
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // --------------------------
      //  Set region approximation
      // --------------------------
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      // ====================================================================
      //  Standard Linear Stiffness
      // ====================================================================
      if( !nonLin_ )
      {
        // either complex material or bloch mode with complex B-matrices
        bool complex = do_bloch | complexMatData_[actRegion];
        BaseBDBInt * stiffInt =  GetStiffIntegrator( actSDMat, actRegion, complex );
        stiffInt->SetName("LinElastInt");
        stiffInt->SetFeSpace( mySpace);
        
        BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );
        
        stiffIntDescr->SetEntities( actSDList, actSDList );
        stiffIntDescr->SetFeFunctions( myFct, myFct );
        
        //check for damping
        if ( dampingList_[actRegion] == RAYLEIGH ) {
          RaylDampingData & actDamp = (regionRaylDamping_[actRegion]);
          stiffIntDescr->SetSecDestMat(DAMPING, actDamp.beta );
        }
        
        assemble_->AddBiLinearForm( stiffIntDescr );
        
        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_[actRegion] = stiffInt;
        
      }
      
      // ====================================================================
      //  Standard Mass Integrator
      // ====================================================================
      PtrCoefFct densCoeff = actSDMat->GetScalCoefFnc( DENSITY,Global::REAL );

      // BLOCH TODO
      // PtrCoefFct densCoeff = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(density));

      BaseBDBInt *massInt = NULL;

      // complex mass integrator due to complex bloch stiffness matrix

      if(dim_== 2 && do_bloch)
        massInt = new BBInt<Complex, Complex>(new IdentityOperator<FeH1,2,2>(), densCoeff, 1.0);
      if(dim_== 2 && !do_bloch)
        massInt = new BBInt<>(new IdentityOperator<FeH1,2,2>(), densCoeff, 1.0);
      if(dim_== 3 && do_bloch)
        massInt = new BBInt<Complex, Complex>(new IdentityOperator<FeH1,3,3>(), densCoeff, 1.0);
      if(dim_== 3 && !do_bloch)
        massInt = new BBInt<>(new IdentityOperator<FeH1,3,3>(), densCoeff, 1.0);

      massInt->SetName("MassInt");
      massInt->SetFeSpace( mySpace );

      BiLinFormContext *massContext =  new BiLinFormContext( massInt, MASS );
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( myFct, myFct );
      
      // Check for damping (mass part)
      if ( dampingList_[actRegion] == RAYLEIGH ) {
        RaylDampingData & actDamp = regionRaylDamping_[actRegion];
        massContext->SetSecDestMat( DAMPING, actDamp.alpha );
      }
      
      // Important: Add mass-integrator to global list, as we need them later
      // for calculation of postprocessing results
      massInts_[actRegion] = massInt;
      assemble_->AddBiLinearForm( massContext );

      
      // in the end, at the region to the feFunction      // to be implemented
      myFct->AddEntityList( actSDList );
    }
    
    // ====================================================================
    //  Concentrated Mechanical Network Elements
    // ====================================================================
    DefineConcentratedElems();
  }
  
  void MechPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        if (dim_ == 2)
          DefineMortarCoupling<2,2>(MECH_DISPLACEMENT, *ncIt);
        else
          DefineMortarCoupling<3,3>(MECH_DISPLACEMENT, *ncIt);
        break;
      case NC_NITSCHE:
      {
        MortarInterface * ncIf = dynamic_cast<MortarInterface*>(ptGrid_
                                    ->GetNcInterface(ncIt->interfaceId).get());
        assert(ncIf);
        
        //check for softening
        bool icModes = false;
        if ( regionSoftening_[ncIf->GetMasterVolRegion()] == "icModesTW" ||
             regionSoftening_[ncIf->GetSlaveVolRegion()]  == "icModesTW" )
               icModes = true;

        if(dim_ == 2)
          DefineNitscheCoupling<2,2>(MECH_DISPLACEMENT, *ncIt, icModes);
        else
          DefineNitscheCoupling<3,3>(MECH_DISPLACEMENT, *ncIt, icModes);

      }
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }
  
  void MechPDE::DefineConcentratedElems() {

     // Get FESpace and FeFunction of mechanical displacement
     shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
     shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
     
      // try to get bcsAndLoads node
      PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
      if( !bcNode )
        return;

      // fetch parameter node specifying spring
      ParamNodeList springNodes =
          bcNode->GetList("concentratedElem");

      // Iterate over all springs
      std::string name, dofName;
      std::string massVal, dampVal, stiffVal;
      for( UInt i = 0; i < springNodes.GetSize(); i++ ) {

        // get data from node
        springNodes[i]->GetValue( "name", name );
        springNodes[i]->GetValue( "dof", dofName );
        springNodes[i]->GetValue( "massValue", massVal );
        springNodes[i]->GetValue( "dampingValue", dampVal );
        springNodes[i]->GetValue( "stiffnessValue", stiffVal );

        UInt dof = results_[0]->GetDofIndex( dofName );

        shared_ptr<EntityList> nodes = 
            ptGrid_->GetEntityList(EntityList::NODE_LIST, name);
        UInt numNodes = nodes->GetSize();

        // Ensure, that only lists with 1 or 2 nodes are present
        if( numNodes > 2 ) {
          WARN( "Concentrated mechanical element on '" 
                 << name << "' is omitted, as it consists of more than "
                 << "2 nodes!"; );
          continue;
        }

        StdVector<FEMatrixType> matrices;
        matrices = STIFFNESS, MASS, DAMPING;
        StdVector<std::string> vals;
        vals = stiffVal, massVal, dampVal;

        if( numNodes == 1 ) {
          // ============================
          //  POINT CONCENTRATED ELEMENT
          // ============================
          
          for( UInt i = 0; i < matrices.GetSize(); ++i ) {
            // if value is zero, just continue
            if( vals[i] == "0" || vals[i] == "0.0" ) {
              continue;
            }
            
            SingleEntryBiLinInt * myInt = new SingleEntryBiLinInt( dim_, vals[i], dof,
                                                                   mp_ );
            BiLinFormContext * intCtx =
                new BiLinFormContext( myInt, matrices[i] );
            intCtx->SetEntities( nodes, nodes );
            intCtx->SetFeFunctions( myFct, myFct );
            
            assemble_->AddBiLinearForm( intCtx );
          } // loop over mass/stiffness/damp

        } else {
          // =================================
          //  PAIR-WISE CONCENTRATED ELEMENTS
          // =================================

          // extract both nodes and put them into a new node list
          shared_ptr<NodeList> node1(new NodeList(ptGrid_));
          shared_ptr<NodeList> node2(new NodeList(ptGrid_));
          StdVector<UInt> tmp(1);
          EntityIterator it = nodes->GetIterator();
          tmp[0] = it.GetNode();
          node1->SetNodes( tmp );
          it++;
          tmp[0] = it.GetNode();
          node2->SetNodes( tmp );

          // loop over stiffness, mass and dampingvalue
          std::string diagVal, offDiagVal;
          for( UInt i = 0; i < matrices.GetSize(); ++ i ) {
            // if value is zero, just continue
            if( vals[i] == "0" || vals[i] == "0.0" ) {
              continue;
            }
            if( matrices[i] == STIFFNESS ) {
              diagVal = "1.0 * (" + vals[i] + ")";
              offDiagVal = "-1.0 * (" + vals[i] + ")";
            } else {
              diagVal = "(" + vals[i] + ") * 1.0 / 3.0";
              offDiagVal = "(" + vals[i] + ") * 1.0 / 6.0";
            }

            // a) diagonal entries
            SingleEntryBiLinInt * diagInt1 = new SingleEntryBiLinInt( dim_, diagVal, dof, mp_);
            BiLinFormContext * diagCtx1 =
                new BiLinFormContext( diagInt1, matrices[i] );
            diagCtx1->SetEntities( node1, node1 );
            diagCtx1->SetFeFunctions( myFct, myFct );
            assemble_->AddBiLinearForm( diagCtx1 );

            SingleEntryBiLinInt * diagInt2 = new SingleEntryBiLinInt( dim_, diagVal, dof, mp_);
            BiLinFormContext * diagCtx2 =
                new BiLinFormContext( diagInt2, matrices[i] );
            diagCtx2->SetEntities( node2, node2 );
            diagCtx2 ->SetFeFunctions( myFct, myFct );
            assemble_->AddBiLinearForm( diagCtx2 );

            // b) off-diagonal entries
            SingleEntryBiLinInt * offDiagInt = 
                new SingleEntryBiLinInt( dim_, offDiagVal, dof, mp_);
            BiLinFormContext * offDiagCtx = 
                new BiLinFormContext( offDiagInt, matrices[i] );
            offDiagCtx->SetEntities( node1, node2 );
            offDiagCtx->SetFeFunctions( myFct, myFct );
            offDiagCtx->SetCounterPart( true );
            assemble_->AddBiLinearForm( offDiagCtx );
          } // loop matrix types
        } // if: 1 / 2 nodes

      } // loop concentrated elements

    }
  
  void MechPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(mechpde) << "Defining rhs load integrators for mechanic PDE";
    
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;

    /* Flag, if coefficient function lives on updated geometry (updated 
     * Lagrangian formulation). For analytically prescribed values
     * (pressure, force, stress), this is in general not the case. 
     * In the coupled case, the magnetic force density however could live
     * on an updated geometry. In this case, the ReadRhsExcitation()
     * method will set the coefUpdateGeo flag to true, so the RHS-
     * integrator also has to use this value.
     * 
     */
    bool coefUpdateGeo = false;
    
    // ========================
    //  FORCES (volume, nodal)
    // ========================
    LOG_DBG(mechpde) << "Reading forces";
    ReadRhsExcitation( "force", dispDofNames, ResultInfo::VECTOR, isComplex_,
                       ent, coef, coefUpdateGeo );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      
      // In case of a total force, we can not have a spatial dependency
      if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
        EXCEPTION("Total forces must not be spatial dependent");
      }
              
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        
        // --------------
        //  Nodal Forces 
        // --------------
        UInt numNodes = ent[i]->GetSize();
        // If more than one node is defined, we divide the total force by the number
        // of nodes to ensure that the total force is applied, independent of the 
        // number of nodes
        if( numNodes > 1 ) {
          Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
          coef[i] = CoefFunction::Generate(mp_, part, 
                                           CoefXprVecScalOp(mp_, coef[i], 
                    boost::lexical_cast<std::string>(numNodes), CoefXpr::OP_DIV) );
        }
        
        lin = new SingleEntryInt(coef[i]);
        lin->SetName("NodalForceInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        
      } else {
        
        // -------------------------
        //  Surface / Volume Forces 
        // -------------------------
        EXCEPTION("Not yet implemented")
    
        // Same issue here as above: We need to "divide" the total force by the
        // area / volume to get the force density.
      }
    } // for
    
    // ===============
    //  FORCE DENSITY 
    // ===============
    LOG_DBG(mechpde) << "Reading force densities";
    
    ReadRhsExcitation( "forceDensity", dispDofNames, ResultInfo::VECTOR, isComplex_, 
                        ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Force density must be defined on elements")
      }
      
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(),
                                            Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(),
                                           1.0, coef[i], coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1,3,3>(),
                                          Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(),
                                           1.0, coef[i], coefUpdateGeo);
        }
      }
      lin->SetName("ForceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
    
    
    // ==================
    //  THERMAL STRAIN
    // ==================
    LOG_DBG(mechpde) << "Reading thermal strain definition";

    ReadRhsExcitation( "thermalStrain", dispDofNames, ResultInfo::SCALAR, isComplex_,
                       ent, coef, coefUpdateGeo );

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("Thermal strain must be defined on elements")
        }

        RegionIdType myRegionId = ent[i]->GetRegion();
        MaterialType typeStrain, typeStress;
        if( subType_ == "axi" ) {
          typeStrain =  MECH_TEC_VECTORAXI;
          typeStress =  MECH_STIFFTENSOR_TEC_VECTORAXI;
        }
        else if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
          typeStrain =  MECH_TEC_VECTORPLANE;
          typeStress =  MECH_STIFFTENSOR_TEC_VECTORPLANE;
        }
        else if( subType_ == "3d") {
          typeStrain = MECH_TEC_VECTOR;
          typeStress = MECH_STIFFTENSOR_TEC_VECTOR;
        }
        else
          EXCEPTION("Mechanical Tensortype not implemented!")

        //get thermal expansion in Voigt notation
        PtrCoefFct vecTEC = materials_[myRegionId]->GetVectorCoefFnc(typeStrain, part);

        //get thermal expansion already multiplied with elasticity tensor
        PtrCoefFct vecMechTEC = materials_[myRegionId]->GetVectorCoefFnc(typeStress, part);

        PtrCoefFct refTemp = materials_[myRegionId]->GetScalCoefFnc(MECH_TEC_REFTEMPERATURE,
                                                                    part);

        PtrCoefFct TminusTref =
                  CoefFunction::Generate( mp_, part,
                                         CoefXprBinOp(mp_,coef[i],refTemp,CoefXpr::OP_SUB));

        //compute the termal strain
        PtrCoefFct thermalStrainCoef =
                    CoefFunction::Generate( mp_, part,
                                           CoefXprBinOp(mp_,TminusTref,vecTEC,CoefXpr::OP_MULT));
        // store the coefFunction for postprocessing
        thermalStrain_->AddRegion( myRegionId, thermalStrainCoef);


        PtrCoefFct thermalStressCoef =
                  CoefFunction::Generate( mp_, part,
                                         CoefXprBinOp(mp_,TminusTref,vecMechTEC,CoefXpr::OP_MULT));

        // store the coefFunction for postprocessing
        thermalStress_->AddRegion( myRegionId, thermalStressCoef);

        BaseBOperator * bOp = GetStrainOperator( isComplex_, false );

        if(isComplex_) {
            EXCEPTION("Complex thermal strain RHS linear form not implemented");
        }
        else {
            lin = new BUIntegrator<Double>(bOp, 1.0, thermalStressCoef, coefUpdateGeo);
        }
        lin->SetName("ThermalStrainInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        myFct->AddEntityList(ent[i]);
      } // for

    // ===============
    //  PRESSURE 
    // ===============
    LOG_DBG(mechpde) << "Reading mechanical pressure";
    StdVector<std::string> empty;
    ReadRhsExcitation( "pressure", empty, ResultInfo::VECTOR, isComplex_, 
                        ent, coef, coefUpdateGeo );
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Mechanical pressure must be defined on elements")
      }
            
      // Factor for the pressure:
      // The pressure is by definition in the opposite direction as the 
      // normal stress, i.e. a positive pressure means a compressive stress
      // (<0). Thus we have to take the minus sign into account
      const Double presFac = -1.0;
      
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,2>(),
                                                  Complex(presFac), coef[i],
                                                  volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true> ( new IdentityOperatorNormalTrans<FeH1,2>(),
                                                presFac, coef[i], volRegions,
                                                coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
                                                 Complex(presFac), coef[i],
                                                 volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
                                                 presFac, coef[i],
                                                 volRegions, coefUpdateGeo);
        }
      }
      lin->SetName("PressureInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for


    
    // ==================
    //  SURFACE TRACTION  
    // ==================
    LOG_DBG(mechpde) << "Reading surface tractions";
      
      ReadRhsExcitation( "traction", dispDofNames, ResultInfo::VECTOR, isComplex_, 
                          ent, coef, coefUpdateGeo );
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("Surface traction must be defined on elements")
        }
        // ensure that list contains only surface elements
        EntityIterator it = ent[i]->GetIterator();
        UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
        if( elemDim != (dim_-1) ) {
          EXCEPTION("Surface traction can only be defined on surface elements");
        }
        
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(),
                                              Complex(1.0), coef[i], coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(),
                                             1.0, coef[i],coefUpdateGeo);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3,3>(),
                                              Complex(1.0), coef[i], coefUpdateGeo);
          } else {
            lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(),
                                             1.0, coef[i], coefUpdateGeo);
          }
        }
        lin->SetName("TractionIntegrator");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        myFct->AddEntityList(ent[i]);
      } // for

      // ==================
      //  RHS NODAL VALUES
      // ==================
      LOG_DBG(mechpde) << "Reading direct right hand side values";

      ReadRhsExcitation( "rhsValues", dispDofNames, ResultInfo::VECTOR, isComplex_,
                          ent, coef, coefUpdateGeo );

      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        //for non-linear simulations we might need a conservative interpolation in each timestep...
        coef[i]->SetConservative(true);
        this->rhsFeFunctions_[MECH_DISPLACEMENT]->AddLoadCoefFunction(coef[i], ent[i]);
      }
  }

  BaseBDBInt*
  MechPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
                               RegionIdType regionId,
                               bool isComplex )
  {

    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );

    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex ) {
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR,
                                          tensorType_, Global::COMPLEX );
    } else {
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR,
                                          tensorType_, Global::REAL );
    }
    
    // TODO: do only for optimizaton
    if(true)
    {
      assert(!isComplex);

      CoefFunctionOpt<double>* tmpFnc = new CoefFunctionOpt<double>(curCoef);
      curCoef.reset(tmpFnc);
    }


    // store coefficient function for later use (e.g. in boundary integrators)
    regionStiffness_[regionId] = curCoef;
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    
    BaseBDBInt* integ = NULL;
    BaseBOperator* bOp = GetStrainOperator(isComplex, false);
    
    if( regionSoftening_[regionId] == "icModesTW") {
      // ===================
      //  ICModes Softening 
      // ===================
      BaseBOperator * gOp = GetStrainOperator( isComplex, true );
      if (isComplex ) 
        integ = new ICModesInt<Complex>(bOp, gOp, curCoef, 1.0);
      else
        integ = new ICModesInt<Double>(bOp, gOp, curCoef, 1.0);
      
    } else {
      // ====================
      //  Standard Stiffness 
      // ====================
      if (isComplex ) 
        integ = new BDBInt<Complex>(bOp, curCoef, 1.0);
      else
        integ = new BDBInt<Double>(bOp, curCoef, 1.0);
    }

    return integ;
  }

  

  BaseBOperator* MechPDE::GetStrainOperator(bool isComplex, bool icModes)
  {
    BaseBOperator* bOp = NULL;
    // determine if we do bloch eigenfrequency analysis
    bool do_bloch = domain->GetDriver()->DoBlochModeEigenfrequency();
    assert(!(do_bloch && !isComplex));


    if(isComplex)
    {
      if( subType_ == "planeStrain" )
      {
        if(do_bloch)
        {
          bOp = new StrainOperatorBloch2D<FeH1,Complex>(icModes);
          EigenFrequencyDriver* efd = dynamic_cast<EigenFrequencyDriver*>(domain->GetSingleDriver());
          dynamic_cast<StrainOperatorBloch2D<FeH1,Complex>* >(bOp)->SetWaveVector(efd->GetCurrentWaveVector());
        }
        else
          bOp = new StrainOperator2D<FeH1,Complex>(icModes);
      }
      if( subType_ == "axi" )
        bOp = new StrainOperatorAxi<FeH1,Complex>(icModes);
      if(subType_ == "planeStress")
        bOp = new StrainOperator2D<FeH1,Complex>(icModes);
      if(subType_ == "3d")
      {
        if(do_bloch)
        {
          bOp = new StrainOperatorBloch3D<FeH1,Complex>(icModes);
          EigenFrequencyDriver* efd = dynamic_cast<EigenFrequencyDriver*>(domain->GetSingleDriver());
          dynamic_cast<StrainOperatorBloch3D<FeH1,Complex>* >(bOp)->SetWaveVector(efd->GetCurrentWaveVector());
        }
        else
          bOp = new StrainOperator3D<FeH1,Complex>(icModes);
      }
    }
    else // not complex
    {
      if(subType_ == "axi")
        bOp = new StrainOperatorAxi<FeH1,Double>(icModes);
      if(subType_ == "planeStrain")
        bOp = new StrainOperator2D<FeH1,Double>(icModes);
      if(subType_ == "planeStress")
        bOp = new StrainOperator2D<FeH1,Double>(icModes);
      if(subType_ == "3d")
        bOp = new StrainOperator3D<FeH1,Double>(icModes);
    }

    if(bOp == NULL)
      EXCEPTION("strain operator not implemented for analysis type");

    return bOp;
  }

  void MechPDE::DefineSolveStep()
  {
		  solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MechPDE::InitTimeStepping()
  {
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(GLMScheme::NEWMARK, 0) );
    feFunctions_[MECH_DISPLACEMENT]->SetTimeScheme(myScheme);
  }

  void MechPDE::DefinePrimaryResults()
  {
    // Check for subType
    StdVector<std::string> dispDofNames, stressDofNames;

    if( subType_ == "3d" ) {
      dispDofNames = "x", "y", "z";
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";

    } else if( subType_ == "planeStrain" ) {
      dispDofNames = "x", "y";
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "planeStress" ) {
      dispDofNames = "x", "y";
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "axi" ) {
      dispDofNames = "r", "z";
      stressDofNames = "rr", "zz", "rz", "phiphi";
    }

    // === MECHANIC DISPLACEMENT ===
    shared_ptr<ResultInfo> disp(new ResultInfo);
    disp->resultType = MECH_DISPLACEMENT;
    disp->dofNames = dispDofNames;
    disp->unit = "m";
    disp->entryType = ResultInfo::VECTOR;
    disp->SetFeFunction(feFunctions_[MECH_DISPLACEMENT]);
    disp->definedOn = ResultInfo::NODE;
    feFunctions_[MECH_DISPLACEMENT]->SetResultInfo(disp);
    
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MECH_DISPLACEMENT] = "fix";
    idbcSolNameMap_[MECH_DISPLACEMENT] = "displacement";
    
    // this defines the primary unknown
    results_.Push_back( disp );
    
    // define functor for interpolation of the field
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
    DefineFieldResult( feFct, disp );
  }
    
  void MechPDE::DefinePostProcResults() {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<std::string> stressComponents;
    if( subType_ == "3d" ) {
      stressComponents = "xx", "yy", "zz", "yz", "xz", "xy";
    } else if( subType_ == "planeStrain" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "planeStress" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "axi" ) {
      stressComponents = "rr", "zz", "rz", "phiphi";
    }
    StdVector<std::string > dispDofNames;
    dispDofNames = feFunctions_[MECH_DISPLACEMENT]->GetResultInfo()->dofNames;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<BaseFeFunction> vFct;
    if ( analysistype_ != STATIC ) {
      // === MECHANIC VELOCITY ===
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = MECH_VELOCITY;
      vel->dofNames = dispDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::NODE;
      availResults_.insert( vel );
      DefineTimeDerivResult( MECH_VELOCITY, 1, MECH_DISPLACEMENT );
      vFct = timeDerivFeFunctions_[MECH_VELOCITY];

      // === MECHANIC ACCELERATION ===
      shared_ptr<ResultInfo> acc(new ResultInfo);
      acc->resultType = MECH_ACCELERATION;
      acc->dofNames = dispDofNames;
      acc->unit = "m/s^2";
      acc->entryType = ResultInfo::VECTOR;
      acc->definedOn = ResultInfo::NODE;
      availResults_.insert( acc );
      DefineTimeDerivResult( MECH_ACCELERATION, 2, MECH_DISPLACEMENT );
    }

    // === MECHANIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MECH_RHS_LOAD;
    rhs->dofNames = dispDofNames;
    rhs->unit = "N";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::NODE;
    rhsFeFunctions_[MECH_DISPLACEMENT]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MECH_DISPLACEMENT], rhs );

    // === MECHANIC STRESS ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressComponents;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    PtrCoefFct stressCoef;

    shared_ptr<CoefFunctionFormBased> sigmaFunc;
    if( isComplex_ ) {
      sigmaFunc.reset(new CoefFunctionFlux<Complex>(feFct, stress));
    } else {
      sigmaFunc.reset(new CoefFunctionFlux<Double>(feFct, stress));
    }
    // check for thermal strains
    bool isThermalStrain = false;
    PtrParamNode bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if( bcsNode ) {
      if( bcsNode->Has("thermalStrain") ) {
        isThermalStrain = true;
      }
    }
    if ( isThermalStrain )  {
      //! Cauchy stress is [c] ( Bu - alpha DeltaT )
      stressCoef =
          CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,sigmaFunc,thermalStress_,CoefXpr::OP_SUB));
    }
    else {
      stressCoef = sigmaFunc;
    }
    DefineFieldResult( stressCoef, stress );
    stiffFormCoefs_.insert(sigmaFunc);

    // === THERMOMECHANICAL STRESS ===
    if ( isThermalStrain )  {
      shared_ptr<ResultInfo> stress(new ResultInfo);
      stress->resultType = MECH_THERMAL_STRESS;
      stress->dofNames = stressComponents;
      stress->unit =  "N/m^2";
      stress->entryType = ResultInfo::TENSOR;
      stress->definedOn = ResultInfo::ELEMENT;
      DefineFieldResult( thermalStress_, stress );
    }

    // === VON_MISES_STRESS  ===
    shared_ptr<ResultInfo> misesStress(new ResultInfo);
    misesStress->resultType = VON_MISES_STRESS;
    misesStress->dofNames = "";
    misesStress->unit = "N/m^2";
    misesStress->entryType = ResultInfo::SCALAR;
    misesStress->definedOn = ResultInfo::ELEMENT;
    if ( !isComplex_ ) {
      shared_ptr<CoefFunctionCompound<Double> >
        misesCoef(new CoefFunctionCompound<Double>(mp_));
      std::map<std::string, PtrCoefFct> var;
      std::string misesStr;

      var["s"] = stressCoef;
      if ( stressDim_==3 ) {
        // 2D plane strain or stress case
        misesStr = "sqrt(   s_0_R^2 + s_1_R^2 \
                          - s_0_R*s_1_R \
                          + 3.0*s_2_R^2 )";
      }
      else if ( stressDim_==4 ) {
        // 2D axi case
        misesStr = "sqrt(   s_0_R^2 + s_1_R^2 + s_3_R^2\
                             - s_0_R*s_1_R - s_0_R*s_3_R - s_1_R*s_3_R  \
                             + 3.0*s_2_R^2 )";
      }
      else if ( stressDim_==6 ) {
        // 3D case
        misesStr = "sqrt(   s_0_R^2 + s_1_R^2 + s_2_R^2\
                           - s_0_R*s_1_R - s_0_R*s_2_R - s_1_R*s_2_R  \
                           + 3.0*(s_3_R^2 + s_4_R^2 + s_5_R^2) )";
      }
      else
        EXCEPTION( "Wrong dimesnion for stress: in DefinePostprocResults");

      misesCoef->SetScalar(misesStr,var);
      DefineFieldResult( misesCoef, misesStress );
    }

    // === MECHANIC STRAIN ===
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = MECH_STRAIN;
    strain->dofNames = stressComponents;
    strain->unit =  "";
    strain->entryType = ResultInfo::TENSOR;
    strain->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionFormBased> strainFunc;
    if( isComplex_ ) {
      strainFunc.reset(new CoefFunctionBOp<Complex>(feFct, strain));
    } else {
      strainFunc.reset(new CoefFunctionBOp<Double>(feFct, strain));
    }

    if ( isThermalStrain )  {
      //add thermal strain to mechanical strain
      shared_ptr<CoefFunction> totalStrain;
      totalStrain =
          CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,strainFunc,thermalStrain_,CoefXpr::OP_ADD));
      DefineFieldResult( totalStrain, strain );
    }
    else {
      DefineFieldResult( strainFunc, strain );
    }
    stiffFormCoefs_.insert(strainFunc);

    if ( isThermalStrain )  {
      shared_ptr<ResultInfo> thermalStrain(new ResultInfo);
      thermalStrain->resultType = MECH_THERMAL_STRAIN;
      thermalStrain->dofNames = stressComponents;
      thermalStrain->unit =  "";
      thermalStrain->entryType = ResultInfo::TENSOR;
      thermalStrain->definedOn = ResultInfo::ELEMENT;
      DefineFieldResult( thermalStrain_, thermalStrain );
    }

    PtrCoefFct intensFct;
    shared_ptr<CoefFunctionFormBased> kedFunc;
    shared_ptr<ResultFunctor> keFunc;
    shared_ptr<CoefFunctionSurf> sNormStructIntens;
    shared_ptr<ResultFunctor> powerFunc;
    if ( analysistype_ != STATIC ) {
      // === MECHANIC STRUCTURAL INTENSTIY ===
      shared_ptr<ResultInfo> intens(new ResultInfo);
      intens->resultType = MECH_STRUCT_INTENSTIY;
      intens->dofNames = dispDofNames ;
      intens->unit =  "N/ms";
      intens->entryType = ResultInfo::VECTOR;
      intens->definedOn = ResultInfo::ELEMENT;

      // The mechanic structural intensity calculates as
      // I = -[stress] * conj(v)
      PtrCoefFct velFnc = this->GetCoefFct( MECH_VELOCITY );
      // define temporary function, without the -1 sign
      PtrCoefFct intensTmp = 
          CoefFunction::Generate(mp_, part,
                                 CoefXprBinOp( mp_,  sigmaFunc, velFnc,
                                               CoefXpr::OP_MULT_VOIGT_TENSOR_VEC_CONJ ) );
      intensFct = 
          CoefFunction::Generate(mp_, part,
                                 CoefXprBinOp(mp_,  "-1.0", intensTmp , CoefXpr::OP_MULT ));
      DefineFieldResult( intensFct, intens );


      // === MECHANIC KINETIC ENERGY DENSITY ===
      shared_ptr<ResultInfo> kinEnergyDens(new ResultInfo);
      kinEnergyDens->resultType = MECH_KIN_ENERGY_DENS;
      kinEnergyDens->dofNames = "";
      kinEnergyDens->unit = "Ws/m^3";
      kinEnergyDens->entryType = ResultInfo::SCALAR;
      kinEnergyDens->definedOn = ResultInfo::ELEMENT;
      if( isComplex_ ) {
        kedFunc.reset(new CoefFunctionBdBKernel<Complex>(vFct, 0.5));
      } else {
        kedFunc.reset(new CoefFunctionBdBKernel<Double>(vFct, 0.5));
      }
      DefineFieldResult( kedFunc, kinEnergyDens );
      massFormCoefs_.insert(kedFunc);

      // === MECHANIC KINETIC ENERGY ===
      shared_ptr<ResultInfo> kinEnergy(new ResultInfo);
      kinEnergy->resultType = MECH_KIN_ENERGY;
      kinEnergy->dofNames = "";
      kinEnergy->unit = "Ws";
      kinEnergy->entryType = ResultInfo::SCALAR;
      kinEnergy->definedOn = ResultInfo::REGION;
      availResults_.insert ( kinEnergy );
      if( isComplex_ ) {
        keFunc.reset(new EnergyResultFunctor<Complex>(vFct, kinEnergy, 0.5));
      } else {
        keFunc.reset(new EnergyResultFunctor<Double>(vFct, kinEnergy, 0.5));
      }
      resultFunctors_[MECH_KIN_ENERGY] = keFunc;
      massFormFunctors_.insert(keFunc);

      // === MECHANIC NORMAL STRUCTURAL INTENSTIY ===
      shared_ptr<ResultInfo> intensNormal(new ResultInfo);
      intensNormal->resultType = MECH_NORMAL_STRUCT_INTENSITY;
      intensNormal->dofNames = "" ;
      intensNormal->unit =  "N/ms";
      intensNormal->entryType = ResultInfo::SCALAR;
      intensNormal->definedOn = ResultInfo::SURF_ELEM;
      sNormStructIntens.reset(new CoefFunctionSurf(true, intensNormal));
      DefineFieldResult( sNormStructIntens, intensNormal );
      surfCoefFcts_[sNormStructIntens] = intensFct;
      
      // === MECHANIC POWER ===
      shared_ptr<ResultInfo> power(new ResultInfo);
      power->resultType = MECH_POWER;
      power->dofNames = "";
      power->unit = "W";
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
      availResults_.insert( power );
      // then, integrate values
      if( isComplex_ ) {
        powerFunc.reset(new ResultFunctorIntegrate<Complex>(sNormStructIntens, 
                                                            feFct, power ) );
      } else {
        powerFunc.reset(new ResultFunctorIntegrate<Double>(sNormStructIntens, 
                                                           feFct, power ) );
      }
      resultFunctors_[MECH_POWER] = powerFunc;
    }


    // === MECHANIC DEFORMATION ENERGY DENSITY ===
    shared_ptr<ResultInfo> defEnergyDens(new ResultInfo);
    defEnergyDens->resultType = MECH_DEFORM_ENERGY_DENS;
    defEnergyDens->dofNames = "";
    defEnergyDens->unit = "Ws/m^3";
    defEnergyDens->entryType = ResultInfo::SCALAR;
    defEnergyDens->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionFormBased> dedFunc;
    if( isComplex_ ) {
      dedFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 0.5));
    } else {
      dedFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 0.5));
    }
    DefineFieldResult( dedFunc, defEnergyDens );
    stiffFormCoefs_.insert(dedFunc);

    // === MECHANIC DEFORMATION ENERGY ===
    shared_ptr<ResultInfo> defEnergy(new ResultInfo);
    defEnergy->resultType = MECH_DEFORM_ENERGY;
    defEnergy->dofNames = "";
    defEnergy->unit = "Ws";
    defEnergy->entryType = ResultInfo::SCALAR;
    defEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert( defEnergy );
    shared_ptr<ResultFunctor> deFunc;
    if( isComplex_ ) {
      deFunc.reset(new EnergyResultFunctor<Complex>(feFct, defEnergy, 0.5));
    } else {
      deFunc.reset(new EnergyResultFunctor<Double>(feFct, defEnergy, 0.5));
    }
    resultFunctors_[MECH_DEFORM_ENERGY] = deFunc;
    stiffFormFunctors_.insert(deFunc);

    // === MECHANIC TOTAL ENERGY DENSITY ===
    shared_ptr<ResultInfo> totEnergyDens(new ResultInfo);
    totEnergyDens->resultType = MECH_TOTAL_ENERGY_DENS;
    totEnergyDens->dofNames = "";
    totEnergyDens->unit = "Ws";
    totEnergyDens->entryType = ResultInfo::SCALAR;
    totEnergyDens->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunction> tedFunc;
    // in static analysis, the total energy density equals the deformation one
    if (analysistype_ == STATIC ) {
      tedFunc = dedFunc;
    } else {
      tedFunc = CoefFunction::Generate(mp_, part, 
                                       CoefXprBinOp( mp_, dedFunc, kedFunc, CoefXpr::OP_ADD) );
    }
    DefineFieldResult(tedFunc, totEnergyDens ); 

    // === MECHANIC TOTALENERGY ===
    shared_ptr<ResultInfo> tEnergy(new ResultInfo);
    tEnergy->resultType = MECH_TOTAL_ENERGY;
    tEnergy->dofNames = "";
    tEnergy->unit = "Ws";
    tEnergy->entryType = ResultInfo::SCALAR;
    tEnergy->definedOn = ResultInfo::REGION;
    availResults_.insert( tEnergy );
    shared_ptr<ResultFunctor> teFunc;
    if( isComplex_ ) {
      teFunc.reset(new ResultFunctorIntegrate<Complex>(tedFunc, feFct, tEnergy ) );
    } else {
      teFunc.reset(new ResultFunctorIntegrate<Double>(tedFunc, feFct, tEnergy ) );
    }
    resultFunctors_[MECH_TOTAL_ENERGY] = teFunc;


    // === MECHANIC DISPLACED SURFACE VOLUME ===
    shared_ptr<ResultInfo> dispNormal, dispVol;
    shared_ptr<CoefFunctionSurf> dispFctNormal;

    //normal mechanical displacement
    dispNormal.reset(new ResultInfo);
    dispNormal->resultType = MECH_NORMAL_DISPLACEMENT;
    dispNormal->dofNames = "";
    dispNormal->unit = "m";
    dispNormal->entryType = ResultInfo::SCALAR;
    dispNormal->definedOn = ResultInfo::SURF_ELEM;

    dispFctNormal.reset(new CoefFunctionSurf(true, dispNormal));
    DefineFieldResult(dispFctNormal, dispNormal);
    surfCoefFcts_[dispFctNormal] = feFct;

    dispVol.reset(new ResultInfo);
    dispVol->resultType = MECH_DEF_SURF_VOLUME;
    dispVol->dofNames = "";
    dispVol->unit = "m^3";
    dispVol->entryType = ResultInfo::SCALAR;
    dispVol->definedOn = ResultInfo::SURF_REGION;
    // Integrate normal displacement
    shared_ptr<ResultFunctor> dispVolFct;
    if( isComplex_ ) {
      dispVolFct.reset(new ResultFunctorIntegrate<Complex>(dispFctNormal,
                                                          feFct, dispVol ) );
    } else {
      dispVolFct.reset(new ResultFunctorIntegrate<Double>(dispFctNormal,
                                                         feFct, dispVol) );
    }
    resultFunctors_[MECH_DEF_SURF_VOLUME] = dispVolFct;
    availResults_.insert(dispVol);


    // === MECHANIC_NORMAL_STRESS ===
    shared_ptr<ResultInfo> normalStressInfo;
    shared_ptr<CoefFunctionSurf> normalStressFct;
    normalStressInfo.reset(new ResultInfo);
    normalStressInfo->resultType = MECH_NORMAL_STRESS;
    normalStressInfo->dofNames = dispDofNames;
    normalStressInfo->unit = "Pa";
    normalStressInfo->entryType = ResultInfo::VECTOR;
    normalStressInfo->definedOn = ResultInfo::SURF_ELEM;
    
    normalStressFct.reset(new CoefFunctionSurf(true, normalStressInfo));
    DefineFieldResult(normalStressFct, normalStressInfo);
    surfCoefFcts_[normalStressFct] = sigmaFunc;
  }
  
  std::map<SolutionType, shared_ptr<FeSpace> >
   MechPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
     
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    
    if( formulation == "default" || formulation == "H1" ){
      PtrParamNode potSpaceNode = infoNode->Get("mechDisplacement");
      crSpaces[MECH_DISPLACEMENT] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[MECH_DISPLACEMENT]->Init(solStrat_);

    }else{
       EXCEPTION( "The formulation " << formulation 
                  << "of the mechanic PDE is not known!" );
     }
     return crSpaces;
   }

} // end namespace CoupledField
