#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "heatCondPDE.hh"

#include "Forms/laplaceInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linNeumannInt.hh"

#include "PDE/trapezoidal.hh"

#include "DataInOut/writeresults.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Domain/ansatzFct.hh"

#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"


namespace CoupledField {


  // ======================================================
  // SET SOLUTION INFORMATION
  // ======================================================
  HeatCondPDE::HeatCondPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, 
                           WriteResults *aptOut)
    :SinglePDE(aptgrid,aptOut,aptTimeFunc) {
    ENTER_FCN( "HeatCondPDE::HeatCondPDE" );

    pdename_          = "heatConduction";
    pdematerialclass_ = THERMIC;

    nonLin_    = false;

    // Create new resultDof object
    shared_ptr<ResultDof> res1(new ResultDof);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = HEAT_TEMPERATURE;
    res1->dofNames = "t";
    res1->definedOn = ResultDof::NODE;
    res1->fctType = fct;
    results_.Push_back( res1 );

    //check, if problem is axisymmetric
    if ( params->HasValue( "type", "axi", "geometry" ) )
      isaxi_ = true;

  }

  void HeatCondPDE::ReadSpecialBCs() {
    ENTER_FCN( "HeatCondPDE::ReadSpecialBCs" );

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";

    // First, delete all of the "normal" boundary conditions
    inBcs_.Clear();
    


    StdVector<std::string> inName, inDof, inValue, inPhase, inDynamics;
    StdVector<std::string> inType;
    StdVector<Double> inHtc, inTSolid, inTFluid;
    
    // Get names of node sets, values and filenames for inhomogenous
    // Neumann boundary conditions
    keyVec = pdename_, "bcsAndLoads", "neumannInhom", "name";
    params->GetList(keyVec, attrVec, valVec, inName);

    keyVec.Last() = "entityType";
    params->GetList(keyVec, attrVec, valVec, inType);

    keyVec.Last() = "value";
    params->GetList(keyVec, attrVec, valVec, inValue);

    keyVec.Last() = "dof";
    params->GetList(keyVec, attrVec, valVec, inDof);
    
    keyVec.Last() = "dynamics";
    params->GetList(keyVec, attrVec, valVec, inDynamics);

    keyVec.Last() = "phase";
    params->GetList(keyVec, attrVec, valVec, inPhase);

    keyVec.Last() = "heatTransferCoefficient";
    params->GetList(keyVec, attrVec, valVec, inHtc);

    keyVec.Last() = "tempSolid";
    params->GetList(keyVec, attrVec, valVec, inTSolid );

    keyVec.Last() = "tempFluid";
    params->GetList(keyVec, attrVec, valVec, inTFluid);

    // Create inhomogeneous heat Neumann boundary conditions
    for( UInt i = 0; i < inName.GetSize(); i++ ) {
      shared_ptr<InhomHeatNeumannBc> actBc ( new InhomHeatNeumannBc );
      shared_ptr<EntityList> actList =
        ptgrid_->GetEntityList( EntityList::StringToType(inType[i]), 
                                inName[i] );
      actBc->entities = actList;
      actBc->result = results_[0];
      actBc->eqnMap = eqnMap_;
      if( inDof.GetSize() == 0 ) {
        actBc->dof = 1;
      } else {
        actBc->dof = results_[0]->GetDofIndex( inDof[i] );
      }
      actBc->value = inValue[i];
      actBc->phase = inPhase[i];
      actBc->dynamics = inDynamics[i];
      actBc->htc = inHtc[i];
      actBc->tSolid = inTSolid[i];
      actBc->tFluid = inTFluid[i];

      // add definition
      heatInBcs_.Push_back( actBc );
    }
  }

  void HeatCondPDE::DefineIntegrators() {
    ENTER_FCN( "HeatCondPDE::DefineIntegerators" );

    Double density, heatCapacity, thermalConductivity;
    Double coeffmass, coeffstiff;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );

      BaseMaterial * actMat = materials_[subdoms_[actSD]];
      actMat->GetScalar(density,DENSITY,REAL);
      actMat->GetScalar(heatCapacity,HEAT_CAPACITY,REAL);
      actMat->GetScalar(thermalConductivity,HEAT_CONDUCTIVITY,REAL);

      // stiffness integrator
      coeffstiff = thermalConductivity;
      BaseForm * bilinearStiff = new LaplaceInt(coeffstiff,isaxi_, true );
      BiLinFormContext * stiffContext = 
        new BiLinFormContext(bilinearStiff, STIFFNESS );

      stiffContext->SetPtPdes(this, this);
      stiffContext->SetResults( results_[0], results_[0],
                                 actSDList, actSDList );

      // mass integrator
      coeffmass = density * heatCapacity;
      BaseForm * bilinearMass  = new MassInt(coeffmass, 1, isaxi_, true );
      BiLinFormContext * massContext = 
        new BiLinFormContext(bilinearMass, MASS );

      massContext->SetPtPdes(this, this);
      massContext->SetResults( results_[0], results_[0],
                           actSDList, actSDList );


      // Finally add the standard integrators
      assemble_->AddBiLinearForm( stiffContext );
      assemble_->AddBiLinearForm( massContext );
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
    }

    // Neumann boundary condition
    for( UInt iBc = 0; iBc < heatInBcs_.GetSize(); iBc++ ) {
      
      // get current Bc
      InhomHeatNeumannBc const & actBc = *heatInBcs_[iBc];

      // we assume the surface normal points out of the domain,
      // but we want to take heat flux into the domain positive
      Double amplitude = -1.0 * String2Double(actBc.value);
      if (actBc.htc  != 0) {
        amplitude = actBc.htc  * ( actBc.tSolid - actBc.tFluid );
        
        Info->PrintF( pdename_, "For inhomogeneous Neumann BC use \
 \n  heat transfer coefficient: %f \n  TempSolid:  %f \n  TempFluid: %f\n\n",
                      actBc.htc, actBc.tSolid , actBc.tFluid );
      }
      LinearSurfForm *neumannBC = new LinNeumannInt( amplitude, 
                                                     HEAT_CONDUCTIVITY,
                                                     isaxi_ );
      neumannBC->SetVoluInfo( materials_ );
      
      LinearFormContext * neumannContext =
        new LinearFormContext( neumannBC );
      neumannContext->SetPtPde( this );
      neumannContext->SetResult( actBc.result, actBc.entities );
      
      assemble_->AddLinearForm( neumannContext );
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *actBc.result, actBc.entities );
    }

  }

  void HeatCondPDE::DefineSolveStep() {
    ENTER_FCN( "HeatCondPDE::DefineSolveStep" );

    solveStep_ = new StdSolveStep(*this);

  }

  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void HeatCondPDE::InitTimeStepping() {
    ENTER_FCN( "HeatCondPDE::InitTimeStepping" );
    
    // Until now no effective mass formulation in the trapezoidal 
    //  integration scheme is implemented!
    TS_alg_ = new Trapezoidal( algsys_ );

  }

  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void HeatCondPDE::InitCoupling(PDECoupling * Coupling) {
    ENTER_FCN( "HeatCondPDE::InitCoupling" );
    
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    
    // Intialize the memory of the coupling values
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      if (ptCoupling_->GetOutputQuantity(i) == HEAT_TEMPERATURE) {
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = true;
      }
    }
  }
  

  void HeatCondPDE::CalcOutputCoupling() {
    ENTER_FCN( "HeatCondPDE::CalcOutputCoupling" );

   //  UInt dof;
//     SolutionType quantity;
//     StdVector<Elem*> * couplingElems = NULL;
//     StdVector<UInt> * couplingNodes = NULL;
//     CFSVector * temp_values = NULL;
  
//     // at first, check if this PDE is iterative coupled
//     if (isIterCoupled_ == false)
//       return;

    // loop over all output coupling quantities
//     for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      
//       quantity = ptCoupling_->GetOutputQuantity(i);
//       ptCoupling_->GetOutputValues(i, temp_values);

//       // hard coded cast, since coupling is only possible with
//       // real valued entries
//       Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

//       switch(ptCoupling_->GetOutputType(i)) {
        
//       case NODE:

//         ptCoupling_->GetOutputNodes(i, couplingNodes);
//         ptCoupling_->GetOutputElements(i, couplingElems);
//         dof = ptCoupling_->GetOutputDof(i);    
  
//         break;
//       case ELEM:
//         Error("No Element coupling output", __FILE__,__LINE__);
//       }
//     }

  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void HeatCondPDE::WriteResultsInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "HeatCondPDE::WriteResultsInFile" );

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;

    NodeStoreSol<Double> * solTransient;    
    solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    if (saveSol_){
      outFile_->WriteNodeSolutionTransient(*solTransient,actStep,actTime);
    }

  }


  void HeatCondPDE::WriteHistoryInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "HeatCondPDE::WriteHistoryInFile" );

    NodeStoreSol<Double> solIm_mesh;
    NodeStoreSol<Double> * solTransient;
    
      
    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;
    
    solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    
    if (saveSolHist_){
      outFile_->WriteNodeHistoryTransient(*solTransient, actStep,actTime);
    }

  }


  void HeatCondPDE::ReadStoreResults() {  
    ENTER_FCN( "HeatCondPDE::ReadStoreResults" );
    
    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string quantity;
    
    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  


    Enum2String(HEAT_TEMPERATURE, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = true;
      hasOutput_ = true;
    }

    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    Enum2String(HEAT_TEMPERATURE, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = true;
      hasOutput_ = true;
      Info->PrintF( pdename_, "Saving temperature for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
      }
    }
  }


} // end of namespace CoupledField
