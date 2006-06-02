#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "heatCondPDE.hh"

#include "Forms/forms_header.hh"
#include "Forms/linNeumannInt.hh"

#include "PDE/scalarnodeEQN.hh"
#include "PDE/trapezoidal.hh"

#include "DataInOut/writeresults.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

#include "Driver/stdSolveStep.hh"


namespace CoupledField {


  // ======================================================
  // SET SOLUTION INFORMATION
  // ======================================================
  HeatCondPDE::HeatCondPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, 
                           WriteResults *aptOut)
    :SinglePDE(aptgrid,aptOut,aptTimeFunc) {
    ENTER_FCN( "HeatCondPDE::HeatCondPDE" );

    dofspernode_ = 1;
    solTypes_ = HEAT_TEMPERATURE;
    solDofs_ = 1;
    pdename_          = "heatConduction";
    pdematerialclass_ = THERMIC;

    nonLin_    = FALSE;

    //check, if problem is axisymmetric
    if ( params->HasValue( "type", "axi", "geometry" ) )
      isaxi_ = TRUE;

  }

  void HeatCondPDE::ReadSpecialBCs() {
    ENTER_FCN( "HeatCondPDE::ReadSpecialBCs" );

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    attrVec = "", "", "";
    valVec  = "", "", "";
    keyVec = pdename_,"bcsAndLoads","neumannInhom","heatTransferCoefficient";
    params->GetList(keyVec, attrVec, valVec, htc_);

    keyVec = pdename_, "bcsAndLoads", "neumannInhom", "tempSolid";
    params->GetList(keyVec, attrVec, valVec, tSolid_);

    keyVec = pdename_, "bcsAndLoads", "neumannInhom", "tempFluid";
    params->GetList(keyVec, attrVec, valVec, tFluid_);

    UInt ctr = htc_.GetSize();
    if ( (tSolid_.GetSize() != ctr) || (tFluid_.GetSize() != ctr) ) {
      (*error) << "Specify heatTransferCoefficient, tempSolid and tempFluid!";
      Error( __FILE__, __LINE__ );
    }
  }

  void HeatCondPDE::DefineIntegrators() {
    ENTER_FCN( "HeatCondPDE::DefineIntegerators" );

    Double density, heatCapacity, thermalConductivity;
    Double coeffmass, coeffstiff;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      BaseMaterial * actMat = materials_[subdoms_[actSD]];
      actMat->GetScalar(density,DENSITY,REAL);
      actMat->GetScalar(heatCapacity,HEAT_CAPACITY,REAL);
      actMat->GetScalar(thermalConductivity,HEAT_CONDUCTIVITY,REAL);

      // stiffness integrator
      coeffstiff = thermalConductivity;
      BaseForm * bilinearStiff = new LaplaceInt(coeffstiff,isaxi_);
      IntegratorDescriptor * stiffIntDescr = 
        new IntegratorDescriptor(bilinearStiff, STIFFNESS);

      stiffIntDescr->SetPDEIds(this, this);

      // mass integrator
      coeffmass = density * heatCapacity;
      BaseForm * bilinearMass  = new MassInt(coeffmass, dofspernode_, isaxi_);
      IntegratorDescriptor * massIntDescr = 
        new IntegratorDescriptor(bilinearMass, MASS);

      massIntDescr->SetPDEIds(this, this);


      // Finally add the standard integrators
      assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);
      assemble_->AddIntegrator(massIntDescr, subdoms_[actSD]);


    }

    // Neumann boundary condition
    StdVector<RegionIdType> IdVec;
    ptgrid_->RegionNameToId( IdVec, bcs_ni_ );
    for (UInt Id = 0; Id < bcs_ni_.GetSize(); Id++) {

      // we assume the surface normal points out of the domain,
      //  but we want to take heat flux into the domain positive
      Double amplitude = -1.0 * val_ni_[Id];
      if (htc_[Id] != 0) {
        amplitude = htc_[Id] * ( tSolid_[Id] - tFluid_[Id] );
        
        Info->PrintF( pdename_, "For inhomogeneous Neumann BC use \
\n  heat transfer coefficient: %f \n  TempSolid:  %f \n  TempFluid: %f\n\n",
                      htc_[Id], tSolid_[Id], tFluid_[Id] );
      }
      LinearSurfForm *neumannBC = new LinNeumannInt( amplitude, 
                                                     HEAT_CONDUCTIVITY,
                                                     isaxi_ );
      neumannBC->SetVoluInfo( materials_ );

      assemble_->AddRhsSrcSurfIntegrator( neumannBC,
                                          IdVec[Id],
                                          fncnames_ni_[Id],
                                          nonLin_ );
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
    
    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();

    // Until now no effective mass formulation in the trapezoidal 
    //  integration scheme is implemented!
    TS_alg_ = new Trapezoidal( algsys_, rhsSize );

  }

  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void HeatCondPDE::InitCoupling(PDECoupling * Coupling) {
    ENTER_FCN( "HeatCondPDE::InitCoupling" );
    
    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;
    
    // Intialize the memory of the coupling values
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
      if (ptCoupling_->GetOutputQuantity(i) == HEAT_TEMPERATURE) {
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = TRUE;
      }
    }
  }
  

  void HeatCondPDE::CalcOutputCoupling() {
    ENTER_FCN( "HeatCondPDE::CalcOutputCoupling" );

    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * temp_values = NULL;
  
    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
      return;

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
    NodeStoreSol<Complex> * solHarmonic;
    
      
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
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
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
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, "Saving temperature for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
      }
    }
  }


} // end of namespace CoupledField
