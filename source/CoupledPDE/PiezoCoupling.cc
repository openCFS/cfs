// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PiezoCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"

// include fespaces
#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/StrainOperator.hh"


namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode  )
    : BasePairCoupling( pde1, pde2, paramNode )
  {
    couplingName_ = "piezoDirect";
    materialClass_ = PIEZO;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();

  }


  // **************
  //   Destructor
  // **************
  PiezoCoupling::~PiezoCoupling() {
  }




  // ***************
  //   CalcResults
  // ***************
  void PiezoCoupling::CalcResults( shared_ptr<BaseResult> result ) {

    REFACTOR
//    // neighborregion for electric charge
//    RegionIdType neighbor = surfNeighborRegions_[result];
//
//    switch (result->GetResultInfo()->resultType ) {
//    case ELEC_CHARGE:
//      if ( isComplex_ ) {
//      	CalcCharges<Complex>( result, neighbor );
//      } else{
//        CalcCharges<Double>( result, neighbor );
//      }
//      break;
//
//    case MECH_STRESS:
//    case MECH_STRAIN:
//      if ( isComplex_ ) {
//        CalcStressStrain<Complex>( result );
//      } else {
//        CalcStressStrain<Double>( result );
//      }
//      break;
//
//    case VON_MISES_STRESS:
//      if(isComplex_)
//        CalcStressStrain<Complex>(result);
//      else
//        CalcStressStrain<Double>(result);
//      break;
//
//
//    case MECH_STRAIN_IRR:
//      if ( isComplex_ ) {
//        EXCEPTION("Ireversible Strain makes no sense in Harmonic analysis");
//      } else {
//        CalcStrainIrr( result );
//      }
//      break;
//
//    case ELEC_POLARIZATION:
//      if ( isComplex_ ) {
//        EXCEPTION("Electric Polarization makes no sense in Harmonic analysis");
//      } else {
//        CalcElecPolarization( result );
//        
//      }
//      break;
//
//    case ELEC_FLUX_DENSITY:
//      if ( isComplex_ ) {
//         CalcElecFluxDensity<Complex>( result );
//      } else {
//        CalcElecFluxDensity<Double>( result );
//      }
//      break;
//
//    default: 
//      WARN( "Resulttype not computable by piezoelectric coupling" );
//    }
  }


 

  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoCoupling::DefineIntegrators() {

    // get hold of both feFunctions
    shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);
    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> elecSpace = elecFct->GetFeSpace();
    
    
    // Global::ComplexPart matType = Global::REAL;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;
      // matType = Global::REAL;

      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      
      // ==================================================================
      //  STANDARD COUPLING INTEGRATOR
      // ==================================================================
      BiLinearForm * stiffInt = 
          GetStiffIntegrator( actSDMat, actRegion, complexMatData_[actRegion] );
      stiffInt->SetName("PiezoCouplingInt");
      BiLinFormContext * stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( dispFct, elecFct );
      stiffIntDescr->SetCounterPart(true);

      assemble_->AddBiLinearForm( stiffIntDescr );

    }

    
    //
//    Global::ComplexPart matType = Global::REAL;
//    RegionIdType actRegion;
//    // BaseMaterial * actSDMat = NULL; // TODO: Unused variable actSDMat
//
//
//    // get material from electrostatics
//    std::map<RegionIdType, BaseMaterial*> elecMat =
//      pde2_->getPDEMaterialData();
//
//
//    // Define integrators for "standard" materials
//    std::map<RegionIdType, BaseMaterial*>::iterator it;
//    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
//      // Set current region and material
//      actRegion = it->first;
//      // actSDMat = it->second;
//      matType = Global::REAL;
//
//      //transform the type
//      SubTensorType tensorType;
//      String2Enum(subType_,tensorType);
//
//      // create new entity list
//      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
//      actSDList->SetRegion( actRegion );
//
//      BiLinFormContext *actContextStiff = NULL;
//
//      if ( nonLinPiezoMicroHF_ ) {
//        // Piezoelectric problem using micro-modeling of Belov-Kreher
//
//        // get material from mechanics
//        std::map<RegionIdType, BaseMaterial*> mechMat = 
//          pde1_->getPDEMaterialData();
//
//        // get material from electrostatics
//        std::map<RegionIdType, BaseMaterial*> elecMat = 
//          pde2_->getPDEMaterialData();
//
//        // get time step size
//        Double dt;
//        dt = dynamic_cast<TransientDriver*>(domain->GetSingleDriver())->GetDeltaT();
//
//        //allocate for micro-modeling
//        StdVector<Elem*> elemssd;
//        ptGrid_->GetElems(elemssd, actRegion);
//        UInt numElSD =  elemssd.GetSize();
//        materials_[actRegion]->InitPiezoMicro(numElSD, actSDList,
//                                              mechMat[actRegion], 
//                                              elecMat[actRegion],
//                                              tensorType, dt);
//
//        // get solutions of PDEs
//        BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
//        BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();
//        
//        NodeStoreSol<Double> * solhelp1 = 
//          dynamic_cast<NodeStoreSol<Double>*>(solPDE1);    
//        NodeStoreSol<Double> * solhelp2 = 
//          dynamic_cast<NodeStoreSol<Double>*>(solPDE2);
//
//        // get previous solution of PDE
//        BaseNodeStoreSol * solPrevPDE1 = pde1_->getPDESolutionPrev();
//        BaseNodeStoreSol * solPrevPDE2 = pde2_->getPDESolutionPrev();
//        
//        NodeStoreSol<Double> * solPrevhelp1 = 
//          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE1);    
//        NodeStoreSol<Double> * solPrevhelp2 = 
//          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE2);
//
//        //
//        //------ add nonlinear coupling stiffness in mechanical equation ---
//        //
//        nLinMicroPiezoCouple* bilinearCoupleMech =  
//          new nLinMicroPiezoCouple( materials_[actRegion], mechMat[actRegion], 
//                                    elecMat[actRegion], tensorType, true);
//           
//        //set solution object for (n+1) and (n)
//        bilinearCoupleMech->SetSolutionMech(*solhelp1);
//        bilinearCoupleMech->SetSolutionElec(*solhelp2);
//        bilinearCoupleMech->SetPrevSolutionMech(*solPrevhelp1);
//        bilinearCoupleMech->SetPrevSolutionElec(*solPrevhelp2);
//
//        bilinearCoupleMech->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
//                                                results2_[0], this);
//
//        bilinearCoupleMech->SetMatDataType( matType );
//        
//        BiLinFormContext *actContextCouplemech =
//          new BiLinFormContext( bilinearCoupleMech, STIFFNESS  );
//        
//        actContextCouplemech->SetEntryType( matType );
//        actContextCouplemech->SetPtPdes( pde1_, pde2_ );
//        actContextCouplemech->SetResults( results1_[0], results2_[0],
//                                          actSDList, actSDList );
//      
//        assemble_->AddBiLinearForm( actContextCouplemech );
//
//        //
//        //------  add nonlinear coupling stiffness in electric equation ---
//        //
//        nLinMicroPiezoCouple* bilinearCoupleElec =  
//          new nLinMicroPiezoCouple( materials_[actRegion], mechMat[actRegion], 
//                                    elecMat[actRegion], tensorType, false);
//           
//        //set solution object for (n+1) and (n)
//        bilinearCoupleElec->SetSolutionMech(*solhelp1);
//        bilinearCoupleElec->SetSolutionElec(*solhelp2);
//        bilinearCoupleElec->SetPrevSolutionMech(*solPrevhelp1);
//        bilinearCoupleElec->SetPrevSolutionElec(*solPrevhelp2);
//
//        bilinearCoupleElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
//                                               results2_[0], this);
//
//        bilinearCoupleElec->SetMatDataType( matType );
//        
//        BiLinFormContext *actContextCoupleElec =
//          new BiLinFormContext( bilinearCoupleElec, STIFFNESS  );
//        
//        actContextCoupleElec->SetEntryType( matType );
//        actContextCoupleElec->SetPtPdes( pde2_, pde1_ );
//        actContextCoupleElec->SetResults( results2_[0], results1_[0],
//                                   actSDList, actSDList );
//      
//        assemble_->AddBiLinearForm( actContextCoupleElec );
//
//        //
//        //------  add nonlinear mechanical stiffness -------------------
//        //
//        nLinMicroPiezoMech* bilinearStiffMech =  
//          new nLinMicroPiezoMech( materials_[actRegion], mechMat[actRegion], 
//                                  elecMat[actRegion], tensorType);
//           
//        //set soltuion object for (n+1) and (n)
//        bilinearStiffMech->SetSolutionMech(*solhelp1);
//        bilinearStiffMech->SetSolutionElec(*solhelp2);
//        bilinearStiffMech->SetPrevSolutionMech(*solPrevhelp1);
//        bilinearStiffMech->SetPrevSolutionElec(*solPrevhelp2);
//
//        bilinearStiffMech->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
//                                              results2_[0], this);
//
//        bilinearStiffMech->SetMatDataType( matType );
//        
//        BiLinFormContext *actContextStiffMech =
//          new BiLinFormContext( bilinearStiffMech, STIFFNESS  );
//        
//        actContextStiffMech->SetEntryType( matType );
//        actContextStiffMech->SetPtPdes( pde1_, pde1_ );
//        actContextStiffMech->SetResults( results1_[0], results1_[0],
//                                         actSDList, actSDList );
//      
//        assemble_->AddBiLinearForm( actContextStiffMech );
//
//        //
//        //------  add nonlinear electrostatic stiffness ---
//        //
//        nLinMicroPiezoElec* bilinearStiffElec =  
//          new nLinMicroPiezoElec( materials_[actRegion], mechMat[actRegion], 
//                                  elecMat[actRegion], tensorType);
//           
//        //set soltuion object for (n+1) and (n)
//        bilinearStiffElec->SetSolutionMech(*solhelp1);
//        bilinearStiffElec->SetSolutionElec(*solhelp2);
//        bilinearStiffElec->SetPrevSolutionMech(*solPrevhelp1);
//        bilinearStiffElec->SetPrevSolutionElec(*solPrevhelp2);
//
//        bilinearStiffElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
//                                              results2_[0], this);
//
//        bilinearStiffElec->SetMatDataType( matType );
//        
//        BiLinFormContext *actContextStiffElec =
//          new BiLinFormContext( bilinearStiffElec, STIFFNESS  );
//        
//        actContextStiffElec->SetEntryType( matType );
//        actContextStiffElec->SetPtPdes( pde2_, pde2_ );
//        actContextStiffElec->SetResults( results2_[0], results2_[0],
//                                         actSDList, actSDList );
//      
//        assemble_->AddBiLinearForm( actContextStiffElec );
//
//        //
//        //------ micro-piezo-model: RHS for electric equation ---
//        //
//        MicroPiezoPolarizationElecRhsInt * elecRHS = 
//          new MicroPiezoPolarizationElecRhsInt( materials_[actRegion],
//                                                mechMat[actRegion],
//                                                elecMat[actRegion],
//                                                tensorType );	 
//        
//        elecRHS->SetSolutionMech(*solhelp1);
//        elecRHS->SetSolutionElec(*solhelp2);
//        elecRHS->SetPrevSolutionMech(*solPrevhelp1);
//        elecRHS->SetPrevSolutionElec(*solPrevhelp2);
//        
//        elecRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
//                                    results2_[0], this);
//        LinearFormContext * rhsContextElec = 
//          new LinearFormContext( elecRHS );
//        rhsContextElec->SetPtPde( pde2_ );
//        rhsContextElec->SetResult( results2_[0], actSDList );
//        assemble_->AddLinearForm( rhsContextElec );
//        
//        //
//        //------ micro-piezo-model: RHS for mechanic equation ---
//        //
//        MicroPiezoPolarizationMechRhsInt * mechRHS = 
//          new MicroPiezoPolarizationMechRhsInt( materials_[actRegion], 
//                                                mechMat[actRegion], 
//                                                elecMat[actRegion], 
//                                                tensorType);
//        mechRHS->SetSolutionMech(*solhelp1);
//        mechRHS->SetSolutionElec(*solhelp2);
//        mechRHS->SetPrevSolutionMech(*solPrevhelp1);
//        mechRHS->SetPrevSolutionElec(*solPrevhelp2);
//        
//        // for computation of electric field
//        mechRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,  
//                                    results2_[0], this);
//        LinearFormContext * rhsContextMech = 
//          new LinearFormContext( mechRHS );
//        rhsContextMech->SetPtPde( pde1_ );
//        rhsContextMech->SetResult( results1_[0], actSDList );
//        assemble_->AddLinearForm( rhsContextMech );
//      }
//
//      else if (  nonLinHysteresis_ ) {
//        //Piezoelectric problem with hysteresis
//
//        //allocate for hystersis modeling
//        StdVector<Elem*> elemssd;
//        ptGrid_->GetElems(elemssd, actRegion);
//        UInt numElSD =  elemssd.GetSize();
//        elecMat[actRegion]->InitHyst(numElSD, actSDList);
//
//        // get material from mechanics
//        std::map<RegionIdType, BaseMaterial*> mechMat =
//          pde1_->getPDEMaterialData();
//
//        // get material from electrostatics
//        std::map<RegionIdType, BaseMaterial*> elecMat =
//          pde2_->getPDEMaterialData();
//
//        // get solutions of PDEs
//        BaseNodeStoreSol * solPDE1 = pde1_->getPDESolution();
//        BaseNodeStoreSol * solPDE2 = pde2_->getPDESolution();
//
//        NodeStoreSol<Double> * solhelp1 =
//          dynamic_cast<NodeStoreSol<Double>*>(solPDE1);
//        NodeStoreSol<Double> * solhelp2 =
//          dynamic_cast<NodeStoreSol<Double>*>(solPDE2);
//
//        // get previous solution of PDE
//        BaseNodeStoreSol * solPrevPDE1 = pde1_->getPDESolutionPrev();
//        BaseNodeStoreSol * solPrevPDE2 = pde2_->getPDESolutionPrev();
//
//        NodeStoreSol<Double> * solPrevhelp1 =
//          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE1);
//        NodeStoreSol<Double> * solPrevhelp2 =
//          dynamic_cast<NodeStoreSol<Double>*>(solPrevPDE2);
//
//        // add nonlinear coupling stiffness in mechanical equation
//        nLinPiezoHystCouple* bilinearStiffNonLin =
//          new nLinPiezoHystCouple( materials_[actRegion], mechMat[actRegion],
//                                   elecMat[actRegion], tensorType, true);
//
//        //set solution object for (n+1) and (n)
//        bilinearStiffNonLin->SetSolutionMech(*solhelp1);
//        bilinearStiffNonLin->SetSolutionElec(*solhelp2);
//        bilinearStiffNonLin->SetPrevSolutionMech(*solPrevhelp1);
//        bilinearStiffNonLin->SetPrevSolutionElec(*solPrevhelp2);
//
//        bilinearStiffNonLin->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
//                                                results2_[0], this);
//
//        bilinearStiffNonLin->SetMatDataType( matType );
//
//        BiLinFormContext *actContextStiffCoupling =
//          new BiLinFormContext( bilinearStiffNonLin, STIFFNESS  );
//
//        // explicit definition of counter part!!
//        actContextStiffCoupling->SetCounterPart( false );
//
//        actContextStiffCoupling->SetEntryType( matType );
//        actContextStiffCoupling->SetPtPdes( pde1_, pde2_ );
//        actContextStiffCoupling->SetResults( results1_[0], results2_[0],
//                                   actSDList, actSDList );
//
//        assemble_->AddBiLinearForm( actContextStiffCoupling );
//
//        // add nonlinear coupling stiffness in electric equation
//        nLinPiezoHystCouple* bilinearCoupleElec =
//          new nLinPiezoHystCouple( materials_[actRegion], mechMat[actRegion],
//                                   elecMat[actRegion], tensorType, false);
//
//        //set solution object for (n+1) and (n)
//        bilinearCoupleElec->SetSolutionMech(*solhelp1);
//        bilinearCoupleElec->SetSolutionElec(*solhelp2);
//        bilinearCoupleElec->SetPrevSolutionMech(*solPrevhelp1);
//        bilinearCoupleElec->SetPrevSolutionElec(*solPrevhelp2);
//
//        bilinearCoupleElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
//                                               results2_[0], this);
//
//        bilinearCoupleElec->SetMatDataType( matType );
//
//        BiLinFormContext *actContextCoupleElec =
//          new BiLinFormContext( bilinearCoupleElec, STIFFNESS  );
//
//        // explicit definition og counter part!!
//        actContextCoupleElec->SetCounterPart( false );
//
//        actContextCoupleElec->SetEntryType( matType );
//        actContextCoupleElec->SetPtPdes( pde2_, pde1_ );
//        actContextCoupleElec->SetResults( results2_[0], results1_[0],
//                                   actSDList, actSDList );
//
//        assemble_->AddBiLinearForm( actContextCoupleElec );
//
//
//        // add nonlinear electrostattic stiffness
//        nLinPiezoHystElec* bilinearStiffElec =
//          new nLinPiezoHystElec( materials_[actRegion], mechMat[actRegion],
//                                 elecMat[actRegion], tensorType);
//
//        //set soltuion object for (n+1) and (n)
//        bilinearStiffElec->SetSolutionMech(*solhelp1);
//        bilinearStiffElec->SetSolutionElec(*solhelp2);
//        bilinearStiffElec->SetPrevSolutionMech(*solPrevhelp1);
//        bilinearStiffElec->SetPrevSolutionElec(*solPrevhelp2);
//
//        bilinearStiffElec->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
//                                              results2_[0], this);
//
//        bilinearStiffElec->SetMatDataType( matType );
//
//        BiLinFormContext *actContextStiffElec =
//          new BiLinFormContext( bilinearStiffElec, STIFFNESS  );
//
//        actContextStiffElec->SetEntryType( matType );
//        actContextStiffElec->SetPtPdes( pde2_, pde2_ );
//        actContextStiffElec->SetResults( results2_[0], results2_[0],
//                                         actSDList, actSDList );
//
//        assemble_->AddBiLinearForm( actContextStiffElec );
//
//
//        //hysteresis formulation: RHS for electric equation
//        PiezoPolarizationElecRhsInt * elecRHS =
//          new PiezoPolarizationElecRhsInt( materials_[actRegion],
//                                           mechMat[actRegion],
//                                           elecMat[actRegion],
//                                           tensorType );
//
//        elecRHS->SetSolutionMech(*solhelp1);
//        elecRHS->SetSolutionElec(*solhelp2);
//        elecRHS->SetPrevSolutionMech(*solPrevhelp1);
//        elecRHS->SetPrevSolutionElec(*solPrevhelp2);
//
//        elecRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
//                                    results2_[0], this);
//        LinearFormContext * rhsContextElec =
//          new LinearFormContext( elecRHS );
//        rhsContextElec->SetPtPde( pde2_ );
//        rhsContextElec->SetResult( results2_[0], actSDList );
//        assemble_->AddLinearForm( rhsContextElec );
//
//
//        // hysteresis formulation: RHS for mechanic equation
//        PiezoPolarizationMechRhsInt * mechRHS =
//          new PiezoPolarizationMechRhsInt( materials_[actRegion],
//                                           mechMat[actRegion],
//                                           elecMat[actRegion],
//                                           tensorType);
//        mechRHS->SetSolutionMech(*solhelp1);
//        mechRHS->SetSolutionElec(*solhelp2);
//        mechRHS->SetPrevSolutionMech(*solPrevhelp1);
//        mechRHS->SetPrevSolutionElec(*solPrevhelp2);
//
//        // for computation of electric field
//        mechRHS->Set4NonLinMaterial(ptGrid_, pde2_, eqnMap2_,
//                                    results2_[0], this);
//        LinearFormContext * rhsContextMech =
//          new LinearFormContext( mechRHS );
//        rhsContextMech->SetPtPde( pde1_ );
//        rhsContextMech->SetResult( results1_[0], actSDList );
//        assemble_->AddLinearForm( rhsContextMech );
//      }
//      else {
//        // add linear stiffness
//        BaseForm *bilinearStiff =
//          new linPiezoCoupling(materials_[actRegion], tensorType);
//        
//        bilinearStiff->SetMatDataType( matType );
//        
//        actContextStiff =
//          new BiLinFormContext( bilinearStiff, STIFFNESS  );
//      
//
//        // We also need to set the transposed of the coupling
//        // matrix to the lower diagonal side
//        actContextStiff->SetCounterPart( true );
//        
//        actContextStiff->SetEntryType( matType );
//        actContextStiff->SetPtPdes( pde1_, pde2_ );
//        actContextStiff->SetResults( results1_[0], results2_[0],
//                                     actSDList, actSDList );
//        
//        assemble_->AddBiLinearForm( actContextStiff );
//        
//        // check for complex valued material parameter
//        if( complexMatData_[actRegion] ) {
//          matType = Global::IMAG;
//          
//          BaseForm * bilinearStiffC =
//            new linPiezoCoupling( materials_[actRegion], tensorType);
//          
//          //GetStiffIntegrator(materialData_[actSD]);
//          BiLinFormContext *actComplexContextStiff =
//            new BiLinFormContext(bilinearStiffC, STIFFNESS );
//          actComplexContextStiff->SetPtPdes(pde1_, pde2_);
//          actComplexContextStiff->SetResults( results1_[0], results2_[0],
//                                              actSDList, actSDList );
//          // We also need to set the transposed of the coupling
//          // matrix to the lower diagonal side
//          actComplexContextStiff->SetCounterPart( true );
//          
//          actComplexContextStiff->SetEntryType(matType);
//          bilinearStiffC->SetMatDataType(matType);
//          assemble_->AddBiLinearForm( actComplexContextStiff );
//        }
//      }
//    }
//
//    // Define integrators for composite materials
//    // (only for subType "flatShell")
//    std::map<RegionIdType, Composite>::iterator compIt;
//    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
//         compIt++ ) {
//      // Get current subdomain and composite material
//      RegionIdType actRegion = compIt->first;
//      Composite * composite = &compIt->second;
//
//      // create new entity list
//      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
//      actSDList->SetRegion( actRegion );
//
//      FlatShellPiezoInt * compPiezoInt = new FlatShellPiezoInt( composite, false);
//      BiLinFormContext * stiffContext =
//        new BiLinFormContext( compPiezoInt, STIFFNESS);
//
//      // We also need to set the transposed of the coupling
//      // matrix to the lower diagonal side
//      stiffContext->SetCounterPart( true );
//
//      stiffContext->SetPtPdes( pde1_, pde2_ );
//      stiffContext->SetResults( results1_[0], results2_[0],
//                                actSDList, actSDList );
//      assemble_->AddBiLinearForm( stiffContext );
//
//      // check for complex valued material parameter
//      if( complexMatData_[actRegion] ) {
//        matType = Global::IMAG;
//        FlatShellPiezoInt * compPiezoIntC = new FlatShellPiezoInt( composite, false);
//        compPiezoIntC->SetMatDataType(matType);
//        BiLinFormContext * stiffContextC =
//            new BiLinFormContext( compPiezoIntC, STIFFNESS);
//
//        // We also need to set the transposed of the coupling
//        // matrix to the lower diagonal side
//        stiffContextC->SetCounterPart( true );
//
//        stiffContextC->SetPtPdes( pde1_, pde2_ );
//        stiffContextC->SetResults( results1_[0], results2_[0],
//                                  actSDList, actSDList );
//        stiffContextC->SetEntryType(matType);
//        assemble_->AddBiLinearForm( stiffContextC );
//
//      }
//    }

  }


  void PiezoCoupling::DefineAvailResults() {

    REFACTOR  
//   // Check for subType
//    StdVector<std::string> stressDofNames;
//    if( subType_ == "3d" ) {
//      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
//
//    } else if( subType_ == "planeStrain" ) {
//      stressDofNames = "xx", "yy", "xy";
//
//    } else if( subType_ == "planeStress" ) {
//      stressDofNames = "xx", "yy", "xy";
//
//    } else if( subType_ == "axi" ) {
//      stressDofNames = "rr", "zz", "rz", "phiphi";
//
//    } else if( subType_ == "flatShell" ) {
//      stressDofNames = "";
//    }
//
//    // Determine vectorial dofNames
//    StdVector<std::string> vecDofNames;
//    if( isaxi_) {
//      vecDofNames = "r", "z";
//    } else if (dim_ == 2) {
//      vecDofNames = "x", "y";
//    } else {
//      vecDofNames = "x", "y", "z";
//    }
//
//    // === MECHANIC STRESS ===
//    shared_ptr<ResultInfo> stress(new ResultInfo);
//    stress->resultType = MECH_STRESS;
//    stress->dofNames = stressDofNames;
//    stress->unit = "N/m^2";
//    stress->entryType = ResultInfo::TENSOR;
//    stress->definedOn = ResultInfo::ELEMENT;
//    stress->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( stress );
//
//    // === MECHANIC VON MISES STRESS (yield criterion, a scalar value)===
//    shared_ptr<ResultInfo> vonMises(new ResultInfo);
//    vonMises->resultType = VON_MISES_STRESS;
//    vonMises->dofNames = "";
//    vonMises->unit =  "";
//    vonMises->entryType = ResultInfo::SCALAR;
//    vonMises->definedOn = ResultInfo::ELEMENT;
//    vonMises->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( vonMises );
//
//    // === MECHANIC STRAIN ===
//    shared_ptr<ResultInfo> strain(new ResultInfo);
//    strain->resultType = MECH_STRAIN;
//    strain->dofNames = stressDofNames;
//    strain->unit = "";
//    strain->entryType = ResultInfo::TENSOR;
//    strain->definedOn = ResultInfo::ELEMENT;
//    strain->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( strain );
//
//    // === MECHANIC STRAIN Irreversibel===
//    shared_ptr<ResultInfo> strainIrr(new ResultInfo);
//    strainIrr->resultType = MECH_STRAIN_IRR;
//    strainIrr->dofNames = stressDofNames;
//    strainIrr->unit = "";
//    strainIrr->entryType = ResultInfo::TENSOR;
//    strainIrr->definedOn = ResultInfo::ELEMENT;
//    strainIrr->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( strainIrr );
//
//    // === ELECTRIC POLARIZATION ===
//    shared_ptr<ResultInfo> pol( new ResultInfo );
//    pol->resultType = ELEC_POLARIZATION;
//    pol->definedOn = ResultInfo::ELEMENT;
//    pol->entryType = ResultInfo::VECTOR;
//    pol->dofNames = vecDofNames;
//    pol->unit = "C/m^2";
//    availResults_.insert( pol );
//
//    // === ELECTRIC Flux Density ===
//    shared_ptr<ResultInfo> flux( new ResultInfo );
//    flux->resultType = ELEC_FLUX_DENSITY;
//    flux->definedOn = ResultInfo::ELEMENT;
//    flux->entryType = ResultInfo::VECTOR;
//    flux->dofNames = vecDofNames;
//    flux->unit = "C/m^2";
//    availResults_.insert( flux );
//
//    // === ELECTRIC CHARGE ===
//    shared_ptr<ResultInfo> charge(new ResultInfo);
//    charge->resultType = ELEC_CHARGE;
//    charge->dofNames = "";
//    charge->unit = "C";
//    charge->definedOn = ResultInfo::SURF_ELEM;
//    charge->entryType = ResultInfo::SCALAR;
//    charge->fctType = shared_ptr<ConstFct>(new ConstFct() );
//    availResults_.insert( charge );
  }

  BiLinearForm *
  PiezoCoupling::GetStiffIntegrator( BaseMaterial* actSDMat,
                                     RegionIdType regionId,
                                     bool isComplex ) {
    
    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );

    SubTensorType tensorType = NO_TENSOR;
    if( subType_ != "3d") {
      tensorType = PLANE_STRAIN;
    } else {
      tensorType = FULL;
    }
    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex ) {
      curCoef = actSDMat->GetCoefFunction(PIEZO_TENSOR,
                                          tensorType, Global::COMPLEX, true );
    } else {
      curCoef = actSDMat->GetCoefFunction(PIEZO_TENSOR,
                                          tensorType, Global::REAL, true );
    }
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    BiLinearForm * integ = NULL;
    if( subType_ == "axi" ) {
      integ = new ADBInt<StrainOperatorAxi<FeH1>,
                         GradientOperator<FeH1,2> >(curCoef, 1.0);
    } else if( subType_ == "planeStrain" ) {
      integ = new ADBInt<StrainOperator2D<FeH1>,
                         GradientOperator<FeH1,2> >(curCoef, 1.0);
    } else if( subType_ == "planeStress" ) {
      EXCEPTION("Not implemented");
    } else if( subType_ == "3d") {
      integ = new ADBInt<StrainOperator3D<FeH1>,
                         GradientOperator<FeH1,3> >(curCoef, 1.0);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }
    return integ;

  }
  
}

