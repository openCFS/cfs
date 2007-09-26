// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "bubblePDE.hh"

#include "Driver/solveStepODE.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/LinearKellerMiksis.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODEDescr/Gilmoredimlos.hh"
#include "ODESolve/ODESolver_RKF45.hh"
#include "ODESolve/ODESolver_Rosenbrock.hh"
#include "Forms/linearForm.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  BubblePDE::BubblePDE( Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) {


    // set flag for algsys
    needsAlgsys_ = false;
    
    // =====================================================================
    // Initialize all private variables
    // =====================================================================
    ptODESolver_ = NULL;
    bubbleDensity_ = 0.0;
    bubbleDynType_ = NOBUBBLETYPE;
    initRadius_ = 0.0;
    initVel_ = 0.0;
    density_ = 0.0;
    sonicVel_ = 0.0;
    pStatic_ = 0.0;
    pVapour_ = 0.0;
    surfaceTension_ = 0.0;
    polytrop_ = 0.0;
    viscosity_ = 0.0;
    t_ = 0.0;
    pressure_ = 0.0;
    pressureDeriv_ = 0.0;  
    pressureAmpl_ = 0.0;
    frequency_ = 0.0;
    dt_ = 0.0;
    dtStep_ = 0.0;
    steptime_ = 0.0;
    y_.Resize(2);


    // =====================================================================
    // set solution information
    // =====================================================================
    //  dofspernode_ = 1;
    //     solTypes_ = BUBBLE_RADIUS;
    //     solDofs_ = 1;
    pdename_          = "bubble";
    pdematerialclass_ = FLUID;
 
    nonLin_    = false;
    writeValues_ = false;
    writeRHS_    = false;

    
    // Create new resultDof object
    shared_ptr<ResultInfo> res1(new ResultInfo);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = BUBBLE_RADIUS;
    res1->dofNames = "";
    res1->definedOn = ResultInfo::ELEMENT;
    res1->fctType = fct;
    results_.Push_back( res1 );
    
    // **************************************************
    //   Check what type of bubble model should be used
    // **************************************************
    
    // vectors for parameter handling
    std::string bubbleDynString = "";
    
    bubbleDynType_ = NOBUBBLETYPE;
    
    ParamNode * bubbleNode = myParam_->Get( "bubbles" );
    
    if( bubbleNode ) {
      bubbleNode->Get( "bubbleType", bubbleDynString );
      String2Enum( bubbleDynString, bubbleDynType_ );
    }
    
    
    // Doppelt siehe unten ***
    //       //Set bubbledensity
    //       keyVec = pdename_, "bubbles", "bubbleNumDensity";
    //       params->Get(keyVec, attrVec, valVec, bubbleDensity_);
  }
  
  BubblePDE::~BubblePDE() {
    
    // close output stream
    if( !isIterCoupled_ ) {
      outValues_.close();

      UInt  numElems = ptgrid_->GetNumElems( subdoms_ );
      for( UInt iElem = 0; iElem < numElems; iElem++ ) {
        delete ptBubble_[iElem];
      }
    } else {
      delete ptBubble_[0];
    }


    delete ptODESolver_;

  }
  


  void BubblePDE::DefineIntegrators()
  {

    // Define entitlists for the eqnMap in order to get
    // a local<->global mapping of results
    for ( UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( subdoms_[actSD] );

      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], actSDList );
    }
    
    // =============================================
    //  Query ParamHandler for material parameters
    // =============================================
    
    // First make sure that there is a section 'bubbleDynamic'
    // We rely on params using a validating Schema parser here.
    ParamNode * bubbleDynNode = myParam_->Get( "bubbleDynamic", false );

    if ( !bubbleDynNode ) {
      EXCEPTION( "Cannot find initRadius! Assuming that section 'bubbleDynamic' "
                 "is missing in xml-file" );
    }

    bubbleDynNode->Get( "initRadius", initRadius_ );
    bubbleDynNode->Get( "initVel", initVel_ );
    bubbleDynNode->Get( "density", density_ );
    bubbleDynNode->Get( "sonicVel", sonicVel_ );
    bubbleDynNode->Get( "pStatic", pStatic_ );
    bubbleDynNode->Get( "pVapour", pVapour_ );
    bubbleDynNode->Get( "surfaceTension", surfaceTension_ );
    bubbleDynNode->Get( "polytrop", polytrop_);
    bubbleDynNode->Get( "viscosity", viscosity_ );

    //Set bubbledensity
    myParam_->Get( "bubbles" )->Get( "bubbleNumDensity", bubbleDensity_ );

    Info->PrintF ("", "The following paramters and methods were used\n");
    Info->PrintF ("", "to compute the bubble behaviour:\n\n");
    Info->PrintF ("", "Initial Radius of the bubbles:\t %10.6e\n", initRadius_);
    Info->PrintF ("", "Initial Velocity of the bubble wall:\t %10.6e\n", initVel_);
    Info->PrintF ("", "Density of bubbles per unit volume:\t %10.6e\n", bubbleDensity_);
    Info->PrintF ("", "Density of the fluid:\t %10.6e\n", density_);
    Info->PrintF ("", "Sound velocity of the fluid:\t %10.6e\n", sonicVel_);
    Info->PrintF ("", "Static pressure:\t %10.6e\n", pStatic_);
    Info->PrintF ("", "Vapour pressure of the fluid:\t %10.6e\n", pVapour_);
    Info->PrintF ("", "Surface tension of the fluid:\t %10.6e\n", surfaceTension_);
    Info->PrintF ("", "Polytropic exponent of the fluid:\t %10.6e\n", polytrop_);
    Info->PrintF ("", "Viscosity of the fluid:\t %10.6e\n\n", viscosity_);

    // Generate ODE solver object
    ptODESolver_ = new ODESolver_RKF45;
    Info->PrintF( pdename_, "Using the Runge-Kutta-Method to solve bubbledynamics\n");

    
    myParam_->Get("excitation")->Get( "frequency", frequency_ );
    
    // Hack for determining iterative coupling:
    // Check the number of PDEs defined in this step. If there is more than one,
    // we assume iterative coupling
    StdVector<ParamNode*> pdeNodes = 
      param->Get("sequenceStep", "indes", GenStr( sequenceStep_ ) )->Get("pdeList")->GetChildren();
    if( pdeNodes.GetSize() > 1 ) {
      isIterCoupled_ = true;
    } else {
      isIterCoupled_ = false;
    }


    // Check for iterative coupling
    UInt numElems = 0;
    if( isIterCoupled_ ) {
      numElems = ptgrid_->GetNumElems( subdoms_ );

      // Resize vectors
      radiusOldStep_.Resize( numElems );
      velocityOldStep_.Resize( numElems );

      // Generate initial data for each element in the iteration workingcopy
      for( UInt iElem = 0; iElem < numElems; iElem++ ) {
      
        // set initial data for each element
        radiusOldStep_[iElem]    = initRadius_;
        velocityOldStep_[iElem]  = initVel_;
      }
      

    } else {
      numElems = 1;

      // Read additional data for pressure and requency
      myParam_->Get( "excitation")->Get( "pressure", pressureAmpl_ );
      //params->Get( "frequency", frequency_,  pdename_, "excitation" );

      //+++++Dimensionless case ++++++++++++
      pressureAmpl_= pressureAmpl_ / pStatic_;
      frequency_   = frequency_ * initRadius_ /(sqrt(pStatic_/ density_));  
      //+++++Dimensionless case ++++++++++++

      // Open output file stream
      outValues_.open("out.dat" );
    } 

    // Resize vectors
    radius_.Resize( numElems );
    velocity_.Resize( numElems );
    hTry_.Resize( numElems );
    ptBubble_.Resize( numElems );

    
    // Generate ODE-Description and initial data for each element
    for( UInt iElem = 0; iElem < numElems; iElem++ ) {
      
      // set initial data for each element
      radius_[iElem]    = initRadius_;
      velocity_[iElem]  = initVel_;
      hTry_[iElem]      = 1e-10;
      
      // Generate ODE-Description for each element
      switch(bubbleDynType_){

      case KELLERMIKSIS:
        ptBubble_[iElem] = new KellerMiksis(initRadius_,density_, sonicVel_, pStatic_, 
					    pVapour_, surfaceTension_, polytrop_, viscosity_);
	if( iElem == 0 ) {
	Info->PrintF( pdename_, "Using Keller-Miksis-Bubble-Model\n");
	}

//         ptBubble_[iElem] = 
//           new LinearKellerMiksis(initRadius_,density_, sonicVel_, pStatic_, 
//                                  pVapour_, surfaceTension_, polytrop_, viscosity_);
//         if( iElem == 0 ) {
//           Info->PrintF( pdename_, "Using Linear-Keller-Miksis-Bubble-Model\n");
//         }
        
        break;
        
        
      case GILMORE:
        //ptBubble_[iElem] = new Gilmore(initRadius_,density_, sonicVel_, pStatic_, 
        //			pVapour_, surfaceTension_, polytrop_, viscosity_);
        // 	  Info->PrintF( pdename_, "Using Gilmore-Bubble-Model\n");
        ptBubble_[iElem] = 
          new Gilmoredimlos(initRadius_,density_, sonicVel_, pStatic_, 
                            pVapour_, surfaceTension_, polytrop_, viscosity_);
        if( iElem == 0 ) {
          Info->PrintF( pdename_, "Using dimensionless Gilmore-Bubble-Model\n");
        }
        break;
        
      default:
        EXCEPTION( "No bubblemethod specified " );
      }
    } // end for
    
  }




  void BubblePDE::DefineSolveStep() {
    solveStep_ = new SolveStepODE(*this);
  }


  // ======================================================
  // SOLVING SECTION
  // ======================================================
  
  void BubblePDE::Solve() {
   
    // Calculate basic data
    mHandle_ =  domain->GetMathParser()->GetNewHandle();
    MathParser * mParser =  domain->GetMathParser();
    mParser->SetExpr( mHandle_, "dt" );
    dt_ = mParser->Eval( mHandle_ );
    dtStep_ = dt_;

    steptime_ =  solveStep_->GetActTime() - dt_;
    
    
    dt_ = dt_ / initRadius_ * (sqrt( pStatic_ / density_));
    steptime_   = steptime_ / initRadius_ * (sqrt(pStatic_/ density_));
    tDim_ = solveStep_->GetActTime();
    t_ = tDim_ / initRadius_ * (sqrt(pStatic_/ density_));

    // Check if we are coupled iterative
    if( isIterCoupled_ == true ) {

      if (iterCoupledCounter_ == 0){
        radiusOldStep_   = radius_;
        velocityOldStep_ = velocity_;
      }

      // NOTE: We are iterating over the couplingElems_ entries, which might have a 
      // different ordering compared to the elements we obtain by calling ptGrid_->GetElems(subdoms_)!!!
      for( UInt iElem = 0; iElem < couplElems_->GetSize(); iElem++ ) {
        
        // get pressure value and derivative of current element
        pressure_ = (*pressBuf_)[iElem];
        pressureDeriv_ = (*pressDerivBuf_)[iElem];



        //Dimensionless case ++++++++++++++++++++++++++++++++++++++++++++++++

        //Diemensionless  
        pressure_ = pressure_ / pStatic_ ;
        pressureDeriv_ = pressureDeriv_ * initRadius_ / pStatic_ / ( sqrt ( pStatic_ / density_));
        
        //Druck�bergabe
        ptBubble_[iElem]->SetP(pressure_);
        ptBubble_[iElem]->SetDpdt(pressureDeriv_);
        
        if ( hTry_ [iElem]> dt_){
          hTry_[iElem] = dt_ / 3.0;
        }

        hTry_[iElem]      = hTry_[iElem] / initRadius_ * (sqrt(pStatic_/ density_));

        //get the current values
        y_[0] = radiusOldStep_[iElem] / initRadius_;
        y_[1] = velocityOldStep_[iElem] / (sqrt( pStatic_/ density_));

        //set iElem to ODE-Solver
        ptODESolver_->SetNumEl(iElem);
		  
        ptODESolver_->Solve(steptime_, t_, y_, *ptBubble_[iElem], hTry_[iElem],
                            0, dt_);
        //set the new values
        radius_[iElem]   = y_[0] * initRadius_;
        velocity_[iElem] = y_[1] * sqrt( pStatic_/ density_);
        hTry_[iElem]     = hTry_[iElem]  * initRadius_ / (sqrt(pStatic_/ density_));


        //Dimensionless case ++++++++++++++++++++++++++++++++++++++++++++++++++++

        iterCoupledCounter_++;
      }
      

    } else {
      
      //+++++Dimensionless case ++++++++++++

  


      hTry_[0]      = hTry_[0] / initRadius_ * (sqrt(pStatic_/ density_));  // Kann man evtl. woanders machen
      
      pressure_ = - pressureAmpl_ * sin( 2 * PI * frequency_ * (t_));
      pressureDeriv_  = - pressureAmpl_ * (2 * PI * frequency_ ) * 
        cos( 2 * PI * frequency_ * (t_));   
      ptBubble_[0]->SetP(pressure_);
      ptBubble_[0]->SetDpdt(pressureDeriv_);
      
      
      //get the current values
      y_[0] = radius_[0] / initRadius_;
      y_[1] = velocity_[0] / (sqrt( pStatic_/ density_));
      
      ptODESolver_->Solve(steptime_, t_, y_, *ptBubble_[0], hTry_[0],
                          0, dt_);
      
      //set the new values
      radius_[0]   = y_[0] * initRadius_;
      velocity_[0] = y_[1] * sqrt( pStatic_/ density_);
      hTry_[0]     = hTry_[0]  * initRadius_ / (sqrt(pStatic_/ density_)); // Kann man evtl. woanders machen
      //+++++Dimensionless case ++++++++++++ 
       
      //+++++Ausgabe von radius_, velocity_, berechnetem Druck, zeit ---- > 
      outValues_ << solveStep_->GetActTime() << '\t' 
                 << pressure_ *pStatic_ << '\t'
                 << radius_[0] << '\t'
                 << velocity_[0] << std::endl;
       
    }
    
  }


  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================


//   // **********************
//   //   WriteResultsInFile
//   // **********************
//   void BubblePDE::WriteResultsInFile( const UInt kstep,
//                                       const Double asteptime,
//                                       UInt stepOffset,
//                                       Double timeOffset ) {
    

//     UInt actStep = kstep + stepOffset;
//     Double actTime = timeOffset + asteptime;

    
//     if( writeRHS_ == true ) {
//       outFile_->WriteElemSolutionTransient( addElemResult_, actStep, actTime);
//     }

//     if( writeValues_ == true ) {
      
      
//       UInt numElems = couplElems_->GetSize();

//       for (UInt el = 0; el < numElems; el++) {
//         Vector<Double> result(3);
 
// 	//zum Testen
// // 	if ( actTime  < 1.0 / frequency_){
	 
// // 	  result[0] = 0.0;
// // 	  result[1] = 0.0;
// // 	  result[2] = 0.0; 
// // 	}
// // 	else{
// // 	  if (el = numElems - 1 ){
// // 	    Info->PrintF( pdename_, "Ausgabe der Bubbleergebnisse %e \n", actTime);
// // 	  }

// 	  result[0] = radius_[el];
// 	  result[1] = velocity_[el];
// 	  result[2] = (4.0/3.0) * PI * bubbleDensity_ * radius_[el]
// 	    * radius_[el] * radius_[el];

// 	//Achtung hier m�sste ich noch das Maximum �ber die Zeit nehmen
// 	//dann habe ich allerdings nur einen Wert, es k�nnte ja nur 
// 	//einmal erreicht werden, was eigentlich zum Reinigen  nicht
// 	//ausreicht
// // 	//Determination of collaps -> cleaning
// // 	Double temp;
// // 	temp = result[0] - initRadius_;
// // 	if ( temp <= 1.0 ){
// // 	  result[2] = 1.0;
// // 	}
// // 	else if( temp <= 1.5 ){
// // 	  result[2] = 2.0;
// // 	}
// // 	else if( temp <= 2.0 ){
// // 	  result[2] = 3.0;
// // 	}
// // 	else if( temp <= 5.0 ){
// // 	  result[2] = 4.0;
// // 	}
// // 	else{
// // 	  result[2] = 5.0;
// // 	}


//         //        if (el == 90)
//         //std::cerr<<actTime<<"   " <<el<< "   " << result[0] << "   " 
//         // result[1] << "     " << result[2] << std::endl;

// 	  //	}
        
//         // Map global element number to local one

//         UInt locElemNum = 
//           eqnMap_->Mesh2PdeElem((*couplElems_)[el]->elemNum );
//         elemResult_.SetElemResult( locElemNum-1,result );
//       }
      
//       outFile_->WriteElemSolutionTransient(elemResult_, actStep, actTime);
//     }

//   }


//   // **********************
//   //   WriteHistoryInFile
//   // **********************
//   void BubblePDE::WriteHistoryInFile( const UInt kstep,
//                                       const Double asteptime,
//                                       UInt stepOffset,
//                                       Double timeOffset ) {



//     if (analysistype_== TRANSIENT) {

//       UInt actStep = kstep + stepOffset;
//       Double actTime = timeOffset + asteptime;

//       if (saveElemBubbleHist_.GetSize() != 0) {


// 	if( writeValues_ == false ) {
// 	  UInt numElems = couplElems_->GetSize();

// 	  for (UInt el = 0; el < numElems; el++) {
// 	    Vector<Double> result(3);
            
// 	    result[0] = radius_[el];
// 	    result[1] = velocity_[el];
// 	    result[2] = (4.0/3.0) * PI * bubbleDensity_ * radius_[el]
// 	      * radius_[el] * radius_[el];

// 	    UInt locElemNum = 
// 	      eqnMap_->Mesh2PdeElem((*couplElems_)[el]->elemNum );
// 	    elemResult_.SetElemResult( locElemNum-1,result );
// 	  }
// 	}
	
// 	outFile_->WriteElemHistoryTransient(elemResult_, actStep, 
// 					    actTime);
//       }
//       if (saveElemRHSHist_.GetSize() != 0) {
	
// 	outFile_->WriteElemHistoryTransient(addElemResult_, actStep, 
// 					    actTime);
//       }



//     }

//   }



  // ***************
  //   PostProcess
  // ***************
  void BubblePDE::CalcResults( shared_ptr<BaseResult> result) {
    EXCEPTION( "Implement" );
 
  }
  


  // ======================================================
  // COUPLING SECTION
  // ======================================================
  void BubblePDE::CalcInputCoupling() {

    CFSVector * tempValues = NULL;
    SolutionType quantity;
    
    // Outer loop over all INPUT coupling terms
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {
      
      ptCoupling_->GetInputValues(i, tempValues);
      quantity = ptCoupling_->GetInputQuantity(i);      

      switch(ptCoupling_->GetInputType(i)) {

      case MAT:

        // get coupling elements
        ptCoupling_->GetInputElements( i, couplElems_ );

        if( quantity == ACOU_PRESSURE ) {
          pressBuf_ = dynamic_cast<Vector<Double>*>(tempValues);
        }
        if( quantity == ACOU_PRESSURE_DERIV_1 ) {
          pressDerivBuf_ = dynamic_cast<Vector<Double>*>(tempValues);
        }

        // Print values on screen
        // for( UInt iElem = 0; iElem < elems->GetSize(); iElem++ ) {
        //std::cerr << "ElemNr: " << (*elems)[iElem]->elemNum
        //          << ", value: " << values[iElem] << std::endl;
        //} 
        
        break;
      default :
        EXCEPTION( "This type of input coupling is not possible for BubblePDE!" );
        break;
        
      }
    }
  }


  void BubblePDE::InitCoupling(PDECoupling * Coupling)
  {
  
    isIterCoupled_ = true;
    ptCoupling_   = Coupling;
    
    StdVector<std::string> * nRegions;
    StdVector<RegionIdType> nRegionIds;
    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();

    std::cerr << "bubbleInitCoupl  numCoupling = " << numCouplings << std::endl;
    std::cerr << "bubbleInitCoupl numOutputCouplings = " << ptCoupling_->GetNumOutputCouplings() << std::endl;
    nonLin_ = false;

    // Initialization of coupling helper arrays
    std::string quantity;
    StdVector<UInt> * couplingnodes = NULL;

    for (UInt i = 0; i < numCouplings; i++) {
      std::cerr << "First time in BubblePDE::InitCoupling\n";
      if (ptCoupling_->GetOutputQuantity(i) == ACOU_BUBBLE_RHS_VAL) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector( i, isComplex_ );
         
      }

      if (ptCoupling_->GetOutputQuantity(i) == BUBBLE_RADIUS) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector( i, isComplex_ );
         
      }

      if (ptCoupling_->GetOutputQuantity(i) == BUBBLE_RADIUS_DERIV_1) {

        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector( i, isComplex_ );
         
      } 
  
       
    } 
  }
  


  void BubblePDE::CalcOutputCoupling()
  {

    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes     = NULL;
    CFSVector * values = 0;
    UInt forcesCount = 0;
    StdVector<std::string> regions;

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // loop over all output coupling quantities
    for (UInt actCoupling=0; 
         actCoupling<ptCoupling_->GetNumOutputCouplings(); 
         actCoupling++) {
      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      ptCoupling_->GetOutputValues(actCoupling, values);
      
      Vector<Double> * temp = dynamic_cast<Vector<Double> *>(values);
      
      switch(ptCoupling_->GetOutputType(actCoupling)) {
        
      case NODE:        
        ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
        ptCoupling_->GetOutputRegions(actCoupling, regions);
        
        // Trigger calculation of RHS values
        CalcAcouRHS( regions, *couplingNodes, *temp );
        
        break;
        
      case ELEM:
        ptCoupling_->GetOutputElements(actCoupling, couplingElems);
        if ( quantity == BUBBLE_RADIUS ||
             quantity == BUBBLE_RADIUS_DERIV_1 ) {
          CalcBubbleRadius( *couplingElems, *temp, quantity );
        }
	//   EXCEPTION("No Element coupling output" );
        break;
      }
    }
  }
  
  void BubblePDE::CalcAcouRHS( StdVector<std::string>& regions, 
                               StdVector<UInt>& nodes,
                               Vector<Double>& values ) {


    Vector<Double> bubbleValues(2);

    // Initialize values vector!
    values.Init();

    // Convert region names to Ids
    StdVector<RegionIdType> regionIds;
    ptgrid_->RegionNameToId( regionIds, regions );

    // create rhs integrator
    std::string dummy = "1.0";
    VolumeSrcInt *rhsForm = new VolumeSrcInt(dummy, isaxi_);        

    for( UInt iRegion = 0; iRegion < regionIds.GetSize(); iRegion++ ) {

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( regionIds[iRegion] );
      
      EntityIterator it = actSDList.GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++ ) {

        BaseFE * ptEl = it.GetElem()->ptElem;
        StdVector<UInt> const & connecth = it.GetElem()->connect;

        // Find correct index for this element in vector couplElems_
        Integer couplIndex = -1;

        for( UInt iCouplElem = 0; iCouplElem < couplElems_->GetSize(); iCouplElem++ ) {
          if( (*couplElems_)[iCouplElem]->elemNum == it.GetElem()->elemNum ) {
            couplIndex = iCouplElem;
            break;
          }
        }
        
        if( couplIndex == -1 ) {
          EXCEPTION( "Found no matching element for coupling element Nr" 
                     << it.GetElem()->elemNum << "!" );
        }
          

        Double beta2 = 0.0;
        Double beta2Com = 0.0;
        Double Rpp = 0.0; 

	//Hier Abfrage wegen Einschwingen
       
        bubbleValues[0] = radius_[couplIndex];
        bubbleValues[1] = velocity_[couplIndex];

	//Compute bubble acceleration
	StdVector<Double> dydt(2);
	y_[0] = bubbleValues[0] / initRadius_;
	y_[1] = bubbleValues[1] / (sqrt( pStatic_/ density_));
	t_    = solveStep_->GetActTime()  / initRadius_ * (sqrt(pStatic_/ density_));
	ptBubble_[couplIndex]->CompDeriv(t_, y_, dydt);
	dydt[1] = dydt[1] * pStatic_ / ( density_ * initRadius_ );
	Rpp=dydt[1];

	//Compute factor for RHS
	//ATTENTION only compatible with Keller-Miksis
	Double temp, temp2, temp3;
	temp = -3.0 / 2.0 * ( 1.0 - (bubbleValues[1] / 3.0 / sonicVel_) ) 
	  * bubbleValues[1] * bubbleValues[1];
	temp -= 2.0 * surfaceTension_ / density_ / bubbleValues[0];
	temp -= 4.0 * viscosity_ * bubbleValues[1] / density_ / bubbleValues[0];
	temp += 1.0 / density_ * ( pStatic_ - pVapour_ + 2.0 * surfaceTension_ / initRadius_)
	  * std::pow((initRadius_/bubbleValues[0]),(3.0*polytrop_)) 
	  * (1.0 + bubbleValues[1] / sonicVel_ * ( 1.0 - 3.0 * polytrop_));
	temp += (1.0 + (bubbleValues[1] / sonicVel_)) /density_ * (pVapour_ - pStatic_);
	
	temp2  = 4.0 * PI * density_ * bubbleDensity_ * bubbleValues[0] * bubbleValues[0] / ( ( 1.0 -( bubbleValues[1] / sonicVel_)) * bubbleValues[0]  + 4.0 * viscosity_ / density_ /sonicVel_ );

	temp3 = 8.0 * PI * density_ * bubbleDensity_ * bubbleValues[0] 
	  * bubbleValues[1] * bubbleValues[1];

	beta2 = ( temp * temp2 + temp3 ) * density_ ;


//         // New rhs for cavitational fluid
//         // rho0 * 4/3 * pi* n * ( 3 * R^2 * d^2R/dt^2 + 6 * R * (dR/dt)^2) 
//         if ( solveStep_->GetActStep() == 1){ 

//           beta2 =density_*density_* 4/3*PI*bubbleDensity_*6*bubbleValues[0]
//             *bubbleValues[1]*bubbleValues[1]; 

//           //Vereinfachte RHS Siehe Commander
//           beta2Com=0.0;
//         }

//         else {
//           StdVector<Double> dydt(2);

//           // dimensionless**************************************
//           y_[0] = bubbleValues[0] / initRadius_;
//           y_[1] = bubbleValues[1] / (sqrt( pStatic_/ density_));
//           t_    = solveStep_->GetActTime()  / initRadius_ * (sqrt(pStatic_/ density_));
//           ptBubble_[couplIndex]->CompDeriv(t_, y_, dydt);
//           dydt[1] = dydt[1] * pStatic_ / ( density_ * initRadius_ );
//           beta2 =density_*density_*4/3*PI*bubbleDensity_*
//             (6*bubbleValues[0]*bubbleValues[1]*bubbleValues[1]
//              + 3*bubbleValues[0]*bubbleValues[0]*dydt[1] ); 

//           //Vereinfachte RHS Siehe Commander
//           beta2Com= density_*density_*4.0*PI*bubbleDensity_*initRadius_ *initRadius_*dydt[1];

//           // dimensionless**************************************

//           //Normal case+++++++++++++++++++++++++++++++++++++++
//           //	  ptBubble_[numEl]->CompDeriv(actTime_, bubbleValues, dydt);
//           //	  beta2 =density_*density_* 4/3*PI*bubbleDensity_*
//           //	    (6*bubbleValues[0]*bubbleValues[1]*bubbleValues[1]
//           //	     + 3*bubbleValues[0]*bubbleValues[0]*dydt[1] ); 
//           //Normal case+++++++++++++++++++++++++++++++++++++++
	  
//           Rpp=dydt[1];
//         }


	if ( tDim_  < 1.0 / frequency_){
	  if (couplIndex == 0  )
	    Info->PrintF( pdename_, "Partial coupling in RHS %e \n",tDim_);

	  rhsForm->SetFactor("0.0");
	}
	else{
	  if (couplIndex == 0  )
	    Info->PrintF( pdename_, "Full coupling in RHS %e \n",tDim_);

	  rhsForm->SetFactor(GenStr(beta2));
	}

        Vector<Double> elemVec, helpVec;
        rhsForm->CalcElemVector( elemVec, it );

        // Get indices in nodes-vector for the node numbers of the current element
        for( UInt iNode = 0; iNode < connecth.GetSize(); iNode++ ) {
          Integer nodePos = nodes.Find( connecth[iNode] );
          values[nodePos] += elemVec[iNode];
        }


        // store element result
        if ( writeRHS_ == true || saveElemRHSHist_.GetSize() > 0 ||  saveElemRHSHistRegion_.GetSize() > 0) {
          helpVec.Resize( 2 );
          helpVec[0] = beta2;
          helpVec[1] = Rpp;
          //helpVec[2] = beta2Com;
          UInt locElemNum = eqnMap_->Mesh2PdeElem(it.GetElem()->elemNum);
          addElemResult_.SetElemResult(locElemNum-1, helpVec);
	  //	  if( locElemNum == 2 ){
	    //  std::cout<<tDim_ << "   "<< beta2<< "   "<<Rpp<<std::endl;
	  //}
        }
      }
    }
    // delete rhs integrator
    delete rhsForm;
    
    // print values on screen
    //std::cerr << "BubbleRHSvalues:\n" << values << std::endl << std::endl;
  }


  void BubblePDE::CalcBubbleRadius( StdVector<Elem*> & elems,
				    Vector<Double>& couplVals,
				    SolutionType solType ) {



    //	std::cerr<<"bubble calcradius radiussize "<<radius_.GetSize()<<std::endl;
    // Initialize coupling values
    couplVals.Init();

    if ( solType == BUBBLE_RADIUS ) {

      for (UInt iElem = 0; iElem < elems.GetSize(); iElem++ ) {

        // Find correct index for this element in vector couplElems_
        Integer couplIndex = -1;
        
        for( UInt iCouplElem = 0; iCouplElem < couplElems_->GetSize(); iCouplElem++ ) {
          if( (*couplElems_)[iCouplElem]->elemNum == elems[iElem]->elemNum ) {
            couplIndex = iCouplElem;
            break;
          }
        }

        if( couplIndex == -1 ) {
          EXCEPTION( "Found no matching element for coupling element Nr" 
                     << elems[iElem]->elemNum << "!" );
        }
        // store back to vector
        couplVals[iElem] = radius_[couplIndex];
      }

    } else {
      for (UInt iElem = 0; iElem < elems.GetSize(); iElem++ ) {

        // Find correct index for this element in vector couplElems_
        Integer couplIndex = -1;
        
        for( UInt iCouplElem = 0; iCouplElem < couplElems_->GetSize(); iCouplElem++ ) {
          if( (*couplElems_)[iCouplElem]->elemNum == elems[iElem]->elemNum ) {
            couplIndex = iCouplElem;
            break;
          }
        }

        if( couplIndex == -1 ) {
          EXCEPTION( "Found no matching element for coupling element Nr" 
                     << elems[iElem]->elemNum << "!" );
        }

        // store back to vector
        couplVals[iElem] = velocity_[couplIndex];
      }
    }


  }  


  bool BubblePDE::HasOutput(SolutionType output)
  {
  
    switch (output) {
    case ACOU_BUBBLE_RHS_VAL:
      std::cerr << "BubblePDE: Return true for for ACOU_BUBBLE_RHS_VAL\n";
      return true;
      break;
    case BUBBLE_RADIUS:
      std::cerr << "BubblePDE: Return true for for BUBBLE_RADIUS\n";
      return true;
      break;
    case BUBBLE_RADIUS_DERIV_1:
      std::cerr << "BubblePDE: Return true for for BUBBLE_RADIUS_DERIV_1\n";
      return true;
      break;
    default:
      return false;
      break;
    }
    return false;
  }

  void BubblePDE::DefineAvailResults() {
    
    EXCEPTION( "Implement -> Contact Andreas!" );
  }

//   // ***********************************************************************
//   //   Obtain information on desired output quantities from parameter file
//   // ***********************************************************************
//   void BubblePDE::ReadStoreResults() {


//     // Construct vectors for restricted parameter search
//     StdVector<std::string> keyVec;
//     StdVector<std::string> attrVec;
//     StdVector<std::string> valVec;
//     StdVector<std::string> regionNames;
//     std::string quantity;

//     // *****************************
//     // Determine element results
//     // *****************************

//     keyVec  = pdename_, "storeResults", "elemResults", "region";
//     attrVec = "", "", "type";

//     // --- Main bubble Result -> MagFluxDensity  ---
//     quantity = "bubbleValues";
//     valVec  = "", "", quantity;
//     params->GetList( keyVec, attrVec, valVec, regionNames );

//     if ( regionNames.GetSize() > 0 ) {
//       hasOutput_ = true;
//       writeValues_ = true;
//       Info->PrintF( pdename_, " Storing bubbleValues to file 1\n" );

//       elemResult_.SetNumSolutions(1);
//       elemResult_.SetSolutionType(MAG_FLUX_DENSITY);
//       elemResult_.SetNumElems(eqnMap_->GetNumLocalElems());
//       elemResult_.SetNumDofs(3);
//       elemResult_.SetPtrEQNData(eqnMap_.get(), ptgrid_);
//       elemResult_.Init(); 

//       keyVec  = pdename_, "storeResults", "elemHistory", "saveRegion";
//       attrVec = "", "", "type";
//       quantity = "bubbleValues";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemBubbleHistRegion_ );

//       if ( saveElemBubbleHistRegion_.GetSize() > 0 ) {
// 	Info->PrintF( pdename_, " Storing bubbleValues to file 2 \n" );
//         for ( UInt k = 0; k < saveElemBubbleHistRegion_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemBubbleHistRegion_[k].c_str() );
//         }
//       }

//       keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
//       attrVec = "", "", "type";
//       quantity = "bubbleValues";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemBubbleHist_ );

//       if ( saveElemBubbleHist_.GetSize() > 0 ) {
// 	Info->PrintF( pdename_, " Storing bubbleValues to file 3 \n" );
//         for ( UInt k = 0; k < saveElemBubbleHist_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemBubbleHist_[k].c_str() );
//         }
//       }

//     }
//     else{

//       keyVec  = pdename_, "storeResults", "elemHistory", "saveRegion";
//       attrVec = "", "", "type";
//       quantity = "bubbleValues";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemBubbleHistRegion_ );

//       keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
//       attrVec = "", "", "type";
//       quantity = "bubbleValues";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemBubbleHist_ );

//       if ( saveElemBubbleHist_.GetSize() > 0 || saveElemBubbleHistRegion_.GetSize() > 0  ) {
// 	hasOutput_ = true;
// 	Info->PrintF( pdename_, " Storing bubbleValues to file 4 \n" );
//         for ( UInt k = 0; k < saveElemBubbleHistRegion_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemBubbleHistRegion_[k].c_str() );
//         }
//         for ( UInt k = 0; k < saveElemBubbleHist_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemBubbleHist_[k].c_str() );
//         }
      
// 	elemResult_.SetNumSolutions(1);
// 	elemResult_.SetSolutionType(MAG_FLUX_DENSITY);
// 	elemResult_.SetNumElems(eqnMap_->GetNumLocalElems());
// 	elemResult_.SetNumDofs(3);
// 	elemResult_.SetPtrEQNData(eqnMap_.get(), ptgrid_);
// 	elemResult_.Init(); 
//       }
//     }

//     // --- RHS bubble Result -> ElecFieldIntensity  ---
//     keyVec  = pdename_, "storeResults", "elemResults", "region";
//     attrVec = "", "", "type";
//     quantity = "bubbleRHS";
//     valVec  = "", "", quantity;
//     params->GetList( keyVec, attrVec, valVec, regionNames );

//     if ( regionNames.GetSize() > 0 ) {
//       hasOutput_ = true;
//       writeRHS_ = true;
//       Info->PrintF( pdename_, " Storing bubbleRHS to file 5\n" );

//       addElemResult_.SetNumSolutions(1);
//       addElemResult_.SetSolutionType(ELEC_FIELD_INTENSITY);
//       addElemResult_.SetNumElems(eqnMap_->GetNumLocalElems());
//       addElemResult_.SetNumDofs(2);
//       addElemResult_.SetPtrEQNData(eqnMap_.get(), ptgrid_);
//       addElemResult_.Init(); 

//       keyVec  = pdename_, "storeResults", "elemHistory", "saveRegion";
//       attrVec = "", "", "type";
//       quantity = "bubbleRHS";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemRHSHistRegion_ );

//       if ( saveElemRHSHistRegion_.GetSize() > 0 ) {
// 	Info->PrintF( pdename_, " Storing bubbleRHS to file 6\n" );
//         for ( UInt k = 0; k < saveElemRHSHistRegion_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemRHSHistRegion_[k].c_str() );
//         }
//       }
//       keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
//       attrVec = "", "", "type";
//       quantity = "bubbleRHS";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemRHSHist_ );

//       if ( saveElemRHSHist_.GetSize() > 0 ) {
// 	Info->PrintF( pdename_, " Storing bubbleRHS to file 7\n" );
//         for ( UInt k = 0; k < saveElemRHSHist_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemRHSHist_[k].c_str() );
//         }
//       }


//     }
//     else{

//       keyVec  = pdename_, "storeResults", "elemHistory", "saveRegion";
//       attrVec = "", "", "type";
//       quantity = "bubbleRHS";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemRHSHistRegion_ );

//       keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
//       attrVec = "", "", "type";
//       quantity = "bubbleRHS";
//       valVec  = "", "", quantity;
//       params->GetList( keyVec, attrVec, valVec, saveElemRHSHist_ );
      
//       if ( saveElemRHSHist_.GetSize() > 0 || saveElemBubbleHistRegion_.GetSize() > 0  ) {
// 	hasOutput_ = true;
// 	Info->PrintF( pdename_, " Storing bubbleRHS to file 8\n" );
//         for ( UInt k = 0; k < saveElemRHSHistRegion_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemRHSHistRegion_[k].c_str() );
//         }
// 	Info->PrintF( pdename_, " Storing bubbleRHS to file 9\n" );
//         for ( UInt k = 0; k < saveElemRHSHist_.GetSize(); k++ ) {
//           Info->PrintF( pdename_, "  %s\n", saveElemRHSHist_[k].c_str() );
//         }

// 	addElemResult_.SetNumSolutions(1);
// 	addElemResult_.SetSolutionType(ELEC_FIELD_INTENSITY);
// 	addElemResult_.SetNumElems(eqnMap_->GetNumLocalElems());
// 	addElemResult_.SetNumDofs(2);
// 	addElemResult_.SetPtrEQNData(eqnMap_.get(), ptgrid_);
// 	addElemResult_.Init();
//       }
//     }
    
//   }
  

} // end of namespace

