// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "assemble.hh"

#include "Domain/entityList.hh"
#include "Forms/baseForm.hh"
#include "Forms/linearForm.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/harmonicDriver.hh"
#include "DataInOut/Logging/cfslog.hh"


namespace CoupledField {


  // declare logging stream
  DECLARE_LOG(assemble)
  DEFINE_LOG(assemble, "assemble")

  Assemble::Assemble( BaseSystem* algsys,
                      AnalysisType analysis,
                      UInt maxTimeDerivOrder ) {
    ENTER_FCN( "Assemble::Assemble" );
    
    // init general params
    algsys_ = algsys;
    analysisType_ = analysis;
    isFirstTime_ = true;
    matrixUpdated_ = false;
    maxTimeDerivOrder_ = maxTimeDerivOrder;

    // init damping specific data
    raylDampFactor_ = 1.0;
    
    // Calculate matrix map from general matrix types to analysis
    // specific ones
    CreateMatrixMap();

    // Initialize reassemble-map
    // Note: in the beginning all matrices have to be "re-"assembled
    // which means the destination matrices have to be cleared and
    // initialized. After the first "AssmebleMatrices" call the
    // correct values for reassembling are determined
    matReassemble_[STIFFNESS] = true;
    matReassemble_[DAMPING] = true;
    matReassemble_[MASS] = true;

    // Obtain a new mathParser handler
    mHandle_ = domain->GetMathParser()->GetNewHandle();
  }

  Assemble::~Assemble() {
    ENTER_FCN( "Assemble::~Assemble" );
    
    // Delete bilinear contexts
    std::set<BiLinFormContext*>::iterator biLinIt;
    
    for( biLinIt = biLinForms_.begin(); 
         biLinIt != biLinForms_.end(); biLinIt++ ) {
      delete *biLinIt;
    }

    // Delete linear contexts
    std::set<LinearFormContext*>::iterator linIt;
    
    for( linIt = linForms_.begin(); 
         linIt != linForms_.end(); linIt++ ) {
      delete *linIt;
    }
  }
 
  void Assemble::AddBiLinearForm( BiLinFormContext* biLinContext ) {
    ENTER_FCN( "Assemble::AddBiLinearForm" );
    
    LOG_DBG(assemble) << "Adding BiLinearForm '" 
                      << biLinContext->GetIntegrator()->GetName()
                      << "' on '" 
                      << biLinContext->GetFirstEntities()->GetName()
                      << "'";
    
    // assert that Integrator is set
    assert( biLinContext->GetIntegrator() != NULL );

    // assert that some entites are set
    assert( biLinContext->GetFirstEntities() != NULL );
    assert( biLinContext->GetSecondEntities() != NULL );

    // assert that the pdes are set
    assert( biLinContext->GetFirstPde() != NULL );
    assert( biLinContext->GetSecondPde() != NULL );

    // If the datatype of the bilinearformcontext is "COMPLEX"
    // we have to ensure that we are in an harmonic case.
    // Otherwise we issue an error
    if( (biLinContext->GetEntryType() == IMAG ||
         biLinContext->GetEntryType() == COMPLEX )
        && analysisType_ != HARMONIC ) {
      EXCEPTION( "Can not add integrator '" 
                 << biLinContext->GetIntegrator()->GetName() 
                 << "' with complex/imaginary entries for a "
                 << "non-harmonic analysis." );
    }
    
    FEMatrixType mappedFEType = matrixMap_[biLinContext->GetDestMat()];
    FEMatrixType mappedSecFEType = 
      matrixMap_[biLinContext->GetSecDestMat()];
    
    // Check if integrator can be assembled in this type of simulation
    if( mappedFEType != NOTYPE ) {

      // Store bilinear form
      biLinForms_.insert( biLinContext );

      // Pass needed matrix type to algebraic system
      PdeIdType id1 = biLinContext->GetFirstPde()->GetPDEId();
      PdeIdType id2 = biLinContext->GetFirstPde()->GetPDEId();
      algsys_->SetFEMatrixType( mappedFEType, id1, id2 );
                                
      // Check for secondary matrix type
      if( mappedSecFEType != NOTYPE ) {
        algsys_->SetFEMatrixType( mappedSecFEType, id1, id2 );
      }
    }
  }

  void Assemble::AddLinearForm( LinearFormContext* linContext ) {
    ENTER_FCN( "Assemble::AddLinearForm" );

    // assert that Integrator is set
    assert( linContext->GetIntegrator() != NULL );

    // assert that the pdes are set
    assert( linContext->GetPde() != NULL );

    // assert that some entites are set
    assert( linContext->GetEntities() != NULL );

    linForms_.insert( linContext );

  }

  void Assemble::AddLoads( LoadList& load ) {
    ENTER_FCN( "Assemble::AddLoad" );
    
    for( UInt i=0; i<load.GetSize(); i++ ) {
      loads_.Push_back( load[i] );
    }
  }

  void Assemble::SetupMatrixGraph(PdeIdType pdeId1, PdeIdType pdeId2 ) {
    ENTER_FCN( "Assemble::SetupMatrixGraph" );

    StdVector<Integer> eqnVec1, eqnVec2;
    PdeIdType id1, id2;
    
    // iterate over all descriptors
    std::set<BiLinFormContext*>::iterator formsIt;
    for ( formsIt = biLinForms_.begin(); 
          formsIt != biLinForms_.end(); formsIt++ ) {

      // get integrator
      BiLinFormContext & actContext = **formsIt;

      // Check if pdeIds of formsContext match
      if( actContext.GetFirstPde()->GetPDEId() == pdeId1 &&
          actContext.GetSecondPde()->GetPDEId() == pdeId2 ) {

        try {

          // get entity iterators
          EntityIterator it1 = actContext.GetFirstEntities()->GetIterator();
          EntityIterator it2 = actContext.GetSecondEntities()->GetIterator();
          UInt size = actContext.GetFirstEntities()->GetSize();
          it1.Begin();
          it2.Begin();
          
          // iterate over all entities
          for ( UInt i = 0; i < size; i++ ) {
            
            // Get equation numbers
            actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, id1, id2 );
            
            // Pass entity eqn-connectivity to algebraic system
            algsys_->
              SetElementPos( id1, eqnVec1.GetPointer(), eqnVec1.GetSize(), 
                             id2, eqnVec2.GetPointer(), eqnVec2.GetSize(), 
                             actContext.IsSetCounterPart());
            
            // increment iterators
            it1++;
            it2++;
          } // loop (entities)
        } catch (Exception& e) {
          RETHROW_EXCEPTION(e, "Could not setup matrix graph for "
                            << "BiLinearForm '"
                            << actContext.GetIntegrator()->GetName() << "' on '"
                            << actContext.GetFirstEntities()->GetName()<< "'" );
        }
        
      } // if
    } // loop (integrators)
  }

  void Assemble::AssembleMatrices() {
    ENTER_FCN( "Assemble::AssembleMatrices" );

    Matrix<Double> elemMatrix;
    StdVector<Integer> eqnVec1, eqnVec2;
    PdeIdType pdeId1, pdeId2;

    // Reset for matrix update
    matrixUpdated_ = false;
    
    // Temporary: Check each time for non-linearities
    CheckNonLinearities();

    // Init all matrices, which have to be reassembled
    std::map<FEMatrixType, bool>::iterator it;
    for( it = matReassemble_.begin(); it != matReassemble_.end(); it++ ) {
      if( it->second == true ) {
        algsys_->InitMatrix( matrixMap_[it->first] );
      }
    }
    
    // iterate over all descriptors
    std::set<BiLinFormContext*>::iterator formsIt;
    for ( formsIt = biLinForms_.begin(); 
          formsIt != biLinForms_.end(); formsIt++ ) {
      
      // get integrator
      BiLinFormContext & actContext = **formsIt;
      
      // get matrix destinations
      FEMatrixType destMat = actContext.GetDestMat(); 
      FEMatrixType secDestMat = actContext.GetSecDestMat();
      
      // If assemble was already called and the current destination
      // matrix must not be reassembled -> continue with next iterator
      if( isFirstTime_ == false && matReassemble_[destMat] == false ) {
        continue;
      }

      // Update flag
      matrixUpdated_ = true;

      BaseForm * form = actContext.GetIntegrator();

      try {

        // get entity iterators
        EntityIterator  it1 = actContext.GetFirstEntities()->GetIterator();
        EntityIterator  it2 = actContext.GetSecondEntities()->GetIterator();
        UInt size = actContext.GetFirstEntities()->GetSize();
        it1.Begin();
        it2.Begin();

        // adjust damping factor 
        AdjustDamping( actContext );

        // iterate over all entities
        for ( UInt i=0; i<size; i++ ) {
          
          // Calc element matrix
          form->CalcElementMatrix( elemMatrix, it1, it2 );
          
          // Map equation numbers
          actContext.MapEqns( it1, it2, eqnVec1, eqnVec2, pdeId1, pdeId2 );

          // Pass element matrix to algebraic system (primary matrix)
          InsertMatrix( destMat, actContext, elemMatrix, eqnVec1, eqnVec2,
                        pdeId1, pdeId2);

          // Check for secondary matrix type
          if (secDestMat != NOTYPE ) {
          
            Double dampFactor = 1.0;
            if ( actContext.getPtDamplayer() != NULL ) {
              actContext.getPtDamplayer()->CalcDampFactor(dampFactor, it1);
              // std::cout << "   dampFactor: " <<  dampFactor << std::endl;
            }

            // Rayleigh damping
            elemMatrix *= raylDampFactor_ * dampFactor * actContext.GetSecMatFac();

            // Pass secondary matrix part to algebraic system
            InsertMatrix( secDestMat, actContext, elemMatrix, eqnVec1, eqnVec2,
                          pdeId1, pdeId2); 
          }

          // increment iterators
          it1++;
          it2++;
        }
      } catch (Exception& e) {
        RETHROW_EXCEPTION(e, "Could not calculate element matrix of "
                          << "BiLinearForm '"
                          << form->GetName() << "' on '"
                          << actContext.GetFirstEntities()->GetName()<< "'" );
      }
    }
    
    // Change flag
    isFirstTime_ = false;
  }
  
  void Assemble::AssembleLinRHS( Double actTimeFreq){
    ENTER_FCN( "Assemble::AssembleLinRHS" );
  
    AssembleRHSLinForms( actTimeFreq, false );
    AssembleRHSLoads( actTimeFreq );
  }
  
  void Assemble::AssembleNonLinRHS( Double actTimeFreq) {
    ENTER_FCN( "Assemble::AssembleNonLinRHS" );
    AssembleRHSLinForms( actTimeFreq, true );
  }

  void Assemble::AssembleRHSLinForms( Double actTimeFreq, bool nonLin ) {
    ENTER_FCN( " Assemble::AssembleRHSLinForms" );
    
    StdVector<Integer> eqnVec;
    PdeIdType pdeId;
    std::set<LinearFormContext*>::iterator formsIt;

    // iterate over all descriptors
    for ( formsIt = linForms_.begin(); 
          formsIt != linForms_.end(); 
          formsIt++ ) {
      
      // get integrator
      LinearFormContext & actContext = **formsIt;

      // Check, if lin/non-lin type of Context matches parameter nonLin
      if( actContext.IsNonLin() != nonLin ) {
        continue;
      }

      LinearForm * form = actContext.GetIntegrator();

      try {

        // get entity iterator
        EntityIterator  entIt = actContext.GetEntities()->GetIterator();

        if ( analysisType_ == HARMONIC ) {

          Vector<Complex> elemVec;
          for ( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {

            // Calculate complex valued element vector
            form->CalcElemVector( elemVec, entIt );

            // Map equation numbers
            actContext.MapEqns( entIt, eqnVec, pdeId );

            // We copy the complex vector 'elemVec' to a double vector
            //  'harmVec' which has real and imaginary part in consecutive 
            //  order and has double length
            Vector<Double> harmVec;
            Integer size = elemVec.GetSize();
            harmVec.Resize(2*size);
    
            Integer k=0;
            //real part
            for (Integer i=0; i<size; i++) {
              harmVec[k] = elemVec[i].real();
              k++;
            }
            //imaginary part
            for (Integer i=0; i<size; i++) {
              harmVec[k] = elemVec[i].imag();
              k++;
            }
          
            algsys_-> SetElementRHS( harmVec.GetPointer(), 
                                     pdeId, eqnVec.GetPointer(),
                                     eqnVec.GetSize() );
          }

        } else {

          // That should be STATIC, TRANSIENT or EIGENFREQUENCY

          Vector<Double> elemVec;
          // iterate over all entities
          for ( entIt.Begin(); !entIt.IsEnd(); entIt++ ) {
          
            // Calculate real valued element vector
            form->CalcElemVector( elemVec, entIt );
          
            // Map equation numbers
            actContext.MapEqns( entIt, eqnVec, pdeId );
          
            // Pass element vector to algebraic system
            algsys_-> SetElementRHS( elemVec.GetPointer(), 
                                     pdeId, eqnVec.GetPointer(),
                                     eqnVec.GetSize() );
          }

        }
      } catch (Exception& e) {
        RETHROW_EXCEPTION(e, "Could not calculate RHS vector for LinearForm '"
                          << form->GetName() << "' on '"
                          << actContext.GetEntities()->GetName() << "'" );
      }
    }  
  }

  void Assemble::AssembleRHSLoads( Double actTimeFreq ) {
    ENTER_FCN( "Assemble::AssembleRHSLoads" );

    Vector<Double> globCoord;
    Integer eqnNr = 0;
    Double phase = 0.0;
    Double val = 0.0;

    // get global math parser and pointer to grid
    MathParser * parser = domain->GetMathParser();
    CoordSystem * coosy = domain->GetCoordSystem();
    Grid * ptGrid = domain->GetGrid();

    // iterate over all load node-lists
    for( UInt iLoad=0; iLoad < loads_.GetSize(); iLoad++ ) {

      // get current load and its equation map
      LoadBc & actLoad = *loads_[iLoad];
      EqnMap & eqnMap = *(actLoad.eqnMap);

      try {
        // Obtain iterator
        EntityIterator it = actLoad.entities->GetIterator();
        
        for( it.Begin(); !it.IsEnd(); it++ ) {
          
          // If iterator points to a node, pass the current coordinate
          // to the parser
          if( it.GetType() == EntityList::NODE_LIST ) {
            
            // Get coordinates of node
            ptGrid->GetNodeCoordinate( globCoord, it.GetNode() );
            parser->SetCoordinates( mHandle_, *coosy, globCoord );
          }
          
          // Evaluate load value
          parser->SetExpr( mHandle_, actLoad.value );
          val = parser->Eval(mHandle_ );
          
          // Obtain equation number
          eqnNr = eqnMap.GetEqn( *(actLoad.result), it, actLoad.dof );
          
          if (analysisType_ == HARMONIC) {
            parser->SetExpr( mHandle_, actLoad.phase  );
            phase = parser->Eval( mHandle_ );
            Complex complexValue( val * cos( phase / 180 * PI ),
                                  val * sin( phase / 180 * PI ) );
            
            algsys_->SetNodeRHS(complexValue, eqnMap.GetPdeId(), eqnNr, 1);    
          } else {
            algsys_->SetNodeRHS(val, eqnMap.GetPdeId(), eqnNr, 1);    
          }
          
          
        } // nodes
      } catch (Exception& e) {
        RETHROW_EXCEPTION(e, "Could not set RHS load with value '"
                          << actLoad.value << "' on '"
                          << actLoad.entities->GetName() << "'" );
      }
    } // loads
    
  }
  
  void Assemble::PrintInfo( std::ostream& out ) {
    ENTER_FCN( "Assemble::PrintInfo" );

    out << "=================================\n"
        << "  Assemble: List of Integrators  \n"
        << "=================================\n\n";

    out << "I) Matrix BiLinearForms\n"
        << "-----------------------\n\n";

    out << std::setw(20) << "integrator name" << " | "
        << std::setw(15) << "region" << " | "
        << std::setw(10) << "matrix type" << "\n";
    out << std::setw(51) << std::setfill('-') << "" 
        << std::setfill(' ') << std::endl;

    // iterate over all descriptors
    std::set<BiLinFormContext*>::iterator it;
    for ( it = biLinForms_.begin(); it != biLinForms_.end(); it++ ) {
      
      // get integrator
      BiLinFormContext & context = **it;
      
      // integrator name
      out << std::setw(20) << context.GetIntegrator()->GetName() << " | ";
      
      // region name of entity list
      std::string regionName;
      if( context.GetFirstEntities()->GetType() == EntityList::ELEM_LIST ) {
        shared_ptr<ElemList> list = 
          boost::dynamic_pointer_cast<ElemList,EntityList> 
          (context.GetFirstEntities());
        regionName = domain->GetGrid()->RegionIdToName( list->GetRegion() );
      } else if ( context.GetFirstEntities()->GetType()
                  == EntityList::SURF_ELEM_LIST ) {
        shared_ptr<SurfElemList> list = 
          boost::dynamic_pointer_cast<SurfElemList,EntityList> 
          (context.GetFirstEntities());
        regionName = domain->GetGrid()->RegionIdToName( list->GetRegion() );
      }
      out << std::setw(15) << regionName << " | " ;
      
      // matrix type
      out << std::setw(10) << OLAS::Enum2String(context.GetDestMat() );
      
      out << std::endl;
    }


    out << "\n\nII) RHS LinearForms\n"
        << "-----------------------\n\n";

    out << std::setw(20) << "integrator name" << " | "
        << std::setw(15) << "region" << "\n";
    out << std::setw(38) << std::setfill('-') << "" 
        << std::setfill(' ') << std::endl;

    // iterate over all descriptors
    std::set<LinearFormContext*>::iterator linIt;
    for ( linIt = linForms_.begin(); linIt != linForms_.end(); linIt++ ) {
      
      // get integrator
      LinearFormContext & context = **linIt;
      
      // integrator name
      out << std::setw(20) << context.GetIntegrator()->GetName() << " | ";
      
      // region name of entity list
      std::string regionName;
      if( context.GetEntities()->GetType() == EntityList::ELEM_LIST ) {
        shared_ptr<ElemList> list = 
          boost::dynamic_pointer_cast<ElemList,EntityList> 
          (context.GetEntities());
        regionName = domain->GetGrid()->RegionIdToName( list->GetRegion() );
      } else if ( context.GetEntities()->GetType()
                  == EntityList::SURF_ELEM_LIST ) {
        shared_ptr<SurfElemList> list = 
          boost::dynamic_pointer_cast<SurfElemList,EntityList> 
          (context.GetEntities());
        regionName = domain->GetGrid()->RegionIdToName( list->GetRegion() );
      }
      out << std::setw(15) << regionName ;
      
      out << std::endl;
    }
    out << "\n\n";

  }

  void Assemble::CheckNonLinearities() {
    ENTER_FCN( "Assemble::CheckNonLinearities" );

    // Clear reassemble mat
    matReassemble_.clear();

    // iterate over all bilinearforms
    std::set<BiLinFormContext*>::iterator it;

    for( it = biLinForms_.begin(); it != biLinForms_.end(); it++ ) {
      BiLinFormContext & actContext = **it;

      if( actContext.IsNonLin() == true 
          || analysisType_ == HARMONIC ) {
        matReassemble_[actContext.GetDestMat()] = true;
        if ( actContext.GetSecDestMat() != NOTYPE ) {
          matReassemble_[actContext.GetSecDestMat()] = true; 
        }
      }
    }
    
    
  }

  
  void Assemble::Matrix2Harmonic(Vector<Double>& harmMat,
                                 Matrix<Double>& origMat,
                                 FEMatrixType matrixType,
                                 DataType entryType,
                                 Double omega ) {
    ENTER_FCN( "Assemble::Matrix2Harmonic" );
    
    Integer numRow = origMat.GetSizeRow();
    Integer numCol = origMat.GetSizeCol();
    harmMat.Resize(2*numRow*numCol);
    harmMat.Init();
    Integer k=0;

    if (entryType == REAL) {

      if (matrixType == STIFFNESS) {
        for (Integer row=0; row<numRow; row++)
          for (Integer col=0; col<numCol; col++) {
            harmMat[k] = origMat[row][col];
            k++;
          }
      }
    
      else if (matrixType == MASS) {

        // NOTE: here we have to distinguish, if the
        //  mass matrix is associated with first
        // or second order time derivative
        if( maxTimeDerivOrder_ == 2 ) {
          Double factor = -omega*omega;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        } else {
          Double factor = omega;
          k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        }
      }
      
      else if (matrixType == DAMPING) {       
        Double factor = omega;
        
        k=numRow*numCol;
        for (Integer row=0; row<numRow; row++)
          for (Integer col=0; col<numCol; col++) {
            harmMat[k] = factor*origMat[row][col];
            k++;
          }
      }
    } // end, if matatType == real...

    else if(entryType == IMAG){  // the "imaginary parts"
   
      if (matrixType == STIFFNESS) {
        k=numRow*numCol;
        for (Integer row=0; row<numRow; row++)
          for (Integer col=0; col<numCol; col++) {
            harmMat[k] = origMat[row][col];
            k++;
          }
      }
      
    
      else if (matrixType == MASS) {
        // NOTE: here we have to distinguish, if the
        //  mass matrix is associated with first
        // or second order time derivative
        if( maxTimeDerivOrder_ == 2 ) {
          Double factor = -omega*omega;
          k=numRow*numCol;
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = factor*origMat[row][col];
              k++;
            }
        } else {
          Double factor = omega;
          k=0;  
          for (Integer row=0; row<numRow; row++)
            for (Integer col=0; col<numCol; col++) {
              harmMat[k] = -factor*origMat[row][col];
              k++;
            }
        }
      }
      
      else if (matrixType == DAMPING) {
        Double factor = omega;
        k=0;  
        for (Integer row=0; row<numRow; row++)
          for (Integer col=0; col<numCol; col++) {
            harmMat[k] = -factor*origMat[row][col];
            k++;
          }
      }
      
    } // end if matType == imag
    
    else {
      (*error) <<"\n DataType" << entryType
               << "not specified "<<std::endl;
      Error( __FILE__, __LINE__ );
    }
    
  }
    

  void Assemble::CreateMatrixMap() {
    ENTER_FCN( "Assemble::CreateMatrixMap()" );

    // Dependent on the type of analysis, only certain matrix types
    // (SYSTEM, STIFFNESS, MASS, DAMPING, CONVECTION) are present.

    if( analysisType_ == STATIC ) {

      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = SYSTEM;
      matrixMap_[DAMPING]   = NOTYPE;
      matrixMap_[MASS]      = NOTYPE;

    } else if( analysisType_ == TRANSIENT ) {
      matrixMap_[SYSTEM]    = SYSTEM;
      matrixMap_[STIFFNESS] = STIFFNESS;
      matrixMap_[DAMPING]   = DAMPING;
      matrixMap_[MASS]      = MASS;
      
    } else if( analysisType_ == HARMONIC ) {
      matrixMap_[SYSTEM]    = SYSTEM; 
      matrixMap_[STIFFNESS] = SYSTEM;
      matrixMap_[DAMPING]   = SYSTEM;
      matrixMap_[MASS]      = SYSTEM;
      
    } else if( analysisType_ == EIGENFREQUENCY) {
      matrixMap_[SYSTEM]    = NOTYPE;
      matrixMap_[STIFFNESS] = STIFFNESS;
      matrixMap_[DAMPING]   = DAMPING;
      matrixMap_[MASS]      = MASS;
      
    } else {
      
      std::string analysisString;
      Enum2String( analysisType_, analysisString );
      (*error) << "Analysistype '" << analysisString << "' not known!";
      Error( __FILE__, __LINE__ );
    }
     
  }

  void Assemble::AdjustDamping( BiLinFormContext& context ) {
    ENTER_FCN( "Assemble::CheckDamping" );
    
    // Check, if damping matrix is present
    if( matReassemble_.find( DAMPING) == matReassemble_.end() )
      return;

    if( analysisType_ == HARMONIC ) {
      
      // obtain current frequency
      MathParser * parser = domain->GetMathParser();
      parser->SetExpr( mHandle_, "f" );
      Double actFreq = parser->Eval( mHandle_ );
      Double actOmega = actFreq * 2.0 * PI;
      
      // obtain matData-freq
      // NOTE: This mechanism should be changed in a way, that the damping
      // is part of the xml-material file, where tanDelta and the
      // related frequecy are specified
      
      HarmonicDriver * harmDriver = 
        dynamic_cast<HarmonicDriver*>(domain->GetSingleDriver() );

      Double matDataOmega = harmDriver->GetMatDataFreq() * 2.0 * PI;

      // get multiplicative pre factor depending on frequency
      if ( matDataOmega > 0 && actOmega > 0 ) {
        
        if ( context.GetDestMat() == STIFFNESS ) {
          raylDampFactor_ = matDataOmega / actOmega;
          Info->PrintF( "", " dampTransform for STIFFNESS ... %e\n",
                        raylDampFactor_ );
        }
        else if ( context.GetDestMat() == MASS ) {
          raylDampFactor_ = actOmega / matDataOmega;
          Info->PrintF( "", " dampTransform for MASS ........ %e\n",
                        raylDampFactor_ );
        }
      }
      
    }
  }

  bool Assemble::IsFEMatSymmetric( FEMatrixType feType ) {
    ENTER_FCN( " Assemble::IsFEMatSymmetric" );
    
    // Run over all bilinearform contexts
    std::map<FEMatrixType, bool> isSymmetric;
    std::set<BiLinFormContext*>::iterator it;

    // Assume at the beginning that all matrices are symmetric
    isSymmetric[SYSTEM] = true;
    isSymmetric[MASS] = true;
    isSymmetric[STIFFNESS] = true;
    isSymmetric[DAMPING] = true;
    
   
    // iterate over all bilinear forms
    for( it = biLinForms_.begin(); it != biLinForms_.end(); it++ ) {

      BiLinFormContext & actCt = (**it);
      
      // Check, where bilinearform gets assembled to
      if( (actCt.GetFirstPde() == actCt.GetSecondPde() )
          && (actCt.GetFirstResultInfo() == actCt.GetSecondResultInfo() )
          && (actCt.GetFirstEntities() == actCt.GetSecondEntities() ) ) {

        // Bilinearform gets assembled to main diagonal.
        // If bilinearform is non-symmetric, so is the related FE-matrix
        if( !actCt.GetIntegrator()->IsSymmetric()  ) {
          FEMatrixType mappedDest = matrixMap_[actCt.GetDestMat()];
          isSymmetric[mappedDest] = false;
          isSymmetric[SYSTEM] = false;
        }
      } else {

        // BiLinearform gets assembled to off-diagonal block.
        
        // If the bilinearorm is also assembled to the transposed block
        // we assume that the matrix still remains symmetric.
        // Otherwise we assume, that we need a non-symmetric matrix.
        if( !actCt.IsSetCounterPart() ) {
          FEMatrixType mappedDest = matrixMap_[actCt.GetDestMat()];
          isSymmetric[mappedDest] = false;
          isSymmetric[SYSTEM] = false;
        } 
      } 
    }

    // return flag for matrix of interest
    return isSymmetric[feType];
  }

   
  void Assemble::InsertMatrix( FEMatrixType dest, BiLinFormContext& context,
                               Matrix<Double>& elemMat, 
                               StdVector<Integer>& eqnVec1,
                               StdVector<Integer>& eqnVec2,
                               PdeIdType pdeId1, PdeIdType pdeId2) {
    ENTER_FCN( "Assemble::InsertMatrix" );
    Vector<Double> harmMat;

    // map original matrix destination to analysis-dependent one
    FEMatrixType mappedDest = matrixMap_[dest]; 
    
    if( mappedDest == NOTYPE ) {
      return;
    }
    
    if( analysisType_ == TRANSIENT 
        || analysisType_ == STATIC
        || analysisType_ == EIGENFREQUENCY) {

      algsys_->
        SetElementMatrix( mappedDest, elemMat.GetDataPointer(), 
                          pdeId1, eqnVec1.GetPointer(), eqnVec1.GetSize(), 
                          pdeId2, eqnVec2.GetPointer(), eqnVec2.GetSize(),
                          context.IsSetCounterPart() );
    } else {
      Double freq = context.GetFirstPde()->GetSolveStep()->GetActFreq();
      Double omega = freq * 2 * PI;
      Matrix2Harmonic( harmMat, elemMat, dest,
                       context.GetEntryType(), omega );
      algsys_->
        SetElementMatrix( mappedDest, harmMat.GetPointer(), 
                          pdeId1, eqnVec1.GetPointer(), eqnVec1.GetSize(), 
                          pdeId2, eqnVec2.GetPointer(), eqnVec2.GetSize(),
                          context.IsSetCounterPart() );
    }
  }
}
