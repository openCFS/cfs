// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <complex>

#include "algsys/entrymanipulatordouble.hh"
#include "algsys/baseidbchandler.hh"
#include "algsys/standardsys.hh"
#include "matvec/matvec.hh"

#include "matvec/vector.hh"


namespace OLAS {


  // ********************
  //   SetElementMatrix
  // ********************
  void EntryManipulatorDouble::
  SetElementMatrix( FEMatrixType matrixID,
                    const PdeIdType pdeID1,
                    const PdeIdType pdeID2,
                    StdMatrix *stdMat,
                    StdMatrix *counterPart,
                    BaseIDBC_Handler *idbcHandler,
                    Double *elemMat,
                    Integer *connect1, UInt length1,
                    Integer *connect2, UInt length2,
                    UInt limit1, UInt limit2,
                    bool setCounterPart ) {

    ENTER_IFCN( "EntryManipulatorDouble::SetElementMatrix" );

    // Clear the arrays
    rowIndList_.clear();
    colIndList1_.clear();
    colIndList2_.clear();

    UInt aux = 0;

    // STEP 1: Generate row index list from first connect array, dropping
    //         equation numbers for dofs fixed by (in)homogeneous Dirichlet
    //         boundary conditions and changing the sign of those fixed by
    //         constraints.
    for ( UInt i = 0; i < length1; i++ ) {
      aux = (UInt)std::abs( connect1[i] );
      if ( aux <= limit1 && aux > 0 ) {

        // Store equation number
        rowIndList_.push_back( aux );

        // Store row index into local element matrix
        rowIndList_.push_back( i );
      }
    }

    // STEP 2: Split the second connect array into two edge lists, one for
    //         the graph and one for the IDBCgraph (which handles the dofs
    //         fixed by inhomogeneous Dirichlet boundary conditions)
    for ( UInt i = 0; i < length2; i++ ) {
      aux = (UInt)std::abs( connect2[i] );
      if ( aux > 0 ) {
        if ( aux > limit2 ) {

          // Store equation number and column index into
          // local element matrix
          // colIndList2_.push_back( aux );
          colIndList2_.push_back( aux - limit2 );
          colIndList2_.push_back(  i  );
        }
        else {

          // Store equation number and column index into
          // local element matrix
          colIndList1_.push_back( aux );
          colIndList1_.push_back(  i  );
        }
      }
    }

#ifdef DEBUG_ASSEMBLE
    // output original connectivity
    (*debug) << "\n ------------------------------------------------------"
             << "\n EntryManipulatorDouble::SetElementMatrix\n"
             << " pdeID1 = " << pdeID1 << '\n'
             << " pdeID2 = " << pdeID2 << '\n'
             << " limit1 = " << limit1 << '\n'
             << " limit1 = " << limit2 << '\n';
    (*debug) << " connect1 ";
    for ( UInt i = 0; i < length1; i++ ) {
      (*debug) << connect1[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " connect2 ";
    for ( UInt i = 0; i < length2; i++ ) {
      (*debug) << connect2[i] << " ";
    }
    (*debug) << std::endl;

    // output new connectivity
    (*debug) << " rowIndList_: ";
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      (*debug) << rowIndList_[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " colIndList1_: ";
    for ( UInt i = 0; i < colIndList1_.size(); i += 2 ) {
      (*debug) << colIndList1_[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " colIndList2_: ";
    for ( UInt i = 0; i < colIndList2_.size(); i += 2 ) {
      (*debug) << colIndList2_[i] << " ";
    }
    (*debug) << std::endl;
#endif


    // STEP 3a: Insert values into matrix for real dofs
    UInt rowInd;
    UInt colInd;
#ifdef DEBUG_ASSEMBLE
    (*debug) << "\n free <-> free matrix:\n";
#endif
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
        colInd = colIndList1_[j+1];

#ifdef DEBUG_ASSEMBLE
        (*debug) << " mat(" << rowIndList_[i] << ", "
                 << colIndList1_[j] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
#endif

        stdMat->AddToMatrixEntry( rowIndList_[i], colIndList1_[j],
                                  elemMat[ rowInd * length2 + colInd ] );

      }
    }


    // STEP 3b: Insert values of transpose counterpart into matrix
    if ( setCounterPart == true ) {

#ifdef DEBUG_ASSEMBLE
      (*debug) << "\n free <-> free matrix (counterpart):\n";
#endif

      for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
        rowInd = rowIndList_[i+1];
        for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
          colInd = colIndList1_[j+1];

#ifdef DEBUG_ASSEMBLE
          (*debug) << " mat(" << colIndList1_[j] << ", "
                   << rowIndList_[i] << ") += "
                   << elemMat[ rowInd * length2 + colInd ] << std::endl;
#endif
          counterPart->AddToMatrixEntry( colIndList1_[j], rowIndList_[i],
                                         elemMat[ rowInd * length2 + colInd]);
        }
      }
    }

    // STEP 4: Insert values for fixed dofs into IDBC_Handler object
#ifdef DEBUG_ASSEMBLE
    (*debug) << "\n free <-> fixed matrix:\n";
#endif
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList2_.size(); j += 2 ) {
        colInd = colIndList2_[j+1];

#ifdef DEBUG_ASSEMBLE
        (*debug) << " mat(" << rowIndList_[i] << ", "
                 << colIndList2_[j] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
#endif
        idbcHandler->AddWeightFixedToFree( matrixID, pdeID1, pdeID2,
                                           rowIndList_[i], colIndList2_[j],
                                           elemMat[ rowInd * length2
                                                    + colInd ] );
      }
    }

#ifdef DEBUG_ASSEMBLE
    (*debug) << "------------------------------------------------------"
             << std::endl;
#endif

  }


  // **********************
  //   SetCounterPartOnly
  // **********************
  void EntryManipulatorDouble::
  SetCounterPartOnly( FEMatrixType matrixID,
                      const PdeIdType pdeID1,
                      const PdeIdType pdeID2,
                      StdMatrix *stdMat,
                      StdMatrix *counterPart,
                      BaseIDBC_Handler *idbcHandler,
                      Double *elemMat,
                      Integer *connect1, UInt length1,
                      Integer *connect2, UInt length2,
                      UInt limit1, UInt limit2 ) {

    ENTER_IFCN( "EntryManipulatorDouble::SetCounterPartOnly" );

    // Clear the arrays
    rowIndList_.clear();
    colIndList1_.clear();
    colIndList2_.clear();

    UInt aux = 0;

    // STEP 1: Generate row index list from first connect array, dropping
    //         equation numbers for dofs fixed by (in)homogeneous Dirichlet
    //         boundary conditions and changing the sign of those fixed by
    //         constraints.
    for ( UInt i = 0; i < length1; i++ ) {
      aux = (UInt)std::abs( connect1[i] );
      if ( aux <= limit1 && aux > 0 ) {

        // Store equation number
        rowIndList_.push_back( aux );

        // Store row index into local element matrix
        rowIndList_.push_back( i );
      }
    }

    // STEP 2: Split the second connect array into two edge lists, one for
    //         the graph and one for the IDBCgraph (which handles the dofs
    //         fixed by inhomogeneous Dirichlet boundary conditions)
    for ( UInt i = 0; i < length2; i++ ) {
      aux = (UInt)std::abs( connect2[i] );
      if ( aux > 0 ) {
        if ( aux > limit2 ) {

          // Store equation number and column index into
          // local element matrix
          // colIndList2_.push_back( aux );
          colIndList2_.push_back( aux - limit2 );
          colIndList2_.push_back(  i  );
        }
        else {

          // Store equation number and column index into
          // local element matrix
          colIndList1_.push_back( aux );
          colIndList1_.push_back(  i  );
        }
      }
    }

#ifdef DEBUG_ASSEMBLE
    // output original connectivity
    (*debug) << "\n ------------------------------------------------------"
             << "\n EntryManipulatorDouble::SetCounterPartOnly\n"
             << " pdeID1 = " << pdeID1 << '\n'
             << " pdeID2 = " << pdeID2 << '\n'
             << " limit1 = " << limit1 << '\n'
             << " limit1 = " << limit2 << '\n';
    (*debug) << " connect1 ";
    for ( UInt i = 0; i < length1; i++ ) {
      (*debug) << connect1[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " connect2 ";
    for ( UInt i = 0; i < length2; i++ ) {
      (*debug) << connect2[i] << " ";
    }
    (*debug) << std::endl;

    // output new connectivity
    (*debug) << " rowIndList_: ";
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      (*debug) << rowIndList_[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " colIndList1_: ";
    for ( UInt i = 0; i < colIndList1_.size(); i += 2 ) {
      (*debug) << colIndList1_[i] << " ";
    }
    (*debug) << std::endl;
    (*debug) << " colIndList2_: ";
    for ( UInt i = 0; i < colIndList2_.size(); i += 2 ) {
      (*debug) << colIndList2_[i] << " ";
    }
    (*debug) << std::endl;
#endif


    // STEP 3: Insert values of transpose counterpart into matrix
    UInt rowInd;
    UInt colInd;

#ifdef DEBUG_ASSEMBLE
    (*debug) << "\n free <-> free matrix (counterpart):\n";
#endif

    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
        colInd = colIndList1_[j+1];

#ifdef DEBUG_ASSEMBLE
        (*debug) << " mat(" << colIndList1_[j] << ", "
                 << rowIndList_[i] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
#endif
        counterPart->AddToMatrixEntry( colIndList1_[j], rowIndList_[i],
                                       elemMat[ rowInd * length2 + colInd] );
      }
    }

    // STEP 4: Insert values for fixed dofs into IDBC_Handler object
#ifdef DEBUG_ASSEMBLE
    (*debug) << "\n free <-> fixed matrix:\n";
#endif
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList2_.size(); j += 2 ) {
        colInd = colIndList2_[j+1];

#ifdef DEBUG_ASSEMBLE
        (*debug) << " mat(" << rowIndList_[i] << ", "
                 << colIndList2_[j] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
#endif
        idbcHandler->AddWeightFixedToFree( matrixID, pdeID1, pdeID2,
                                           rowIndList_[i], colIndList2_[j],
                                           elemMat[ rowInd * length2
                                                    + colInd ] );
      }
    }

#ifdef DEBUG_ASSEMBLE
    (*debug) << "------------------------------------------------------"
             << std::endl;
#endif

  }


  // *****************
  //   SetElementRHS
  // *****************
  void EntryManipulatorDouble::SetElementRHS( StdVector *rhs,
                                              Double *elemRHS,
                                              Integer *connect, 
                                              UInt elemSize,
                                              UInt limit ) {

    ENTER_IFCN( "EntryManipulatorDouble::SetElementRHS" );

#ifdef DEBUG_ASSEMBLE
    (*debug) << "SetElementRHS (" << elemSize << "): ";
    for ( UInt i = 0; i < elemSize; i++ ) {
      (*debug) << "(" << connect[i] << "," << elemRHS[i] << ") ";
    }
    (*debug) << std::endl;
#endif

    for ( UInt i = 0; i < elemSize; i++ ) {

      // For constraints: Add values of slaves to master
      if ( connect[i] < 0 ) {
        rhs->AddToVectorEntry( -connect[i], elemRHS[i] );
      }

      // Insert values related to free degrees of freedom
      // Note: Cast to UInt is save, since connect[i] is
      // not negative
      else if ( connect[i] > 0 && (UInt)connect[i] <= limit ) {
        rhs->AddToVectorEntry( connect[i], elemRHS[i] );
      }
    }

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }
  
  // ***************************
  //   GetMatrixEntry (Double)
  // ***************************
  void EntryManipulatorDouble::
  GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Double & entry, Integer eqnNr1, 
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2 ) {
    ENTER_IFCN( "EntryManipulatorDouble::GetMatrixEntry" );
    
    Double dummy = 0.0;

    // Calculate absolute value of equation numbers (for constraints)
    UInt aux1 = abs(eqnNr1);
    UInt aux2 = abs(eqnNr2);
    
    // Check if entry is for a free or a fixed dof
    if ( eqnNr1 <= limit1 && aux1 > 0 ) {
      
      // check if entry depends on a free or a fixed dof
      if ( eqnNr2 <= limit2 && aux2 > 0 ) {
        stdMat->GetMatrixEntry( eqnNr1, eqnNr2, entry );
      } else {
        idbcHandler->GetWeightFixedToFree( matrix_id, pdeID1, pdeID2,
                                           eqnNr1, aux2 - limit2, 
                                           entry, dummy );
      } 
    }else {
      entry = 0.0;
    }
  }
  
  // ***************************
  //   GetMatrixEntry (Complex)
  // ***************************
  void EntryManipulatorDouble::
  GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Complex & entry, Integer eqnNr1, 
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2 ) {
    ENTER_IFCN( "EntryManipulatorDouble::GetMatrixEntry" );
    (*error) << "EntryManipulatorDouble::GetMatrixEntry: "
             << "Interface for complex not implemented, since this "
             << "makes no sense!";
    Error( __FILE__, __LINE__ );

  }

  // ***************************
  //   SetMatrixEntry (Double)
  // ***************************
  void EntryManipulatorDouble::
  SetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Double  entry, Integer eqnNr1, 
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2, bool setCounterPart ) {
    ENTER_IFCN( "EntryManipulatorDouble::SetMatrixEntry" );

    Double dummy = 0.0;
    
    // Calculate absolute value of equation numbers (for constraints)
    UInt aux1 = abs(eqnNr1);
    UInt aux2 = abs(eqnNr2);

    // Check if entry is for a free or a fixed dof 
    if ( eqnNr1 <= limit1 && aux1 > 0 ) {

      // check if entry depends on a free or a fixed dof
      if ( eqnNr2 <= limit2 && aux2 > 0 ) {
        stdMat->SetMatrixEntry( eqnNr1, eqnNr2, entry, setCounterPart );
      } else {
        idbcHandler->SetWeightFixedToFree( matrix_id, pdeID1, pdeID2,
                                           eqnNr1, aux2 - limit2, 
                                           entry, dummy );
      } 
    }else {
      (*error) << "Can not set an matrix entry for position eqnNr1 = "
               << eqnNr1 << ", eqnNr2 = " << eqnNr2 << " for pdeID = "
               << pdeID1 << "!\n";
      Error( __FILE__, __LINE__ );
    }

  }

  // ***************************
  //   SetMatrixEntry (Complex)
  // ***************************
  void EntryManipulatorDouble::
  SetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Complex entry, Integer eqnNr1, 
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2,bool setCounterPart ) {
    ENTER_IFCN( "EntryManipulatorDouble::SetMatrixEntry" );
    (*error) << "EntryManipulatorDouble::SetMatrixEntry: "
             << "Interface for complex not implemented, since this "
             << "makes no sense!";
    Error( __FILE__, __LINE__ );

  }

  // ***********************
  //   SetNodeRHS (Double)
  // ***********************
  void EntryManipulatorDouble::SetNodeRHS( StdVector *rhs, Double val,
                                           Integer node ) {

    ENTER_IFCN( "EntryManipulatorDouble::SetNodeRHS" );

    if ( node != 0 ) {
      rhs->AddToVectorEntry( std::abs(node), val );
    }

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }

  // ***********************
  //   SetNodeRHS (Complex)
  // ***********************
  void EntryManipulatorDouble::SetNodeRHS( StdVector *rhs, Complex val,
                                           Integer node ) {
    ENTER_IFCN( "EntryManipulatorDouble::SetNodeRHS" );
    (*error) << "EntryManipulatorDouble::SetNodeRHS: "
             << "Interface for complex not implemented, since this "
             << "makes no sense!";
    Error( __FILE__, __LINE__ );
  }


  // ******************
  //   SetVectorEntry
  // ******************
  void EntryManipulatorDouble::SetVectorEntry( StdVector *vec,
                                               UInt index,
                                               Double &newVal ) {
    ENTER_IFCN( "EntryManipulatorDouble::SetVectorEntry" );
    vec->SetVectorEntry( index, newVal );
  }

  void EntryManipulatorDouble::SetVectorEntry( StdVector *vec,
                                               UInt index,
                                               Complex &newVal ) {
    ENTER_IFCN( "EntryManipulatorDouble::SetVectorEntry" );
    (*error) << "EntryManipulatorDouble::SetVectorEntry: "
             << "Interface for complex not implemented, since this "
             << "makes no sense!";
    Error( __FILE__, __LINE__ );
  }


  // ***********
  //   InitRHS
  // ***********
  void EntryManipulatorDouble::InitRHS( StdVector *rhs,
                                        const Double *newRHS ) {

    ENTER_FCN( "EntryManipulatorDouble::InitRHS" );

    // For the standard vector class we use fast access
    if ( Vector<Double> *myRHS = dynamic_cast<Vector<Double>*>(rhs) ) {
      for ( UInt i = 1; i <= rhs->GetSize(); i++ ) {
        (*myRHS)[i] = (newRHS-1)[i];
      }
    }

    // otherwise we use the slow access
    // i.e. in the parallel case this doesn't hurt too much
    // since CFS does everything sequentially anyway
    //
    // What makes you think this is slower than the above?
    else {
      for ( UInt i = 1; i <= rhs->GetSize(); i++ ) {
        rhs->SetVectorEntry(i, (newRHS-1)[i]);
      }
    }
      
    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }


  // *******************************
  //   UpdateRHS (f = f + A * fup)
  // *******************************
  void EntryManipulatorDouble::UpdateRHS( StdVector *rhs, StdMatrix *stdMat,
                                          Double* fup ) { 

    ENTER_FCN( "EntryManipulatorDouble::UpdateRHS" );


    // Support for the Vector<Double> case was implemented for
    // HyreMatrix::MultAdd some time ago. Status not quite clear.
    // Commented this out on 21.12.2004 (Marcus)

    // Currently we can only do this with the Vector<T> class,
    // i.e. with OLAS/LAPACK matrices
    // MatrixStorageType matType = sysmat_[matrix_id]->GetStorageType();

    // Generate an empty vector object to use as container
    Vector<Double> UpdVec;

    // Associate data array with vector object, but retain
    // responsibility for memory management
    UpdVec.Replace( stdMat->GetNcols(), fup - 1, false );

    // Perform arithmetic update
    stdMat->MultAdd( UpdVec, *rhs );

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }


  // ======================================================================
  // METHODS FOR HANDLING IDBCs
  // ======================================================================


  // *********************
  //   AdaptSystemMatrix
  // *********************
  void EntryManipulatorDouble::AdaptSystemMatrix( StdMatrix &stdMat,
                                                  UInt *dirichletEQN,
                                                  UInt numIDBC,
                                                  Double &penaltyTerm ) {

    ENTER_IFCN( "EntryManipulatorDouble::AdaptSystemMatrix" );

    // Loop over all Dirichlet equation numbers and replace
    // diagonal matrix entry by penalty term (no components
    // required in this scalar case)
    for ( UInt i = 1; i <= numIDBC; i++ ) {
      stdMat.SetDiagEntry( dirichletEQN[i], penaltyTerm );
    }
  }


  // *******************
  //   AdaptRHSForIDBC
  // *******************
  void EntryManipulatorDouble::AdaptRHSForIDBC( StdVector &rhs,
                                                UInt *dirichletEQN,
                                                StdVector &dirichletValue,
                                                Double &penaltyTerm,
                                                UInt numIDBC ) {

    ENTER_IFCN( "EntryManipulatorDouble::AdaptRHSForIDBC" );

    Double entry;
      
    for( UInt i = 1; i <= numIDBC; i++ ) {
      dirichletValue.GetVectorEntry( i, entry );
      entry *= penaltyTerm;
      rhs.SetVectorEntry( dirichletEQN[i], entry );
    }
  }


  // ******************
  //   GetSolutionVal
  // ******************

  void EntryManipulatorDouble::GetSolutionVal( StdVector *sol, Double* &ptSol,
                                               const PdeIdType identifierPDE){

    ENTER_FCN( "EntryManipulatorDouble::GetSolutionVal" );

    // only update the solution buffer if it is invalid
    if ( solBufferIsValid_ == false ) {
      solBuffer_ = *sol;
      // set validation flag
      solBufferIsValid_ = true;
    } 

    // Adapt 1-based array to 0-based array
    // for CFS++
    Double * ptr = solBuffer_.GetPointer();
    ptr++;
    ptSol = ptr;
  }


  void
  EntryManipulatorDouble::GetSolutionVal( StdVector *sol, Complex* &ptSol,
                                          const PdeIdType identifierPDE ) {
    (*error) << "GetSolutionVal(Complex) not defined for real entries";
    Error( __FILE__, __LINE__ );
  }


  // *************
  //   GetRHSVal
  // *************

  // this function needs some special treatment for parallel case
  void EntryManipulatorDouble::GetRHSVal( StdVector *rhs, Double* &ptRhs,
                                          const PdeIdType identifierPDE ) {

    ENTER_FCN( "EntryManipulatorDouble::GetRHSVal" );

    // only update the rhs buffer if it is invalid
    if ( rhsBufferIsValid_ == false ) {

      rhsBuffer_ = *rhs;
    
      // set validation flag
      rhsBufferIsValid_ = true;    
    }

    // Adapt 1-based array to 0-based array
    // for CFS++
    ptRhs = rhsBuffer_.GetPointer();
    ptRhs++;
  }

  void EntryManipulatorDouble::GetRHSVal( StdVector *rhs, Complex* &ptRhs, 
                                             const PdeIdType identifierPDE ){
    (*error) << "GetRHSVal(Complex) not defined for real entries";
    Error( __FILE__, __LINE__ );
  }


  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  void EntryManipulatorDouble::AddToDiagMatrixEntry( StdMatrix *stdMat,
                                                     Integer eqnNum,
                                                     Double *val ) {

    ENTER_IFCN( "EntryManipulatorDouble::AddToDiagMatrixEntry" );

    // Perform consistency checks
#ifdef DEBUG_ASSEMBLE
    if ( val == NULL ) {
      (*error) << "EntryManipulatorDouble::AddToDiagMatrixEntry: "
               << "Parameter val is a NULL pointer";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Get hold of matrix entry and adapt it
    Double entry;
    stdMat->GetDiagEntry( std::abs(eqnNum), entry );
    entry += val[0];
    stdMat->SetDiagEntry( std::abs(eqnNum), entry );
  }

}
