// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <complex>
#include <assert.h>

#include "OLAS/algsys/baseentrymanipulator.hh"
#include "OLAS/algsys/baseidbchandler.hh"
#include "OLAS/algsys/standardsys.hh"

#include "MatVec/vector.hh"
#include "MatVec/stdmatrix.hh"
#include "MatVec/matrix.hh"


namespace CoupledField {


  // ********************
  //   SetElementMatrix
  // ********************
  template<typename T>
  void BaseEntryManipulator::
  SetElementMatrix( FEMatrixType matrixID,
                    const PdeIdType pdeID1,
                    const PdeIdType pdeID2,
                    StdMatrix *stdMat,
                    StdMatrix *counterPart,
                    BaseIDBC_Handler *idbcHandler,
                    const Matrix<T>& elemMat,
                    const StdVector<Integer>& connect1,
                    const StdVector<Integer>& connect2,
                    UInt limit1, UInt limit2,
                    bool setCounterPart) {

    // Clear the arrays
    rowIndList_.clear();
    colIndList1_.clear();
    colIndList2_.clear();

    UInt aux = 0;
    UInt length1 = connect1.GetSize();
    UInt length2 = connect2.GetSize();

    // STEP 1: Generate row index list from first connect array, dropping
    //         equation numbers for dofs fixed by (in)homogeneous Dirichlet
    //         boundary conditions and changing the sign of those fixed by
    //         constraints.
    for ( UInt i = 0; i < length1; i++ ) {
      aux = (UInt)std::abs( connect1[i] );
      if ( aux <= limit1 && aux > 0 ) {

        // Store equation number
        rowIndList_.push_back( aux - 1 );

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
          colIndList2_.push_back( aux - limit2 - 1 );
          colIndList2_.push_back(  i  );
        }
        else {

          // Store equation number and column index into
          // local element matrix
          colIndList1_.push_back( aux -1 );
          colIndList1_.push_back(  i  );
        }
      }
    }

    /*
    // output original connectivity
    (*debug) << "\n ------------------------------------------------------"
             << "\n BaseEntryManipulator::SetElementMatrix\n"
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
    */


    // STEP 3a: Insert values into matrix for real dofs
    UInt rowInd;
    UInt colInd;

    // (*debug) << "\n free <-> free matrix:\n";
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
        colInd = colIndList1_[j+1];

        /*(*debug) << " mat(" << rowIndList_[i] << ", "
                 << colIndList1_[j] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
        */
        stdMat->AddToMatrixEntry( rowIndList_[i], colIndList1_[j],
                                  elemMat[rowInd][colInd] );

      }
    }


    // STEP 3b: Insert values of transpose counterpart into matrix
    if ( setCounterPart == true ) {

      // (*debug) << "\n free <-> free matrix (counterpart):\n";

      for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
        rowInd = rowIndList_[i+1];
        for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
          colInd = colIndList1_[j+1];

          /*(*debug) << " mat(" << colIndList1_[j] << ", "
                   << rowIndList_[i] << ") += "
                   << elemMat[ rowInd * length2 + colInd ] << std::endl;
          */         
          counterPart->AddToMatrixEntry( colIndList1_[j], rowIndList_[i],
                                         elemMat[rowInd][colInd]);
        }
      }
    }

    // STEP 4: Insert values for fixed dofs into IDBC_Handler object
    // (*debug) << "\n free <-> fixed matrix:\n";
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList2_.size(); j += 2 ) {
        colInd = colIndList2_[j+1];

        /*(*debug) << " mat(" << rowIndList_[i] << ", "
                 << colIndList2_[j] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
        */
        idbcHandler->AddWeightFixedToFree( matrixID, pdeID1, pdeID2,
                                           rowIndList_[i], colIndList2_[j],
                                           elemMat[rowInd][colInd]);
      }
    }

    // (*debug) << "------------------------------------------------------"  << std::endl;

  }
  
  // explicit template instantiation
  template void BaseEntryManipulator::
  SetElementMatrix( FEMatrixType matrixID,
                    const PdeIdType pdeID1,
                    const PdeIdType pdeID2,
                    StdMatrix *stdMat,
                    StdMatrix *counterPart,
                    BaseIDBC_Handler *idbcHandler,
                    const Matrix<Double>& elemMat,
                    const StdVector<Integer>& connect1,
                    const StdVector<Integer>& connect2,
                    UInt limit1, UInt limit2,
                    bool setCounterPart);
  template  void BaseEntryManipulator::
  SetElementMatrix( FEMatrixType matrixID,
                    const PdeIdType pdeID1,
                    const PdeIdType pdeID2,
                    StdMatrix *stdMat,
                    StdMatrix *counterPart,
                    BaseIDBC_Handler *idbcHandler,
                    const Matrix<Complex>& elemMat,
                    const StdVector<Integer>& connect1,
                    const StdVector<Integer>& connect2,
                    UInt limit1, UInt limit2,
                    bool setCounterPart);


  // **********************
  //   SetCounterPartOnly
  // **********************
  template<typename T>
  void BaseEntryManipulator::
  SetCounterPartOnly( FEMatrixType matrixID,
                      const PdeIdType pdeID1,
                      const PdeIdType pdeID2,
                      StdMatrix *stdMat,
                      StdMatrix *counterPart,
                      BaseIDBC_Handler *idbcHandler,
                      const Matrix<T>& elemMat,
                      const StdVector<Integer>& connect1,
                      const StdVector<Integer>& connect2,
                      UInt limit1, UInt limit2 ) {


    // Clear the arrays
    rowIndList_.clear();
    colIndList1_.clear();
    colIndList2_.clear();

    UInt aux = 0;
    UInt length1 = connect1.GetSize();
    UInt length2 = connect2.GetSize();

    // STEP 1: Generate row index list from first connect array, dropping
    //         equation numbers for dofs fixed by (in)homogeneous Dirichlet
    //         boundary conditions and changing the sign of those fixed by
    //         constraints.
    for ( UInt i = 0; i < length1; i++ ) {
      aux = static_cast<UInt>(std::abs( connect1[i] ));
      if ( aux <= limit1 && aux > 0 ) {

        // Store equation number
        rowIndList_.push_back( aux - 1 );

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
          colIndList2_.push_back( aux - limit2 - 1 );
          colIndList2_.push_back(  i  );
        }
        else {

          // Store equation number and column index into
          // local element matrix
          colIndList1_.push_back( aux - 1 );
          colIndList1_.push_back(  i  );
        }
      }
    }

    /*
    // output original connectivity
    (*debug) << "\n ------------------------------------------------------"
             << "\n BaseEntryManipulator::SetCounterPartOnly\n"
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
    */


    // STEP 3: Insert values of transpose counterpart into matrix
    UInt rowInd;
    UInt colInd;

    // (*debug) << "\n free <-> free matrix (counterpart):\n";

    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
        colInd = colIndList1_[j+1];

        /*
        (*debug) << " mat(" << colIndList1_[j] << ", "
                 << rowIndList_[i] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
        */         
        counterPart->AddToMatrixEntry( colIndList1_[j], rowIndList_[i],
                                       elemMat[rowInd][colInd] );
      }
    }

    // STEP 4: Insert values for fixed dofs into IDBC_Handler object
    //(*debug) << "\n free <-> fixed matrix:\n";
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList2_.size(); j += 2 ) {
        colInd = colIndList2_[j+1];

        /*(*debug) << " mat(" << rowIndList_[i] << ", "
                 << colIndList2_[j] << ") += "
                 << elemMat[ rowInd * length2 + colInd ] << std::endl;
        */         
        idbcHandler->AddWeightFixedToFree( matrixID, pdeID1, pdeID2,
                                           rowIndList_[i], colIndList2_[j],
                                           elemMat[rowInd][colInd] );
      }
    }

  }

  // explicit template instantiation
  template
  void BaseEntryManipulator::
  SetCounterPartOnly( FEMatrixType matrixID,
                      const PdeIdType pdeID1,
                      const PdeIdType pdeID2,
                      StdMatrix *stdMat,
                      StdMatrix *counterPart,
                      BaseIDBC_Handler *idbcHandler,
                      const Matrix<Double>& elemMat,
                      const StdVector<Integer>& connect1,
                      const StdVector<Integer>& connect2,
                      UInt limit1, UInt limit2 );
  template 
  void BaseEntryManipulator::
  SetCounterPartOnly( FEMatrixType matrixID,
                      const PdeIdType pdeID1,
                      const PdeIdType pdeID2,
                      StdMatrix *stdMat,
                      StdMatrix *counterPart,
                      BaseIDBC_Handler *idbcHandler,
                      const Matrix<Complex>& elemMat,
                      const StdVector<Integer>& connect1,
                      const StdVector<Integer>& connect2,
                      UInt limit1, UInt limit2 );


  // *****************
  //   SetElementRHS
  // *****************
  template<typename T>
  void BaseEntryManipulator::SetElementRHS( SingleVector *rhs,
                                            const Vector<T>& elemRhs,
                                            const StdVector<Integer>& connect,
                                            UInt limit ) {


    /* (*debug) << "SetElementRHS (" << elemSize << "): ";
    for ( UInt i = 0; i < elemSize; i++ ) {
      (*debug) << "(" << connect[i] << "," << elemRHS[i] << ") ";
    }
    (*debug) << std::endl;
    */  
    UInt size = connect.GetSize();
    for ( UInt i = 0; i < size; i++ ) {

      // For constraints: Add values of slaves to master
      if ( connect[i] < 0 ) {
        rhs->AddToEntry( -(connect[i])-1, elemRhs[i] );
      }

      // Insert values related to free degrees of freedom
      // Note: Cast to UInt is save, since connect[i] is
      // not negative
      else if ( connect[i] > 0 && (UInt)connect[i] <= limit ) {
        rhs->AddToEntry( connect[i]-1, elemRhs[i] );
      }
    }

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }
  // explicit template instantiation
  template void BaseEntryManipulator::
  SetElementRHS( SingleVector *rhs, const Vector<Double>& elemRhs,
                 const StdVector<Integer>& connect, UInt limit );
  template void BaseEntryManipulator::
  SetElementRHS( SingleVector *rhs, const Vector<Complex>& elemRhs,
                 const StdVector<Integer>& connect, UInt limit );


  // *****************
  //   GetMatrixEntry 
  // *****************
  template<typename T>
  void BaseEntryManipulator::
  GetMatrixEntry(FEMatrixType matrix_id,
                 const PdeIdType pdeID1,
                 const PdeIdType pdeID2,
                 StdMatrix *stdMat,
                 BaseIDBC_Handler *idbcHandler,
                 T & entry,
                 Integer eqnNr1, 
                 Integer eqnNr2, 
                 UInt limit1, UInt limit2 ) {

    // Calculate absolute value of equation numbers (for constraints)
    UInt aux1 = abs(eqnNr1);
    UInt aux2 = abs(eqnNr2);

    // Check if entry is for a free or a fixed dof
    if ( eqnNr1 <= static_cast<Integer>(limit1) && aux1 > 0 ) {

      // check if entry depends on a free or a fixed dof
      if ( eqnNr2 <= static_cast<Integer>(limit2) && aux2 > 0 ) {
        stdMat->GetMatrixEntry( eqnNr1-1, eqnNr2-1, entry );
      } else {
        idbcHandler->GetWeightFixedToFree( matrix_id, pdeID1, pdeID2,
                                           eqnNr1 - 1, aux2 - limit2 - 1,
                                           entry );
      }
    }else {
      entry = 0.0;
    }
  }
  template void BaseEntryManipulator::
  GetMatrixEntry<Double>(FEMatrixType matrix_id, const PdeIdType pdeID1,  const PdeIdType pdeID2,
                 StdMatrix *stdMat, BaseIDBC_Handler *idbcHandler, Double & entry,
                 Integer eqnNr1, Integer eqnNr2, UInt limit1, UInt limit2 );
  template void BaseEntryManipulator::
  GetMatrixEntry<Complex>(FEMatrixType matrix_id, const PdeIdType pdeID1,  const PdeIdType pdeID2,
                 StdMatrix *stdMat, BaseIDBC_Handler *idbcHandler, Complex& entry,
                 Integer eqnNr1, Integer eqnNr2, UInt limit1, UInt limit2 );


  // *****************
  //   SetMatrixEntry
  // *****************
  template<typename T>
  void BaseEntryManipulator::
  SetMatrixEntry(FEMatrixType matrix_id,
                 const PdeIdType pdeID1,
                 const PdeIdType pdeID2,
                 StdMatrix *stdMat,
                 BaseIDBC_Handler *idbcHandler,
                 const T& entry,
                 Integer eqnNr1, 
                 Integer eqnNr2, 
                 UInt limit1, UInt limit2,
                 bool setCounterPart ) {

    // Calculate absolute value of equation numbers (for constraints)
    UInt aux1 = abs(eqnNr1);
    UInt aux2 = abs(eqnNr2);

    // Check if entry is for a free or a fixed dof
    if ( eqnNr1 <= static_cast<Integer>(limit1) && aux1 > 0 ) {

      // check if entry depends on a free or a fixed dof
      if ( eqnNr2 <= static_cast<Integer>(limit2) && aux2 > 0 ) {
        stdMat->SetMatrixEntry( eqnNr1 - 1, eqnNr2 - 1, entry, setCounterPart );
      } else {
        idbcHandler->SetWeightFixedToFree( matrix_id, pdeID1, pdeID2,
                                           eqnNr1 - 1, aux2 - limit2 - 1,
                                           entry );
      }
    }else {
      EXCEPTION( "Can not set an matrix entry for position eqnNr1 = "
                 << eqnNr1 << ", eqnNr2 = " << eqnNr2 << " for pdeID = "
                 << pdeID1 );
    }

  }
  // explicit template instantiaton
  template void BaseEntryManipulator::
  SetMatrixEntry(FEMatrixType matrix_id, const PdeIdType pdeID1, const PdeIdType pdeID2,
                 StdMatrix *stdMat, BaseIDBC_Handler *idbcHandler,const Double& entry,
                 Integer eqn1, Integer eqnNr2, UInt limit1, UInt limit2,
                 bool setCounterPart ) ;
  template void BaseEntryManipulator::
  SetMatrixEntry(FEMatrixType matrix_id, const PdeIdType pdeID1, const PdeIdType pdeID2,
                 StdMatrix *stdMat, BaseIDBC_Handler *idbcHandler, const Complex& entry,
                 Integer eqn1, Integer eqnNr2, UInt limit1, UInt limit2,
                 bool setCounterPart ) ;
  
  // *************
  //   SetNodeRHS
  // *************
  template<typename T>
  void BaseEntryManipulator::SetNodeRHS( SingleVector *rhs, 
                                         const T& val,
                                         Integer node ) {
  //since the graph will be 0-based, we substract from the
  //equation number always minus one, when accessing the vector
    if ( node != 0 ) {
      node -= 1;
      rhs->AddToEntry( std::abs(node), val );
    }

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }
  
  // explicit template instantiation
  template void BaseEntryManipulator::
  SetNodeRHS( SingleVector *rhs, const Double& val, Integer node );
  template void BaseEntryManipulator::
  SetNodeRHS( SingleVector *rhs, const Complex& val, Integer node );

  // ******************
  //   SetVectorEntry
  // ******************
  template<typename T>
  void BaseEntryManipulator::SetVectorEntry( SingleVector *vec,
                                             UInt index,
                                             const T& newVal) {
    //since the graph will be 0-based, we substract from the
    //equation number always minus one, when accessing the vector
    vec->SetEntry( index - 1, newVal );
  }
  // explicit template instantiation
  template void BaseEntryManipulator::
  SetVectorEntry( SingleVector *vec, UInt index, const Double& newVal);
  template void BaseEntryManipulator::
  SetVectorEntry( SingleVector *vec, UInt index, const Complex& newVal);

  // ***********
  //   InitSol
  // ***********
  template<typename T>
  void BaseEntryManipulator::InitSol( SingleVector *sol,
                                      const Vector<T>& newSol) {


    // For the standard vector class we use fast access
    if ( Vector<T> *mySol = dynamic_cast<Vector<T>*>(sol) ) {
      for ( UInt i = 0; i < sol->GetSize(); i++ ) {
        (*mySol)[i] = newSol[i];
      }
    }

    // otherwise we use the slow access
    // i.e. in the parallel case this doesn't hurt too much
    // since CFS does everything sequentially anyway
    //
    // What makes you think this is slower than the above?
    else {
      for ( UInt i = 0; i < sol->GetSize(); i++ ) {
        sol->SetEntry(i, newSol[i]);
      }
    }

    // Set flag for rhsBuffer invaldiation
    //rhsBufferIsValid_ = false;
  }
  
  // explicit template instantiation
  template void BaseEntryManipulator::
  InitSol( SingleVector *sol, const Vector<Double>& newSol);
  template void BaseEntryManipulator::
  InitSol( SingleVector *sol, const Vector<Complex>& newSol);

  // ***********
  //   InitRHS
  // ***********
  template<typename T>
  void BaseEntryManipulator::InitRHS( SingleVector *rhs,
                                      const Vector<T>& newRHS) {


    // For the standard vector class we use fast access
    if ( Vector<T> *myRHS = dynamic_cast<Vector<T>*>(rhs) ) {
      for ( UInt i = 0; i < rhs->GetSize(); i++ ) {
        (*myRHS)[i] = newRHS[i];
      }
    }

    // otherwise we use the slow access
    // i.e. in the parallel case this doesn't hurt too much
    // since CFS does everything sequentially anyway
    //
    // What makes you think this is slower than the above?
    else {
      for ( UInt i = 0; i < rhs->GetSize(); i++ ) {
        rhs->SetEntry(i, newRHS[i]);
      }
    }
  }
  template void BaseEntryManipulator::
  InitRHS( SingleVector *rhs, const Vector<Double>& newRHS);
  template void BaseEntryManipulator::
  InitRHS( SingleVector *rhs, const Vector<Complex>& newRHS);

  // *******************************
  //   UpdateRHS (f = f + A * fup)
  // *******************************
  template<typename T>
  void BaseEntryManipulator::UpdateRHS( SingleVector *rhs, StdMatrix *stdMat,
                                        const Vector<T>& fup ) {



    // Support for the Vector<Double> case was implemented for
    // HyreMatrix::MultAdd some time ago. Status not quite clear.
    // Commented this out on 21.12.2004 (Marcus)

    // Currently we can only do this with the Vector<T> class,
    // i.e. with OLAS/LAPACK matrices
    // MatrixStorageType matType = sysmat_[matrix_id]->GetStorageType();

    // Perform arithmetic update
    stdMat->MultAdd( fup, *rhs );
  }
  
  // explicit template instantiation
  template void BaseEntryManipulator::
  UpdateRHS( SingleVector *rhs, StdMatrix *stdMat, const Vector<Double>& fup );
  template void BaseEntryManipulator::
  UpdateRHS( SingleVector *rhs, StdMatrix *stdMat, const Vector<Complex>& fup );


  // ======================================================================
  // METHODS FOR HANDLING IDBCs
  // ======================================================================


  // *********************
  //   AdaptSystemMatrix
  // *********************
  void BaseEntryManipulator::AdaptSystemMatrix( StdMatrix &stdMat,
                                                UInt *dirichletEQN,
                                                UInt numIDBC,
                                                Double &penaltyTerm ) {

    // Loop over all Dirichlet equation numbers and replace
    // diagonal matrix entry by penalty term (no components
    // required in this scalar case)
    if( stdMat.GetEntryType() == BaseMatrix::COMPLEX ) { 
      Complex penaltyC(penaltyTerm, 0.0);
      for ( UInt i = 0; i < numIDBC; i++ ) {
        stdMat.SetDiagEntry( dirichletEQN[i]-1, penaltyC );
      } 
    } else {
      for ( UInt i = 0; i < numIDBC; i++ ) {
        stdMat.SetDiagEntry( dirichletEQN[i]-1, penaltyTerm );
      }
    }
  }


  // *******************
  //   AdaptRHSForIDBC
  // *******************
  void BaseEntryManipulator::AdaptRHSForIDBC( SingleVector &rhs,
                                              UInt *dirichletEQN,
                                              SingleVector &dirichletValue,
                                              Double &penaltyTerm,
                                              UInt numIDBC ) {

    if( rhs.GetEntryType() == BaseMatrix::COMPLEX ) { 
      Complex penaltyC(penaltyTerm, 0.0);
      Complex entry;

      for( UInt i = 0; i < numIDBC; i++ ) {
        dirichletValue.GetEntry( i, entry );
        entry *= penaltyC;
        rhs.SetEntry( dirichletEQN[i] - 1, entry );
      }
    } else {
      Double entry;

      for( UInt i = 0; i < numIDBC; i++ ) {
        dirichletValue.GetEntry( i, entry );
        entry *= penaltyTerm;
        rhs.SetEntry( dirichletEQN[i] - 1, entry );
      }
    }
  }


  

  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  template<typename T>
  void BaseEntryManipulator::AddToDiagMatrixEntry( StdMatrix *stdMat,
                                                   Integer eqnNum,
                                                   const T& val ) {

    // Get hold of matrix entry and adapt it
    Double entry;
    stdMat->GetDiagEntry( std::abs(eqnNum)-1, entry );
    entry += val;
    stdMat->SetDiagEntry( std::abs(eqnNum)-1, entry );
  }
  
  
  // instantiate template methods
  

}
