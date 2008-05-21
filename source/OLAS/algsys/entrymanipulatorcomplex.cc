// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <complex>

#include "algsys/entrymanipulatorcomplex.hh"
#include "algsys/baseidbchandler.hh"
#include "algsys/standardsys.hh"
#include "matvec/matvec.hh"
#include "matvec/vector.hh"

#include "DataInOut/Logging/cfslog.hh"
#include "Utils/StdVector.hh"
#include "General/exception.hh"

using CoupledField::Exception;
using CoupledField::StdVector;

DECLARE_LOG(entryManip)
DEFINE_LOG(entryManip, "entryManip")

namespace OLAS {


  // ********************
  //   SetElementMatrix
  // ********************
  void EntryManipulatorComplex::
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
                    bool setCounterPart) {


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
    
    // output original connectivity
    LOG_DBG(entryManip) << "limit1=" << limit1 << " limit2 << " << limit2;
    LOG_DBG3(entryManip) << "connect1=" << StdVector<int>::ToString(length1, connect1);
    LOG_DBG3(entryManip) << "connect2=" << StdVector<int>::ToString(length2, connect2);
    // LOG_DBG3(entryManip) << "rowIndList_=" << StdVector<unsigned int>::ToString(rowIndList_.size(), rowIndList_[0], 0, 2);
    LOG_DBG3(entryManip) << "colIndList1_=" << StdVector<unsigned int>::ToString(colIndList1_.size(), &colIndList1_[0], 0, 2);
    LOG_DBG3(entryManip) << "colIndList2_=" << StdVector<unsigned int>::ToString(colIndList2_.size(), &colIndList2_[0], 0, 2);

    // STEP 3a: Insert values into matrix for complex dofs
    UInt rowInd;
    UInt colInd;
    Double real;
    Double imag;

    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
        colInd = colIndList1_[j+1];
        real = elemMat[ rowInd * length2 + colInd ];
        imag = elemMat[ rowInd * length2 + colInd + length1 * length2 ];
        Complex cVal( real, imag );

        LOG_DBG3(entryManip) << "mat[" << rowIndList_[i] << "][" << colIndList1_[j] 
                             << "] += (" << cVal;
        stdMat->AddToMatrixEntry( rowIndList_[i], colIndList1_[j], cVal );

      }
    }

    // STEP 3b: Insert values of transposed counter part into matrix
    if ( setCounterPart == true ) {
      for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
        rowInd = rowIndList_[i+1];
        for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
          colInd = colIndList1_[j+1];
          real = elemMat[ rowInd * length2 + colInd ];
          imag = elemMat[ rowInd * length2 + colInd + length1 * length2 ];
          Complex cVal( real, imag );

          LOG_DBG2(entryManip) << "counterpart: mat[" << colIndList1_[j] << "][" << colIndList1_[i] 
                               << "] += (" << cVal;

          counterPart->AddToMatrixEntry( colIndList1_[j], rowIndList_[i],
                                         cVal );
        }
      }
    }

    // STEP 4: Insert values for fixed dofs into IDBC_Handler object
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList2_.size(); j += 2 ) {
        colInd = colIndList2_[j+1];

        LOG_DBG2(entryManip) << "idbcHandler->SetWeightFixedToFree mat[" << rowIndList_[i]
                             << ", " << colIndList2_[j] << "] += "
                             << elemMat[ rowInd * length2 + colInd ] << ", "
                             << elemMat[ rowInd * length2 + colInd + length1 * length2];

        idbcHandler->SetWeightFixedToFree( matrixID, pdeID1, pdeID2, rowIndList_[i], colIndList2_[j],
                                           elemMat[ rowInd * length2 + colInd ],
                                           elemMat[ rowInd * length2 + colInd  + length1 * length2 ]);
      }
    }
  }


  // **********************
  //   SetCounterPartOnly
  // **********************
  void EntryManipulatorComplex::
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
             << "\n EntryManipulatorComplex::SetCounterPartOnly\n"
             << " limit1 = " << limit1 << '\n'
             << " limit2 = " << limit2 << '\n';
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

    // STEP 3: Insert values of transposed counter part into matrix
    UInt rowInd;
    UInt colInd;
    Double real;
    Double imag;

#ifdef DEBUG_ASSEMBLE
    (*debug) << "\n free <-> free matrix (counterpart):\n";
#endif
    for ( UInt i = 0; i < rowIndList_.size(); i += 2 ) {
      rowInd = rowIndList_[i+1];
      for ( UInt j = 0; j < colIndList1_.size(); j += 2 ) {
        colInd = colIndList1_[j+1];
        real = elemMat[ rowInd * length2 + colInd ];
        imag = elemMat[ rowInd * length2 + colInd + length1 * length2 ];
        Complex cVal( real, imag );

#ifdef DEBUG_ASSEMBLE
        (*debug) << "mat(" << colIndList1_[j] << ", "
                 << rowIndList_[j] << ") += " << cVal << std::endl;
#endif

        counterPart->AddToMatrixEntry( colIndList1_[j], rowIndList_[i],
                                       cVal );
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
                 << elemMat[ rowInd * length2 + colInd ]
                 << ", "
                 << elemMat[ rowInd * length2 + colInd + length1 * length2]
                 << std::endl;
#endif
        idbcHandler->SetWeightFixedToFree( matrixID, pdeID1, pdeID2,
                                           rowIndList_[i], colIndList2_[j],
                                           elemMat[ rowInd * length2
                                                    + colInd ],
                                           elemMat[ rowInd * length2
                                                    + colInd
                                                    + length1 * length2 ]);
      }
    }
  }


  // *****************
  //   SetElementRHS
  // *****************
  void EntryManipulatorComplex::SetElementRHS( SparseVector *rhs,
                                               Double *elemRHS,
                                               Integer *connect, 
                                               UInt elemSize,
                                               UInt limit) {


    for ( UInt i = 0; i < elemSize; i++ ) {

      // Insert values related to free degrees of freedom
      if ( connect[i] > 0 && connect[i] <= limit ) {
        Complex val( elemRHS[i], elemRHS[ i + elemSize ] );
        rhs->AddToVectorEntry( connect[i], val );
      }

      // For constraints: Add values of slaves to master
      else if ( connect[i] < 0 ) {
        Complex val( elemRHS[i], elemRHS[ i + elemSize ] );
        rhs->AddToVectorEntry( -connect[i], val );
      }
    }

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }

  // ***************************
  //   GetMatrixEntry (Double)
  // ***************************
  void EntryManipulatorComplex::
  GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Double & entry, Integer eqnNr1,
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2 ) {
    (*error) << "EntryManipulatorComplex::GetMatrixEntry: "
             << "Interface for double not implemented, since this "
             << "makes no sense!";
    Error( __FILE__, __LINE__ );
  }

  // ***************************
  //   GetMatrixEntry (Complex)
  // ***************************
  void EntryManipulatorComplex::
  GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Complex & entry, Integer eqnNr1,
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2 ) {

    Double real, imag = 0.0;

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
                                           real, imag );
        entry = Complex( real, imag );
      } 
    }else {
      entry = Complex( 0.0, 0.0 );
    }

  }

  // ***************************
  //   SetMatrixEntry (Double)
  // ***************************
  void EntryManipulatorComplex::
  SetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Double  entry, Integer eqnNr1,
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2, bool setCounterPart ) {
    (*error) << "EntryManipulatorComplex::SetMatrixEntry: "
             << "Interface for double not implemented, since this "
             << "makes no sense!";
    Error( __FILE__, __LINE__ );
  }

  // ***************************
  //   SetMatrixEntry (Complex)
  // ***************************
  void EntryManipulatorComplex::
  SetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                  const PdeIdType pdeID2, StdMatrix *stdMat,
                  BaseIDBC_Handler *idbcHandler,
                  Complex entry, Integer eqnNr1, 
                  Integer eqnNr2, 
                  UInt limit1, UInt limit2, bool setCounterPart ) {

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
                                           entry.real(), entry.imag() );
      } 
    }else {
      (*error) << "Can not set an matrix entry for position eqnNr1 = "
               << eqnNr1 << ", eqnNr2 = " << eqnNr2 << " for pdeID = "
               << pdeID1 << "!\n";
      Error( __FILE__, __LINE__ );
    }

  }

  // ***********************
  //   SetNodeRHS (Double)
  // ***********************
  void EntryManipulatorComplex::SetNodeRHS( SparseVector *rhs, Double val,
                                           Integer node  ) {
    
    if ( node != 0 ) {
      Complex aux = Complex(val, 0.0);
      rhs->AddToVectorEntry( std::abs(node), val );
      
      // Set flag for rhsBuffer invaldiation
      rhsBufferIsValid_ = false;
    }
  }

  // **************
  //   SetNodeRHS
  // **************
  void EntryManipulatorComplex::SetNodeRHS( SparseVector *rhs, Complex val,
                                            Integer node ) {


    if ( node != 0 ) {
      rhs->AddToVectorEntry( std::abs(node), val );
      
      // Set flag for rhsBuffer invaldiation
      rhsBufferIsValid_ = false;
    }
  }


  // ******************
  //   SetVectorEntry
  // ******************
  void EntryManipulatorComplex::SetVectorEntry( SparseVector *vec, UInt index,
                                                Double &newVal ) {
    (*error) << "EntryManipulatorComplex::SetVectorEntry: "
             << "Interface for Double not implemented!";
    Error( __FILE__, __LINE__ );
  }


  void EntryManipulatorComplex::SetVectorEntry( SparseVector *vec, UInt index,
                                                Complex &newVal ) {
    vec->SetVectorEntry( index, newVal );
  }

  
  // ***********
  //   InitRHS
  // ***********
  void EntryManipulatorComplex::InitSol( SparseVector *sol,
                                         const Double *newSol ){

	    UInt size = sol->GetSize();

	    Complex aux = Complex(0.0, 0.0);

	    Warning ( "EntryManipulatorComplex::InitRHS() has not been tested yet",
	              __FILE__, __LINE__ );

	    // For the standard vector class we use fast access
	    if ( Vector<Complex> *mySol = dynamic_cast<Vector<Complex>*>(sol) ) {
	      for ( Integer i = 1; i <= size; i++ ) {
	        aux = Complex( (newSol-1)[i], (newSol-1)[i+size] );
	        (*mySol)[i] = aux;
	      }
	    }

	    // otherwise we use the slow access
	    // i.e. in the parallel case this doesn't hurt too much
	    // since CFS does everything sequentially anyway
	    //
	    // What makes you think this is slower than the above?
	    else {
	      for ( Integer i = 1; i <= size; i++ ) {
	        aux = Complex( (newSol-1)[i], (newSol-1)[i+size] );
	        sol->SetVectorEntry(i, aux);
	      }
	    }

	    // Set flag for rhsBuffer invaldiation
	    //rhsBufferIsValid_ = false;
  }
  
  // ***********
  //   InitRHS
  // ***********
  void EntryManipulatorComplex::InitRHS( SparseVector *rhs,
                                         const Double *newRHS ) {


    UInt size = rhs->GetSize();

    Complex aux = Complex(0.0, 0.0);

    // For the standard vector class we use fast access
    if ( Vector<Complex> *myRHS = dynamic_cast<Vector<Complex>*>(rhs) ) {
      for ( unsigned int i = 1; i <= size; i++ ) {
        aux = Complex( (newRHS-1)[i], (newRHS-1)[i+size] );
        (*myRHS)[i] = aux;
      }
    }

    // otherwise we use the slow access
    // i.e. in the parallel case this doesn't hurt too much
    // since CFS does everything sequentially anyway
    //
    // What makes you think this is slower than the above?
    else {
      for (unsigned int i = 1; i <= size; i++ ) {
        aux = Complex( (newRHS-1)[i], (newRHS-1)[i+size] );
        rhs->SetVectorEntry(i, aux);
      }
    }

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }


  // ***********
  //   InitRHS
  // ***********
  void EntryManipulatorComplex::InitRHS( SparseVector *rhs, const Vector<Complex>* newRHS ) 
  {
    assert(rhs->GetSize() == newRHS->GetSize());
    
    // we assuse our rhs is complex
    Vector<Complex> *myRHS = dynamic_cast<Vector<Complex>*>(rhs);
    assert(myRHS != NULL);

    for ( unsigned int i = 1; i <= rhs->GetSize(); i++ ) 
      (*myRHS)[i] = (*newRHS)[i];

    // Set flag for rhsBuffer invaldiation
    rhsBufferIsValid_ = false;
  }
  
  
  // *******************************
  //   UpdateRHS (f = f + A * fup)
  // *******************************
  void EntryManipulatorComplex::UpdateRHS( SparseVector *rhs,
                                           StdMatrix *stdMat, Double *fup ) { 
    (*error) << "EntryManipulatorComplex::UpdateRHS not implemented yet!";
    Error( __FILE__, __LINE__ );
  }


  // ======================================================================
  // METHODS FOR HANDLING IDBCs
  // ======================================================================


  // *********************
  //   AdaptSystemMatrix
  // *********************
  void EntryManipulatorComplex::AdaptSystemMatrix( StdMatrix &stdMat,
                                                   UInt *dirichletEQN,
                                                   UInt numIDBC,
                                                   Double &penaltyTerm ) {


    // Loop over all Dirichlet equation numbers and replace
    // diagonal matrix entry by penalty term (no components
    // required in this scalar case)
    Complex aux( penaltyTerm, 0.0 );
    for ( UInt i = 1; i <= numIDBC; i++ ) {
      stdMat.SetDiagEntry( dirichletEQN[i], aux );
    }
  }


  // *******************
  //   AdaptRHSForIDBC
  // *******************
  void EntryManipulatorComplex::AdaptRHSForIDBC( SparseVector &rhs,
                                                 UInt *dirichletEQN,
                                                 SparseVector &dirichletValue,
                                                 Double &penaltyTerm,
                                                 UInt numIDBC ) {


    Complex entry;

    for( UInt i = 1; i <= numIDBC; i++ ) {
      dirichletValue.GetVectorEntry( i, entry );
      entry *= penaltyTerm;
      rhs.SetVectorEntry( dirichletEQN[i], entry );
    }
  }


  // ******************
  //   GetSolutionVal
  // ******************
  void EntryManipulatorComplex::GetSolutionVal( SparseVector *sol,
                                                Double* &ptSol, 
                                                const PdeIdType pdeID ) {
    (*error) << "GetSolutionVal(Double) not defined for complex entries";
    Error( __FILE__, __LINE__ );
  }


  void EntryManipulatorComplex::GetSolutionVal( SparseVector *sol,
                                                Complex* &ptSol, 
                                                const PdeIdType pdeID ) {


    // only update the solution buffer if it is invalid
    if ( solBufferIsValid_ == false ) {
      // copy complex data in return buffer of Double-type
      solBuffer_ = *sol;
      
      // Set invalidation flag
      solBufferIsValid_ = true;
    }
      
    // Adapt 1-based array to 0-based array
    // for CFS++
    ptSol = solBuffer_.GetPointer();
    ptSol++;
  }


  // *************
  //   GetRHSVal
  // *************
  void EntryManipulatorComplex::GetRHSVal( SparseVector *rhs, Double* &ptRhs, 
                                           const PdeIdType identifierPDE ) 
  {
    EXCEPTION("GetRHSVal(Double) not defined for complex entries");
  }

  void EntryManipulatorComplex::GetRHSVal( SparseVector *rhs, Complex* &ptRhs, 
                                           const PdeIdType identifierPDE ) {


    // only update the rhs buffer if it is invalid
    if ( rhsBufferIsValid_ == false ) {

#ifdef DEBUG_ASSEMBLE
      (*debug) << "GetRHSVal: case sequential OLAS" << std::endl;
#endif
      // copy vector into buffer
      rhsBuffer_ = *rhs;
      
      // Set invalidation flag
      rhsBufferIsValid_ = true;
    }
      
    // Adapt 1-based array to 0-based array
    // for CFS++
    ptRhs = rhsBuffer_.GetPointer();
    ptRhs++;
  }


  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  void EntryManipulatorComplex::AddToDiagMatrixEntry( StdMatrix *stdMat, 
                                                      Integer eqnNum,
                                                      Double *val ) {


    // Perform consistency checks
#ifdef DEBUG_ASSEMBLE
    if ( dof != 1 ) {
      (*error) << "EntryManipulatorComplex::AddToDiagMatrixEntry: "
               << "Parameter dof = "
               << dof << ", but should be 1";
      Error( __FILE__, __LINE__ );
    }
    if ( val == NULL ) {
      (*error) << "EntryManipulatorComplex::AddToDiagMatrixEntry: "
               << "Parameter val is a NULL pointer";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Get hold of matrix entry and adapt it
    Complex entry;
    Complex cVal( val[0], val[1] );
    stdMat->GetDiagEntry( std::abs(eqnNum), entry );
    entry += cVal;
    stdMat->SetDiagEntry( std::abs(eqnNum), entry );
  }

}
