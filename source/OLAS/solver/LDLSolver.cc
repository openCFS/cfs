#include <iterator>
#include <list>

#include "MatVec/opdefs.hh"
#include "MatVec/SCRS_Matrix.hh"

#include "LDLSolver.hh"

#ifdef DEBUG_LDLSOLVER
#define DEBUG_LDLSOLVER_ANALYSE
#define DEBUG_LDLSOLVER_FACTORISE
#endif

// temporary (to be deleted once class is finished)
// #define DEBUG_LDLSOLVER_ANALYSE
// #define DEBUG_LDLSOLVER_FACTORISE


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  LDLSolver<T>::LDLSolver( PtrParamNode solverNode, PtrParamNode olasInfo ) {

    xml_ = solverNode;
    infoNode_ = olasInfo->Get("directLDL");

    // No factorisation was performed yet
    amFactorised_ = false;

    // Initialise all pointers to zero
    denseVec_   = NULL;
    firstU_     = NULL;
    scanList_   = NULL;
    activeList_ = NULL;
    dataD_      = NULL;
    dataU_      = NULL;
    cidxU_      = NULL;
    rptrU_      = NULL;

    // No problem size known yet
    sysMatDim_  = 0;
    auxVecSize_ = 0;

    itRefSteps_ = 2;
    xml_->GetValue("itRefSteps", itRefSteps_, ParamNode::INSERT);

    // NOTE: The class currently gives wrong results, if used with
    //       block entries instead of scalar entries. Thus, we do
    //       a check here to avoid problems with production runs.
    if ( BlockSize<T>::size != 1 ) {
      EXCEPTION( "LDLSolver class currently gives wrong results for "
               << "matrices with non-scalar entries! Refusing to "
               << "generate a solver for matrix with block size = "
               << BlockSize<T>::size << "." );
    }

  }


  // **************
  //   Destructor
  // **************
  template<typename T>
  LDLSolver<T>::~LDLSolver() {


    // Free dynamically allocated memory
    delete [] ( denseVec_ );  denseVec_  = NULL;
    delete [] ( firstU_ );  firstU_  = NULL;
    delete [] ( scanList_ ); scanList_ = NULL;
    delete [] ( activeList_ ); activeList_ = NULL;
    delete [] ( dataD_ );  dataD_  = NULL;
    delete [] ( dataU_ );  dataU_  = NULL;
    delete [] ( cidxU_ );  cidxU_  = NULL;
    delete [] ( rptrU_ );  rptrU_  = NULL;

  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void LDLSolver<T>::Setup( BaseMatrix &sysMat ) {


    // Check that we have a StdMatrix
    if ( sysMat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "LDLSolver cannot deal with matrices other than "
               << "SCRS_Matrix/sparseSym" );
    }

    StdMatrix& stdMat = dynamic_cast<StdMatrix&>(sysMat);

      // Now test the storage layout
      BaseMatrix::StorageType sType = stdMat.GetStorageType();
      if ( sType != BaseMatrix::SPARSE_SYM ) {
        EXCEPTION( "LDLSolver::Setup: The LDLSolver requires the system "
                 << "matrix to be an SCRS_Matrix, i.e. sparseSym. The system "
                 << "matrix you supplied is a matrix in '"
                 << MatrixStorageTypeEnum.ToString( sType )
                 << "' format." );
      }

      // Down-cast to SCRS_Matrix
      SCRS_Matrix<T>& scrsMat = dynamic_cast<SCRS_Matrix<T>&>(sysMat);

      // Get new problem size and perform consistency check
      if ( amFactorised_ == false ) {
        sysMatDim_ = scrsMat.GetNumCols();
      }
      else {
        // TODO: newMatrixPattern was previously only set to false in olasparams.cc
        // and never changed. Is it really necessary? Investigate further!
        bool newMatrixPattern = false;
        if ( !newMatrixPattern &&
             sysMatDim_ != scrsMat.GetNumCols() ) {
          EXCEPTION( "LDLSolver::Setup: newMatrixPattern = false, but "
                   << "matrix dimension changed from " << sysMatDim_ << " to "
                   << scrsMat.GetNumCols() );
        }
        else {
          sysMatDim_ = scrsMat.GetNumCols();
        }
      }

      // Check dimension of problem
      if ( sysMatDim_ <= 2 ) {
        EXCEPTION( "Sorry but your problem is too small for the LDLSolver! "
                 << "We only solve problems with n > 2!");
      }

      // TODO: newMatrixPattern was previously only set to false in olasparams.cc
      // and never changed. Is it really necessary? Investigate further!
      bool newMatrixPattern = false;

      // If this is the first time we are asked for a factorisation, or
      // if there is a new matrix pattern, we start an analyse phase
      if ( amFactorised_ == false ||
           newMatrixPattern ) {

        // Clear "old" vectors
        delete [] ( dataD_ );  dataD_  = NULL;
        delete [] ( dataU_ );  dataU_  = NULL;
        delete [] ( cidxU_ );  cidxU_  = NULL;
        delete [] ( rptrU_ );  rptrU_  = NULL;

        // Perform analysis
        Analyse( scrsMat );

        // If problem size has increased free old dense auxilliary vectors
        if ( auxVecSize_ < scrsMat.GetNumCols() ) {
          delete [] ( denseVec_ );  denseVec_  = NULL;
          delete [] ( firstU_ );  firstU_  = NULL;
	        delete [] ( scanList_ ); scanList_ = NULL;
	        delete [] ( activeList_ ); activeList_ = NULL;
          denseVec_   = NULL;
          firstU_     = NULL;
	        scanList_   = NULL;
	        activeList_ = NULL;
	        auxVecSize_ = 0;
        }

        // Allocate memory for dense auxilliary vectors and lists
        // (and initialise with zeros)
        if ( denseVec_ == NULL ) {
          auxVecSize_ = scrsMat.GetNumCols();
          NEWARRAY( denseVec_, T, auxVecSize_ );
          NEWARRAY( firstU_, UInt, auxVecSize_ );
	        scanList_   = new Integer[auxVecSize_ + 1];
	        activeList_ = new Integer[auxVecSize_ + 1];

          for ( UInt j = 0; j < auxVecSize_; j++ ) 
            denseVec_[j]   = 0.0;

          for ( UInt j = 0; j <= auxVecSize_; j++ ) {
            scanList_[j]   = -1;
            activeList_[j] = -1;
          }
	        scanList_[0]   = auxVecSize_+1;
	        activeList_[0] = auxVecSize_+1;
        }
      }

      // Now trigger factorisation
      Factorise( scrsMat );
      amFactorised_ = true;

      // Convert D to D^{-1}
      for ( UInt i = 0; i < sysMatDim_; i++ ) {
        dataD_[i] = OpType<T>::invert( dataD_[i] );
      }

      // ===============
      //  Export Factor
      // ===============

      // if the user desires it, we will export the matrix factor
      // to a file in Matrix Market format
      bool saveFacToFile = false;
      std::string facFileName = "fac.out";
      bool savePatternOnly = false;
      
      if(xml_->Has("saveFacFile")) {
        saveFacToFile = true;
        xml_->GetValue("saveFacFile", facFileName, ParamNode::INSERT);
      }

      xml_->GetValue("savePatternOnly", savePatternOnly, ParamNode::INSERT);
      
      if( saveFacToFile &&  (!savePatternOnly) ) {
        ExportFactorisation( facFileName.c_str(), false );
      }
  }


  // *********
  //   Solve
  // *********
  template<typename T>
  void LDLSolver<T>::Solve( const BaseMatrix  &sysMat,
                            const BaseVector  &rhs,
                            BaseVector &sol ) {

    // Test that a factorisation is available, if not issue an error
    if ( amFactorised_ == false ) {
      EXCEPTION( "LDLSolver::Solve: No factorisation available. "
	       << "Call Setup() first!" );
    }

    // Solve the problem
    const Vector<T>& myRHS = dynamic_cast<const Vector<T>&>(rhs);
    Vector<T>& mySol = dynamic_cast<Vector<T>&>(sol);
    
      this->SolveLDLSystem( &(cidxU_[0]), &(rptrU_[0]), &(dataU_[0]),
                      &(dataD_[0]), mySol, myRHS, sysMatDim_ );


      // If desired perform iterative refinement
      UInt logLevel = 2;
      UInt numSteps = itRefSteps_;
      
      xml_->GetValue("itRefVerbosity", logLevel, ParamNode::INSERT);

      if ( itRefSteps_ > 0 ) {
        // Avoid recursion, we do not want to do refinement on the
        // solution in the refinement step
        itRefSteps_ = 0;

        // Refine
        iterativeRefiner_.Refine( (*this), sysMat, sol, rhs, numSteps,
                                  logLevel );

        // Reset number of iterations after refinement
        itRefSteps_ = numSteps;
      }

  }


  // ***********
  //   Analyse
  // ***********
  template<typename T>
  void LDLSolver<T>::Analyse( const SCRS_Matrix<T> &sysMat ) {

		// COMPWARNING: unused variable UInt nnzInRowOfU
    UInt i, j, k, nnzInRowOfA, colInd, succ, pred, lastNZ, tmp;

    // Get hold of column index array
    const UInt *cidxA = sysMat.GetColPointer();

    // Get hold of row pointer index array
    const UInt *rptrA = sysMat.GetRowPointer();

    // Query number of matrix rows and columns
    UInt nRows = sysMat.GetNumRows();
    UInt nCols = sysMat.GetNumCols();

    // Allocate auxilliary index array
    Integer *auxVec = new Integer[nCols+1];
    ASSERTMEM( auxVec, sizeof( Integer ) * nCols );
    UInt auxVecNumEntries = 0;

#ifdef DEBUG_LDLSOLVER_ANALYSE
    (*debug) << " LDLSOLVER:\n Phase: ANALYSE\n" << std::endl;
#endif


    // ========================
    //  Determine the profile
    // ========================

    // We use a bottom up approach in order to determine for each matrix
    // column the smallest row index of a non-zero entry in that column.
    i = nRows;
    while ( i > 0 ) {
      i--;
      //    for ( i = nRows; i > 0; i-- ) {
      for ( j = (UInt)rptrA[i]; j < (UInt)rptrA[i+1]; j++ ) {
        auxVec[ cidxA[j] ] = i;
      }
    }

#ifdef DEBUG_LDLSOLVER_ANALYSE
    (*debug) << " column | smallest row index" << std::endl;
    for ( i = 0; i < nCols; i++ ) {
      (*debug) << " " << i << " | " << auxVec[i] << std::endl;
    }

    // Maximal fill-in will generate a column with all non-zero entries
    // between the diagonal and the row determined above. Summing these
    // differences up gives the profile
    UInt profile = 1;
    for ( i = 1; i < nCols; i++ ) {
      profile += i - auxVec[i] + 1;
    }

    (*debug) << " column | max number of entries in factor L" << std::endl;
    for ( i = 0; i < nCols; i++ ) {
      (*debug) << " " << i << " | " << i - auxVec[i] + 1 << std::endl;
    }
    (*debug) << " Factor (L + D) will contain " << profile << " entries"
             << " (at most)" << std::endl;
#endif

    // ============================================
    //  Memory Allocation (Phase 1):
    // ============================================

    // We allocate memory for those vectors of the
    // factorisation whose length is already fixed
    NEWARRAY( dataD_, T   , nCols     );
    NEWARRAY( rptrU_, UInt, nRows + 1 );

    // We allocate memory for the auxilliary 2D data structure that we use
    // to dynamically build the column index information.
    //
    // cIndex[j] is an array with admissible indices
    // 1 ... number of non-zeros in j-th row of D
    // containing the column indices of these non-zero entries
    UInt **cIndex = NULL;
    UInt *auxPtr = NULL;
    NEWARRAY( cIndex, UInt*, nRows );

#ifdef DEBUG_LDLSOLVER_ANALYSE
    (*debug) << " Allocated memory (phase 1)"
             << std::endl;
#endif

    // =========================================
    //  Determine fill-pattern of factor U=L^T
    // =========================================

    // We need one auxilliary dense index vector to store a linked list which
    // contains the column indices of the fill-pattern of row k. We can re-use
    // auxVec for this purpose, but must zero it now. We store the header of
    // the linked list, i.e. the index of the first element in auxVec[0] and
    // use the value nCols + 1 as terminator
    for ( i = 0; i <= nCols; i++ ) {
      auxVec[i] = -1;
    }
    auxVecNumEntries = 0;

    // Initialise row pointer array
    rptrU_[0] = 0;

    // We keep in an STL list the indices of all rows that have not yet
    // contributed their fill-patterns to rows with larger index and must
    // thus be checked for non-zero entry in column k, when determining
    // the fill-pattern of row k. The list entry encodes the index of the row
    std::list<UInt> rowList;
    std::list<UInt>::iterator listIt;

    // If the first row has not only a diagonal entry (unlikely of course in
    // FEM) Initially only the first row is in the row list.
    if ( rptrA[1] - rptrA[0] > 1 ) {
      rowList.push_back( 0 );
    }

    // Pattern of first row is identical to that of first row of U = L^T
    // minus the diagonal entry
    NEWARRAY( cIndex[0], UInt, (rptrA[1]-1) );
    auxPtr = cIndex[0];
    for ( i = 1; i < (UInt)rptrA[1]; i++ ) {
      auxPtr[i-1] = cidxA[i];
    }
    rptrU_[1] = (UInt)rptrA[1] - 1;

    // Needed for writing progress report of factorisation
    UInt percentDone = 0;
    Double actDone = 0.0;

    // ----------------------------------------
    //  main loop for determining fill-pattern
    // ----------------------------------------

    // The last row of L does not contain off-diagonal entries and the first
    // row was already treated, so we only go from 2 <-> n - 1
    for ( k = 1; k < nRows-1; k++ ) {

      // Keep user informed on progress
      actDone = (double)(k*100) / (double)nRows;
      actDone = (UInt)(actDone/10.0)*10;
      if ( actDone > percentDone ) {
        percentDone = (UInt)actDone;
      }

      // Determine number of non-zero entries in this row
      nnzInRowOfA = rptrA[k+1] - rptrA[k];

      // Insert non-zero pattern of row k of A into dense auxilliary vector
      // (omitting the diagonal entry) as linked list
      if ( nnzInRowOfA > 1 ) {
        auxVec[ 0 ] = cidxA[ rptrA[k]+1 ];
        for ( i = (UInt)rptrA[k] + 2; i < (UInt)rptrA[k+1]; i++ ) {
          auxVec[ cidxA[i-1] ] = cidxA[i];
        }
        auxVec[ cidxA[i-1] ] = nCols+1;
      }
      auxVecNumEntries = nnzInRowOfA - 1;

#ifdef DEBUG_LDLSOLVER_ANALYSE
      (*debug) << " Pattern of row " << k << " of A:";
      if ( auxVec[0] != 0 ) {
        succ = auxVec[0];
        while( succ != nCols + 1 ) {
          (*debug) << " " << succ;
          succ = auxVec[succ];
        }
      }
      (*debug) << std::endl;
      listIt = rowList.begin();
      (*debug) << " --> Row list contains row(s):";
      while ( listIt != rowList.end() ) {
        (*debug) << " " << (*listIt);
        listIt++;
      }
      (*debug) << std::endl;
#endif

      // For all rows in the current row list, check whether the first
      // off-diagonal entry has index k, so that it contributes to the
      // fill-pattern for this row
      listIt = rowList.begin();
      while ( listIt != rowList.end() ) {

        // Check that row has non-zero entry in position k
        // (if row has not yet contributed to pattern, then
        // maybe its first off-diagonal entry is in column k)
        if ( rptrU_[(*listIt)+1] - rptrU_[(*listIt)] > 1 &&
             cIndex[(*listIt)][0] == k ) {

#ifdef DEBUG_LDLSOLVER_ANALYSE
          (*debug) << " --> Adding pattern of row " << (*listIt)
                   << " to pattern of row " << k << std::endl;
#endif

          // Initialise variable storing index of last non-zero entry
          // in pattern / linked list
          lastNZ = 0;

          // Loop over pattern of this row adding new entries
          // at correct positions in linked list (omitting diagonal entry)
          auxPtr = cIndex[(*listIt)];
          for ( j = 1; j < rptrU_[(*listIt)+1] - rptrU_[(*listIt)]; j++ ) {

            // Get corresponding column index
            colInd = auxPtr[j];
            // (*debug) << (*listIt) << ": colInd = " << colInd << std::endl;

            // If in pattern update lastNZ index
            if ( auxVec[ colInd ] != -1 ) {
              lastNZ = colInd;
            }

            // Not yet in pattern, so insert it into list
            else {

              // First we have to find the correct lastNZ
              // by walking through the linked list
              while( auxVec[ lastNZ ] < static_cast<Integer>(colInd) ) {
                lastNZ = auxVec[ lastNZ ];
              }

              // Now insert fill-in into the list
              tmp = auxVec[ lastNZ ];
              auxVec[ lastNZ ] = colInd;
              auxVec[ colInd ] = tmp;

              // And update list length
              auxVecNumEntries++;

            }
          }

          // Now we can delete this row from the row list
          // (increments iterator at the same time)
          listIt = rowList.erase( listIt );

        }

        // Nothing to do for this row, so proceed to next row
        else {
          listIt++;
        }
      }

      // The current row must now be added to the list
      if ( auxVec[0] < static_cast<Integer>(nCols + 1) ) {
        rowList.push_back( k );
      }


      // Report fill-pattern of row in factor L^T = U
#ifdef DEBUG_LDLSOLVER_ANALYSE
      (*debug) << " Pattern of row " << k << " of L^T = U:";
      if ( auxVec[0] != 0 ) {
        succ = auxVec[0];
        while( succ != nCols + 1 ) {
          (*debug) << " " << succ;
          succ = auxVec[succ];
        }
      }
      (*debug) << '\n' << std::endl;
#endif

      // Allocate space for column indices for this row
      if ( auxVec[0] != -1 ) {
        NEWARRAY( cIndex[k], UInt, auxVecNumEntries );
      }
      else {
        cIndex[k] = NULL;
      }

      // Now we must insert the column index information into the
      // auxilliary data structure
      if ( auxVec[0] != -1 ) {
        auxPtr = cIndex[k];
        succ = auxVec[0];
        auxVec[0] = nCols + 1;
        i = 0;
        while( succ != nCols + 1 ) {
          auxPtr[i] = succ;
          pred = succ;
          succ = auxVec[succ];
          auxVec[pred] = -1;
          i++;
        }
        rptrU_[k+1] = rptrU_[k] + auxVecNumEntries;
        //(*debug) << "rptrU_[" << k << "] = " << rptrU_[k] << std::endl;
        //(*debug) << "rptrU_[" << k+1 << "] = " << rptrU_[k+1] << std::endl;
      }
      else {
        rptrU_[k+1] = rptrU_[k];
      }
    }

    // The upper loop only went to n - 1, so we finalise the row pointer
    // array explicitely
    rptrU_[sysMatDim_] = rptrU_[sysMatDim_-1];

    // ============================================
    //  Memory Allocation (Phase 2):
    // ============================================

    // Now we know how many entries the pattern will have, so we can
    // do the final allocation
    NEWARRAY( cidxU_, UInt, rptrU_[sysMatDim_] );
    NEWARRAY( dataU_, T   , rptrU_[sysMatDim_] );


#ifdef DEBUG_LDLSOLVER_ANALYSE
    (*debug) << " Allocated memory (phase 2)"
             << std::endl;
#endif

    // ============================================
    //  Insert column indices into final structur
    // ============================================
    for ( i = 0; i < sysMatDim_-1; i++ ) {
      auxPtr = cIndex[i];
      k = rptrU_[i];
      for ( j = 0; j < rptrU_[i+1] - rptrU_[i]; j++ ) {
        cidxU_[k+j] = auxPtr[j];
      }
    }

    // ================
    //  Export pattern
    // ================

    // if the user desires it, we will export the pattern of the matrix
    // factor to a file in Matrix Market format
    bool saveFacToFile = false;
    std::string facFileName = "fac.out";
    bool savePatternOnly = false;
      
    if(xml_->Has("saveFacFile")) {
      saveFacToFile = true;
      xml_->GetValue("saveFacFile", facFileName, ParamNode::INSERT);
    }

    xml_->GetValue("savePatternOnly", savePatternOnly, ParamNode::INSERT);
    
    if( saveFacToFile && savePatternOnly ) {
      ExportFactorisation( facFileName.c_str(), false );
    }


    // ==========
    //  Clean up
    // ==========

    // Free auxilliary index array
    delete[] auxVec;
    for ( k = 0; k < nRows-1; k++ ) {
		  delete [] ( cIndex[k] );  cIndex[k] = NULL;
		}
    delete [] ( cIndex );  cIndex = NULL;


#ifdef DEBUG_LDLSOLVER_ANALYSE
    (*debug) << " Row pointer array:";
    for ( i = 0; i < nCols + 1; i++ ) {
      (*debug) << " " << rptrU_[i];
    }
    (*debug) << "\n Column index array:";
    for ( i = 0; i < rptrU_[nCols]; i++ ) {
      (*debug) << " " << cidxU_[i];
    }
    (*debug) << '\n' << std::endl;
#endif

  }


  // *************
  //   Factorise
  // *************
  template<typename T>
  void LDLSolver<T>::Factorise( const SCRS_Matrix<T> &sysMat ) {


    UInt k, i, j, numOffD;
    T elim;

#ifdef DEBUG_LDLSOLVER_FACTORISE
    (*debug) << " LDLSOLVER:\n Phase: FACTORISE\n" << std::endl;
#endif

    // Needed for writing progress report of factorisation
    UInt percentDone = 0;
    Double actDone = 0.0;

    // ==================
    //  Get Matrix Info
    // ==================

    // Query number of matrix rows and columns
    UInt nRows = sysMat.GetNumRows();
    UInt nCols = sysMat.GetNumCols();

    // Get hold of column index array
    const UInt *cidxA = sysMat.GetColPointer();

    // Get hold of row pointer index array
    const UInt *rptrA = sysMat.GetRowPointer();

    // Get hold of data array
    const T *dataA = sysMat.GetDataPointer();

    // =======================
    //  Init diagonal factor
    // =======================
    for ( i = 0; i < nCols; i++ ) {
      dataD_[i] = dataA[ rptrA[i] ];
    }

    // ==================
    //  Treat first row
    // ==================

    // Copy and scale first row
    for ( i = 0; i < rptrU_[1]; i++ ) {
      dataU_[i] = OpType<T>::invert( dataD_[0] ) * dataA[i+1];
    }

    // Perform update of diagonal factor for row one
    for ( i = 0; i < rptrU_[1]; i++ ) {
      dataD_[ cidxU_[i] ] -= dataU_[i] * dataD_[0] * dataU_[i];
    }


    // =================================
    //  Init auxilliary data structures
    // =================================


    // Generate markers for linked lists
    //take care: this elements still assume row and col indices
    //           starting from 1!!
    UInt scanListElem, scanListPrevElem;
    UInt activeListElem, activeListPrevElem;
    UInt listEnd = nCols+1;
    Integer listNoElem = -1;

    // Check that linked lists are empty
    if ( scanList_[0] != static_cast<Integer>(listEnd) ) {
      EXCEPTION( "LDLSolver::Factorise: Internal error: scanList_[0] = "
	       << scanList_[0] << " and not " << listEnd );
    }
    if ( activeList_[0] != static_cast<Integer>(listEnd) ) {
      EXCEPTION( "LDLSolver::Factorise: Internal error: activeList_[0] = "
	       << activeList_[0] << " and not " << listEnd );
    }

    // list iterator
    std::list<UInt>::iterator listIt;

    // Check, if the first row has an off-diagonal entry
    // (everything else would be very suspicious in the FEM setting)
    if ( rptrU_[1] > 1 ) {

      // add first row to scan list
      scanList_[0] = 1;
      scanList_[1] = listEnd;

      // Determine index of first off-diagonal entry in first row
      // in data array
      firstU_[0] = rptrU_[0];

      // If the corresponding column index is equal to one, we must add
      // the first row to the active list
      if ( cidxU_[ firstU_[0] ] == 1 ) {
        activeList_[0] = 1;
        activeList_[1] = listEnd;
      }
    }
    else {
      WARN("First row has no off-diagonal entry! I'm frightened! "
           << "Can this be correct my master?");
    }

    // ============
    //  Main Loop
    // ============
    for ( k = 1; k < nRows-1; k++ ) {

#ifdef DEBUG_LDLSOLVER_FACTORISE
      (*debug) << " ==================\n  Treating row " << k
               << "\n ==================\n";
      (*debug) << " + scan list contains row(s):";
      scanListElem = scanList_[0];
      while ( scanListElem != listEnd ) {
        (*debug) << " " << scanListElem;
        scanListElem = scanList_[scanListElem];
      }
      (*debug) << std::endl;
      (*debug) << " + active list contains row(s):";
      activeListElem = activeList_[0];
      while ( activeListElem != listEnd ) {
        (*debug) << " " << activeListElem;
        activeListElem = activeList_[activeListElem];
      }
      (*debug) << std::endl;
#endif

      // Keep user informed on progress
      actDone = (double)(k*100) / (double)nRows;
      actDone = (UInt)(actDone/10.0)*10;
      if ( actDone > percentDone ) {
        percentDone = (UInt)actDone;
      }

      // Copy k-th row of A into dense vector (omitting diagonal)
      for ( i = rptrA[k]+1; i < rptrA[k+1]; i++ ) {
        denseVec_[ cidxA[i] ] = dataA[i];
      }

      // All rows in active list contribute to current row
      activeListElem = activeList_[0];
      while ( activeListElem != listEnd ) {

	// Get data array index of row entry for column k
	i = firstU_[ activeListElem-1  ];

	// Get value of row entry for column k
	elim = dataU_[i];

	// Scale value with corresponding entry in diagonal factor
	elim = elim * dataD_[ activeListElem-1 ];

	// Add multiple of active row to current row
	for ( j = i + 1; j < rptrU_[activeListElem]; j++ ) {
	  denseVec_[ cidxU_[j] ] -= elim * dataU_[j];
	}

#ifdef DEBUG_LDLSOLVER_FACTORISE
        (*debug) << " + updated row by row " << activeListElem
                 << std::endl;
#endif

	// Proceed to next row
	activeListElem = activeList_[activeListElem];
      }

      // Fusion of two loops:
      //
      // - Insert dense vector entries into sparse structure
      //   (scaling them at the same time with d_kk)
      // - Update all entries in diagonal factor with index larger k
      //   (if the U_ki entry is not zero)
      for ( i = rptrU_[k]; i < rptrU_[k+1]; i++ ) {
        dataU_[i] = OpType<T>::invert( dataD_[k] ) * denseVec_[ cidxU_[i] ];
        denseVec_[ cidxU_[i] ] = 0.0;
        dataD_[ cidxU_[i] ] -= dataU_[i] * dataD_[k] * dataU_[i];
      }

      // -----------------------------------
      //  Update auxilliary data structures
      // -----------------------------------

      // Prepare list element pointers (both lists are in ascending
      // order, activeList_ is a subset of scanList_ and we keep
      // activeListElem >= scanListElem)
      scanListPrevElem = 0;
      scanListElem = scanList_[0];
      activeListPrevElem = 0;
      activeListElem = activeList_[0];

      // Loop over scan list
      while ( scanListElem < listEnd ) {

#ifdef DEBUG_LDLSOLVER_FACTORISE
        (*debug) << " + scanListElem = " << scanListElem;
#endif

	// Check if column index of first entry is equal to k
	// (in this case the row must also be in the active list
        // and activeListElem = scanListElem)
	if ( cidxU_[ firstU_[scanListElem-1] ] == k ) {

#ifdef DEBUG_LDLSOLVER_FACTORISE
        (*debug) << ": incr. firstU_ from " << firstU_[scanListElem-1]
                 << " / " << cidxU_[firstU_[scanListElem-1]]
                 << " --> ";
#endif

	  // If this was the last row entry, delete row from scan list
	  // and also from active list
	  if ( firstU_[scanListElem-1] + 1 == rptrU_[scanListElem] ) {

#ifdef DEBUG_LDLSOLVER_FACTORISE
            (*debug) << "nil (del from both lists)" << std::endl;
#endif

	    // Deletion from scan list
	    scanList_[scanListPrevElem] = scanList_[scanListElem];
	    scanList_[scanListElem] = listNoElem;
	    scanListElem = scanList_[scanListPrevElem];

	    // Deletion from active list
	    activeList_[activeListPrevElem] = activeList_[activeListElem];
	    activeList_[activeListElem] = listNoElem;
	    activeListElem = activeList_[activeListPrevElem];

            // Incrementing was done implicitely in deletion and
            // activeListElem >= scanListElem must still hold,
            // so procceed with loop
            continue;
	  }

	  // This was not the last entry in the row
	  else {

            // Increment the column index pointer
	    firstU_[scanListElem-1]++;

#ifdef DEBUG_LDLSOLVER_FACTORISE
            (*debug) << firstU_[scanListElem-1]
                     << " / " << cidxU_[firstU_[scanListElem-1]];
#endif

	    // Check whether we must delete row from active list now
	    if ( cidxU_[ firstU_[scanListElem-1] ] > k + 1 ) {
	      activeList_[activeListPrevElem] = activeList_[activeListElem];
	      activeList_[activeListElem] = listNoElem;
	      activeListElem = activeList_[activeListPrevElem];

#ifdef DEBUG_LDLSOLVER_FACTORISE
              (*debug) << " (active list del)";
#endif
	    }


	    // If there is no other off-diagonal entry, then we can delete
	    // the row from both lists
	    else if ( firstU_[scanListElem-1] + 1 == rptrU_[scanListElem] ) {

#ifdef DEBUG_LDLSOLVER_FACTORISE
	      (*debug) << " (del from both lists)" << std::endl;
#endif

	      // Deletion from scan list
	      scanList_[scanListPrevElem] = scanList_[scanListElem];
	      scanList_[scanListElem] = listNoElem;
	      scanListElem = scanList_[scanListPrevElem];

	      // Deletion from active list
	      activeList_[activeListPrevElem] = activeList_[activeListElem];
	      activeList_[activeListElem] = listNoElem;
	      activeListElem = activeList_[activeListPrevElem];

	      // Incrementing was done implicitely in deletion and
	      // activeListElem >= scanListElem must still hold,
	      // so procceed with loop
	      continue;
	    }

            // We did not delete the row, so we must advance the element
            // in the activeList_ (we know that here activeListElem < listEnd)
            else {
              activeListPrevElem = activeListElem;
              activeListElem = activeList_[activeListElem];
            }

#ifdef DEBUG_LDLSOLVER_FACTORISE
            (*debug) << std::endl;
#endif

            // Proceed to next row in scan list
            scanListPrevElem = scanListElem;
            scanListElem = scanList_[scanListElem];
            continue;
	  }
	}

	// So the column index of first entry was > k (the row was not in
	// the active list, i.e. scanListElem < activeListElem).
        // Check, whether it is (k+1)
	else if ( cidxU_[ firstU_[scanListElem-1] ] == k + 1 ) {

          // If there is still another off-diagonal entry with larger column
          // index available add row to active list
          if ( firstU_[scanListElem-1] + 1 < rptrU_[scanListElem] ) {

#ifdef DEBUG_LDLSOLVER_FACTORISE
            if ( activeListPrevElem > scanListElem ) {
              EXCEPTION("Error in active list add: (activeListPrevElem = "
                       << activeListPrevElem << ") < (scanListElem = "
                        << scanListElem << ")");
            }
#endif

            activeList_[ scanListElem       ] = activeListElem;
            activeList_[ activeListPrevElem ] = scanListElem;
            activeListPrevElem = scanListElem;

#ifdef DEBUG_LDLSOLVER_FACTORISE
            (*debug) << ": (active list add)";
#endif
          }

          // There is no other off-diagonal entry, so we do not
          // add the row to the activeList_ and can remove it
          // from the scanList_
          else {
            scanList_[scanListPrevElem] = scanList_[scanListElem];
            scanList_[scanListElem] = listNoElem;
            scanListElem = scanList_[scanListPrevElem];

#ifdef DEBUG_LDLSOLVER_FACTORISE
            (*debug) << ": (scan list del)" << std::endl;
#endif

            // Incrementing was done implicitely in deletion and
            // activeListElem >= scanListElem must still hold,
            // so procceed with loop
            continue;
          }
        }

        // Proceed to next row in scan list
        scanListPrevElem = scanListElem;
        scanListElem = scanList_[scanListElem];

#ifdef DEBUG_LDLSOLVER_FACTORISE
	(*debug) << std::endl;
#endif

      }


#ifdef DEBUG_LDLSOLVER_FACTORISE
      // We claim that activeListElem >= scanListElem always holds, so
      // now we should have activeListElem = listEnd
      if ( activeListElem != listEnd ) {
        EXCEPTION( "LDLSolver::Factorise: activeListElem = "
                   << activeListElem << ", but should be " << listEnd);
      }
#endif
      


      // --------------------------------------------
      //  Add k-th row to auxilliary data structures
      // --------------------------------------------

      // Determine number off off-diagonal entries in row
      numOffD = rptrU_[k+1] - rptrU_[k];

      // Make sure that row has any off-diagonal elements at all
      if ( numOffD > 0 ) {

        // Set pointer to first off-diagonal entry
        firstU_[k] = rptrU_[k];

	// If there are two off-diagonal entries, at least one must
	// have a column index > (k+1), so we must add row to scan list
	if ( numOffD > 1 ) {

	  // Addition is performed at the end (so the list will
	  // always be ordered with increasing row index)
	  // scanList_[k] = scanList_[0];
	  // scanList_[0] = k;
	  scanList_[scanListPrevElem] = k+1;
	  scanList_[k+1] = listEnd;

#ifdef DEBUG_LDLSOLVER_FACTORISE
	  (*debug) << " + added row " << k << " to scan list" << std::endl;
#endif

	  // If the column index of the first off-diagonal entry is (k+1)
	  // we must add it to the active list
	  if ( cidxU_[ firstU_[k] ] == (k+1) ) {

#ifdef DEBUG_LDLSOLVER_FACTORISE
	    (*debug) << " + added row " << k << " to active list" << std::endl;
#endif
	    // Addition is performed at the end
	    // activeList_[k] = activeList_[0];
	    // activeList_[0] = k;
	    activeList_[activeListPrevElem] = k+1;
	    activeList_[k+1] = listEnd;
	  }
        }
      }

#ifdef DEBUG_LDLSOLVER_FACTORISE
      (*debug) << " --> scanList_:";
      for ( UInt ii = 0; ii <= sysMatDim_; ii++ ) {
        (*debug) << " " << scanList_[ii];
      }
      (*debug) << std::endl;
      (*debug) << " --> activeList_:";
      for ( UInt ii = 0; ii <= sysMatDim_; ii++ ) {
        (*debug) << " " << activeList_[ii];
      }
      (*debug) << std::endl;
#endif
    }

  }


  // =======================
  //   ExportFactorisation
  // =======================
  template <typename T>
  void LDLSolver<T>::ExportFactorisation( const char *fname,
                                          bool patternOnly ) {


    UInt i, j, k;
    T aux;

    // ====================
    //   Open output file
    // ====================
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      EXCEPTION( "LDLSolver::ExportFactorisation: Cannot open file "
               << fname << " for writing!" );
    }

    // =====================
    //   Write file header
    // =====================

    // Matrix Market Format Specification
    std::string myFormat;
    if ( patternOnly == true ) {
      myFormat = "pattern";
    }
    else {
      if ( EntryType<T>::M_EntryType == BaseMatrix::DOUBLE ) {
        myFormat = "real";
      }
      else if ( EntryType<T>::M_EntryType == BaseMatrix::COMPLEX ) {
        myFormat = "complex";
      }
      else {
        EXCEPTION( "LDLSolver::ExportFactorisation: Cannot identify "
                 << "template parameter" );
      }
    }
    fprintf( fp, "%%%%MatrixMarket matrix coordinate %s general\n",
             myFormat.c_str() );

    // Comment
    if ( patternOnly == true ) {
      fprintf( fp, "%%\n%% Sparsity pattern of LDL^T factorisation matrix " );
      fprintf( fp, "F = D + L^T\n" );
      fprintf( fp, "%% computed by LDLSolver\n%%\n" );
    }
    else {
      fprintf( fp, "%%\n%% LDL^T factorisation matrix F = D + L^T " );
      fprintf( fp, " computed by LDLSolver\n%%\n" );
    }

    // Information on number of rows, columns and entries
    Integer dof = BlockSize<T>::size;
    fprintf( fp, "%d\t%d\t%d\n", sysMatDim_ * dof, sysMatDim_ * dof,
             (rptrU_[sysMatDim_] + sysMatDim_) * dof * dof );

    // ======================
    //   Write entries of D
    // ======================
    for ( i = 1; i <= sysMatDim_; i++ ) {
      for ( j = 1; (int) j <= dof; j++ ) {
        for ( k = 1; (int) k <= dof; k++ ) {
          fprintf( fp, "%6d\t%6d\t", (i-1) * dof + j, (i-1) * dof + k );
          if ( patternOnly == false ) {
            aux = OpType<T>::invert( dataD_[i-1] );
            OpType<T>::ExportEntry( aux, j-1, k-1, fp );
          }
          fprintf( fp, "\n" );
        }
      }
    }

    // ======================
    //   Write entries of U
    // ======================
    Integer ib, jb;
    UInt nblocks;
    
    // loop over all block rows
    for ( i = 0; i < sysMatDim_; i++ ) {

      // get number of blocks in i-th row
      nblocks = rptrU_[i+1] - rptrU_[i];

      // loop over all blocks in this row
      for ( j = 0; j < nblocks; j++ ) {

        // loop over block entries
        for ( ib = 0; ib < dof; ib++ ) {
          for ( jb = 0; jb < dof; jb++ ) {

            // store row and column index
            fprintf( fp, "%6d\t%6d\t", i * dof + ib + 1,
                     ( cidxU_[rptrU_[i]+j] ) * dof + jb + 1);

            // store non-zero entry
            OpType<T>::ExportEntry( dataU_[rptrU_[i]+j], ib, jb, fp );
            fprintf( fp, "\n" );
          }
        }
      }
    }

    // =====================
    //   Close output file
    // =====================
    if ( fclose( fp ) == EOF ) {
      WARN("LDLSolver::ExportFactorisation: Could not close file "
           << fname << " after writing!");
    }
  }

  // Explicit template instantiation
  template class LDLSolver<Double>;
  template class LDLSolver<Complex>;
  
}
