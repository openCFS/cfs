#ifndef CROUT_LU_CC
#define CROUT_LU_CC

#include <algorithm>
#include <iterator>

#include "General/Environment.hh"
#include "Utils/tools.hh"
#include "MatVec/opdefs.hh"
#include "MatVec/BaseMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "OLAS/utils/math/CroutLU.hh"

// Be commenting in or out the macros below debugging of the CroutLU class
// can be made more or less fine-grain
#ifdef DEBUG_CROUTLU
#define DEBUG_CROUTLU_ILU
#define DEBUG_CROUTLU_SOLVE
#endif

namespace CoupledField {

  // =======================
  //   Default Constructor
  // =======================
  template <typename T>
  CroutLU<T>::CroutLU() {


    // We do not yet have performed a factorisation
    storingFactors_ = false;

    // Set estimate for growth of sparsity pattern
    // This is quite static of course. The value can
    // be adapted by the constructors of the derived
    // classes. ILUTP_Precond e.g. can do this much
    // better.
    memGrowthEstimate_ = 5;

    // By default we do not set fill-values in GetColOfA()
    // or GetRowOfA() or compute the level of fill-in,
    // since only ILUK needs this
    compFillLevels_ = false;

  }


  // ======================
  //   Default Destructor
  // ======================
  template <typename T>
  CroutLU<T>::~CroutLU() {
    DeleteDynamicDataStructures();
  }


  // ===============================
  //   DeleteDynamicDataStructures
  // ===============================
  template <typename T>
  void CroutLU<T>::DeleteDynamicDataStructures() {
  }


  // =========
  //   Reset
  // =========
  template <typename T>
  void CroutLU<T>::Reset() {


    // Free dynamically allocate memory
    DeleteDynamicDataStructures();

    // Delete old factors
    // Note that this does not free the memory reserved by the STL vectors
    cptrL_.clear();
    ridxL_.clear();
    entryL_.clear();

    rptrU_.clear();
    cidxU_.clear();
    entryU_.clear();

    // Delete size information for system matrix
    sysMatDim_ = 0;
    sysMatNNZ_ = 0;

  }


  // =============
  //   Factorise
  // =============
  template <typename T>
  void CroutLU<T>::Factorise( CRS_Matrix<T> &sysMat ) {

    
    // Reset the object before over-writing an existing factorisation
    if ( storingFactors_ == true ) {
      Reset();
      storingFactors_ = false;
    }

    // Obtain size information of the matrix
    sysMatDim_ = sysMat.GetNumRows();
    if ( sysMatDim_ != sysMat.GetNumCols() ) {
      EXCEPTION( "CroutLU::Factorise encountered a "
               << sysMatDim_ << " x " << sysMat.GetNumRows()
               << "system matrix. Should be square, however." );
    }
    sysMatNNZ_ = sysMat.GetNnz();

    // If the matrix is not in LEX or LEX_DIAG_FIRST sub-format, we prefer
    // to have it in LEX, which is the default if it comes from openCFS.
    if ( sysMat.GetCurrentLayout() == CRS_Matrix<T>::UNSORTED ) {
      sysMat.ChangeLayout( CRS_Matrix<T>::LEX );
    }

    // ***************************************
    //   Initialise internal data structures
    // ***************************************

    // In order to prevent constant re-allocation we reserve space a priori
    // for the representation of the L and U factors and the other data
    // structures. For some we know how large they will be. For others
    // we use an initial estimate based on memGrowthEstimate_.

    // Note that we use one-based operations and thus need to reserve an
    // additional element. And insert a dummy value up front.

    // Compute an estimate for the number of entries in a row / column
    // and per factor L and U
    UInt estNZPerRow = (sysMatNNZ_ * memGrowthEstimate_) / sysMatDim_;
    UInt estNZPerFac = (sysMatNNZ_ + sysMatDim_ ) / 2 * memGrowthEstimate_;

    // Avoid doing over-allocation
    estNZPerRow = estNZPerRow > sysMatDim_ ? sysMatDim_ : estNZPerRow;
    estNZPerFac = estNZPerFac > sysMatDim_ * ( sysMatDim_ + 1 ) / 2 ?
      sysMatDim_ * ( sysMatDim_ + 1 ) / 2 : estNZPerFac;

#ifdef DEBUG_CROUTLU_ILU
    (*debug) << "\nCroutLU: Memory pre-allocation:\n"
             << "estNZPerRow = " << estNZPerRow << "\n"
             << "estNZPerFac = " << estNZPerFac << "\n"
             << std::endl;
#endif

    // Initialise the data structure for storing the matrix factor U,
    // we use compressed row storage format for this
    cidxU_.reserve( estNZPerFac + 1 );
    //    cidxU_.push_back( 0 );

    rptrU_.reserve( sysMatDim_ + 2 );
    rptrU_.push_back( 0 );
    //rptrU_.push_back( 1 );

    entryU_.reserve( estNZPerFac + 1 );
    //    entryU_.push_back( 0 );

    // Initialise the data structure for storing the matrix factor L,
    // we use compressed column storage format for this
    ridxL_.reserve( estNZPerFac + 1 );
    //ridxL_.push_back( 0 );

    cptrL_.reserve( sysMatDim_ + 2 );
    cptrL_.push_back( 0 );
    //cptrL_.push_back( 1 );

    entryL_.reserve( estNZPerFac + 1 );
    //    entryL_.push_back( 0 );

    // The following containers are the two dense vectors we need for our
    // Crout implementation and two vectors that are used to keep track of
    // the non-zero entries in these two dense vectors. We initialise the
    // dense vectors with zeros
    std::vector<T> vecZ;
    std::vector<T> vecW;
    std::vector<UInt> vecZFill;
    std::vector<UInt> vecWFill;

    vecZ.reserve( sysMatDim_ + 1 );
    vecW.reserve( sysMatDim_ + 1 );

    // Pre-allocate memory for index set of non-zero entries
    vecZFill.reserve( estNZPerRow );
    vecWFill.reserve( estNZPerRow );

    T zero = 0.0;
    for ( UInt k = 0; k <= sysMatDim_; k++ ) {
      vecZ.push_back( zero );
      vecW.push_back( zero );
    }

    // Setup the data structures for accessing the columns of U and the rows
    // of L. We use what Saad calls incomplete linked lists:
    //
    // For each column i of U we have a vector listU_[i] that stores the
    // indices of all rows of U that have non-zero entries in column i.
    // These lists are build dynamically during the update of the firstU
    // vector. L is treated analogously
    std::vector<Integer> *listL;
    std::vector<Integer> *listU;
    NEWARRAY( listL, std::vector<Integer>, sysMatDim_ );
    NEWARRAY( listU, std::vector<Integer>, sysMatDim_ );
    for ( UInt i = 0; i < sysMatDim_; i++ ) {
      listL[i].reserve( estNZPerRow );
      listU[i].reserve( estNZPerRow );
      listL[i].push_back( -1 );
      listU[i].push_back( -1 );
    }

    // Initialise the index vectors showing where the non-zero entries with
    // column (for U) / row (for L) indices larger then k (see outer most loop
    // below) start. We can update these vectors sucessively during the outer
    // most loop and need not initialise all here now.
    std::vector<Integer> firstL;
    std::vector<Integer> firstU;

    firstL.reserve( sysMatDim_ + 1 );
    firstU.reserve( sysMatDim_ + 1 );

//     firstL.push_back( 0 );
//     firstU.push_back( 0 );

    // Needed for writing progress report of factorisation
    UInt percentDone = 0;
    Double actDone = 0.0;

    // *************************
    //   Extra things for ILUK
    // *************************
    if ( compFillLevels_ == true ) {

      // Pre-allocate memory for the fill-levels
      vecZLevel_.reserve( sysMatDim_ + 1 );
      vecWLevel_.reserve( sysMatDim_ + 1 );

      // Set initial fill-levels (to maximal value that can
      // never be attained)
      if ( vecZLevel_.empty() == true ) {
        vecZLevel_.resize( sysMatDim_ + 1, sysMatDim_ + 1 );
        vecWLevel_.resize( sysMatDim_ + 1, sysMatDim_ + 1 );
      }
      else {
        for ( UInt i = 0; i < vecZLevel_.size(); i++ ) {
          vecZLevel_[i] = sysMatDim_ + 1;
          vecWLevel_[i] = sysMatDim_ + 1;
        }
        for ( UInt i = vecZLevel_.size(); i <= sysMatDim_; i++ ) {
          vecZLevel_.push_back( sysMatDim_ + 1 );
          vecWLevel_.push_back( sysMatDim_ + 1 );
        }
      }

      // Pre-allocate memory for fill-levels of factor entries
      levelU_.reserve( estNZPerFac + 1 );
      levelL_.reserve( estNZPerFac + 1 );
      if ( levelU_.empty() ) {
        levelU_.push_back( 0 );
        levelL_.push_back( 0 );
      }

    }

    // ***********************
    //   Crout Factorisation
    // ***********************

    UInt k, i, j, memIndex, idx, auxLvl, lvlParent1, lvlParent2;
    UInt rightMostU, lowestL;
    T aux;

    // Outermost loop
    for ( k = 0; k < sysMatDim_; k++ ) {

#ifdef DEBUG_CROUTLU_ILU
      (*debug) << "--------------------------------------\n"
               << " K = " << k << '\n'
               << "--------------------------------------" << std::endl;

      PrintVector( firstU, "firstU" );
      PrintVector( firstL, "firstL" );
#endif

      // Keep user informed on progress
      actDone = (double)(k*100) / (double)sysMatDim_;
      actDone = (UInt)(actDone/10.0)*10;
      if ( actDone > percentDone ) {
        percentDone = (UInt)actDone;
      }


      // *********************************
      //   COMPUTATION of k-th ROW of U
      //         (as dense vector)
      // *********************************

      // Copy the non-zero entries of the k-th row of the system matrix right
      // of and including the diagonal into a dense vector
      rightMostU = GetRowOfA( sysMat, vecZ, k );

      if ( compFillLevels_ == false ) {

        // Compute row k of U by performing elimination with all rows
        // above row k. However, we only need those rows idx, for which the
        // elimination factor l(k,idx) is not zero. The indices of the columns
        // of L for which this is the case are stored in listL[k]
        for ( i = 1; i < listL[k].size(); i++ ) {

          // Obtain row index from list
          idx = listL[k][i];

          // Obtain the non-zero elimination factor. This is the first non-
          // zero entry in column idx of L with row index greater equal to k
          aux = entryL_[firstL[idx]];

          // Elimination factor l(k,idx) is not zero, so we must treat all
          // non-zero entries in row idx of U, with column indices starting
          // from k. The first one is stored in firstU[idx]
          for ( j = firstU[idx]; j < rptrU_[idx+1]; j++ ) {
            vecZ[ cidxU_[j] ] -= aux * entryU_[j];
          }

	  // Check if inded of maximal (i.e. right-most) non-zero entry
	  // in row of U has increased
	  rightMostU = cidxU_[j-1] > rightMostU ? cidxU_[j-1] : rightMostU;
        }
      }

      else {

        // Compute row k of U by performing elimination with all rows
        // above row k. However, we only need those rows idx, for which the
        // elimination factor l(k,idx) is not zero. The indices of the columns
        // of L for which this is the case are stored in listL[k]
        for ( i = 1; i < listL[k].size(); i++ ) {

          // Obtain row index from list
          idx = listL[k][i];

          // Obtain the non-zero elimination factor. This is the first non-
          // zero entry in column idx of L with row index greater equal to k
          aux = entryL_[ firstL[idx] ];
          lvlParent1 = levelL_[ firstL[idx] ];

          // Elimination factor l(k,idx) is not zero, so we must treat all
          // non-zero entries in row idx of U, with column indices starting
          // from k. The first one is stored in firstU[idx]
          lvlParent2 = sysMatDim_+1;
          for ( j = firstU[idx]; j < rptrU_[idx+1]; j++ ) {

            // Add to entry
            vecZ[ cidxU_[j] ] -= aux * entryU_[j];
 
            // Correct fill-level
            lvlParent2 = levelU_[j];
            auxLvl = lvlParent1 + lvlParent2;
            if ( auxLvl < vecZLevel_[ cidxU_[j] ] ) {
              vecZLevel_[ cidxU_[j] ] = auxLvl;
            }
          }

	  // Check if inded of maximal (i.e. right-most) non-zero entry
	  // in row of U has increased
	  rightMostU = cidxU_[j-1] > rightMostU ? cidxU_[j-1] : rightMostU;

        }
      }


      // ***********************************
      //   COMPUTATION of k-th COLUMN of L
      //         (as dense vector)
      // ***********************************

      // Copy the non-zero entries of the k-th column of the system matrix
      // below the diagonal into a dense vector
      lowestL = GetColOfA( sysMat, vecW, k );

      if ( compFillLevels_ == false ) {

        // Compute column k of L by performing elimination with all columns
        // left of column k. However, we only need those columns idx, for
        // which there is a non-zero entry u(idx,k) in the k-th column of U.
        // elimination factor l(k,idx) is not zero. The indices of the rows
        // of U for which this is the case are stored in listU[k]
        for ( i = 1; i < listU[k].size(); i++ ) {

          // Obtain column index from list
          idx = listU[k][i];

          // Obtain the non-zero entry u(idx,k)
          aux = entryU_[firstU[idx]];

          // The entry u(k,idx) is not zero, so we must treat all non-zero
          // entries in column idx of L, with with row indices starting
          // from k. The first of these is stored in firstL[idx]
          for ( j = firstL[idx]+1; j < cptrL_[idx+1]; j++ ) {
            vecW[ ridxL_[j] ] -= aux * entryL_[j];
          }

	  // Check if index of maximal (i.e. bottom-most) non-zero entry
	  // in column of L has increased
         // std::cerr << "j-1: "<<j-1<<" size of ridxL_:"<< ridxL_.size() << std::endl; 
         // std::cerr << ridxL_[j-1]<<" > "<<lowestL << std::endl;   
	  lowestL = ridxL_[j-1] > lowestL ? ridxL_[j-1] : lowestL;
        }
      }

      else {

        // Compute column k of L by performing elimination with all columns
        // left of column k. However, we only need those columns idx, for
        // which there is a non-zero entry u(idx,k) in the k-th column of U.
        // elimination factor l(k,idx) is not zero. The indices of the rows
        // of U for which this is the case are stored in listU[k]
        for ( i = 1; i < listU[k].size(); i++ ) {

          // Obtain column index from list
          idx = listU[k][i];

          // Obtain the non-zero entry u(idx,k)
          aux = entryU_[ firstU[idx] ];
          lvlParent1 = levelU_[ firstU[idx] ];

          // The entry u(k,idx) is not zero, so we must treat all non-zero
          // entries in column idx of L, with with row indices starting from
          // k. The first of these is stored in firstL[idx]
          for ( j = firstL[idx]+1; j < cptrL_[idx+1]; j++ ) {

            // Add to entry
            vecW[ ridxL_[j] ] -= aux * entryL_[j];

            // Correct fill-level
            lvlParent2 = levelL_[j];
            auxLvl = lvlParent1 + lvlParent2;
            if ( auxLvl < vecWLevel_[ ridxL_[j] ] ) {
              vecWLevel_[ ridxL_[j] ] = auxLvl;
            }
          }

	  // Check if inded of maximal (i.e. bottom-most) non-zero entry
	  // in column of L has increased
          lowestL = ridxL_[j-1] > lowestL ? ridxL_[j-1] : lowestL;

        }
      }


      // **********************
      //   GENERATE INDEX SET
      // **********************
      for ( i = k; i <= rightMostU; i++ ) {
        if ( vecZ[i] != 0.0 ) {
          vecZFill.push_back( i );
        }
      }
      for ( i = k + 1; i <= lowestL; i++ ) {
        if ( vecW[i] != 0.0 ) {
          vecWFill.push_back( i );
        }
      }


      // ************************************
      //   APPLICATION of DROPPING STRATEGY
      // ************************************
      DropEntries( k, vecZ, vecZFill, vecW, vecWFill );


      // ****************************************
      //   INSERT NEW ROW and COLUMN into L / U
      // ****************************************

      if ( compFillLevels_ == false ) {

        // Insert new row into U (zeroing at the same time vecZ)
        memIndex = 0;
        for ( std::vector<UInt>::iterator it = vecZFill.begin();
              it != vecZFill.end(); it++ ) {
          entryU_.push_back( vecZ[ *it ] );
          cidxU_.push_back( *it );
          vecZ[ *it ] = 0.0;
          memIndex++;
        }
        rptrU_.push_back( rptrU_[k] + memIndex );

        // Insert new column into L (zeroing at the same time vecW)
        // we must devide by u(k,k) at this point
        memIndex = 0;
        for ( std::vector<UInt>::iterator it = vecWFill.begin();
              it != vecWFill.end(); it++ ) {
          entryL_.push_back( vecW[ *it ] / entryU_[ rptrU_[k] ] );
          ridxL_.push_back( *it );
          vecW[ *it ] = 0.0;
          memIndex++;
        }
        cptrL_.push_back( cptrL_[k] + memIndex );
      }

      else {

        // Insert new row into U (zeroing at the same time vecZ)
        memIndex = 0;
        for ( std::vector<UInt>::iterator it = vecZFill.begin();
              it != vecZFill.end(); it++ ) {
          entryU_.push_back( vecZ[ *it ] );
          levelU_.push_back( vecZLevel_[ *it ] );
          cidxU_.push_back( *it );
          vecZ[ *it ] = 0.0;
          vecZLevel_[ *it ] = sysMatDim_+1;
          memIndex++;
        }
        rptrU_.push_back( rptrU_[k] + memIndex );

        // Insert new column into L (zeroing at the same time vecW)
        // we must devide by u(k,k) at this point
        memIndex = 0;
        for ( std::vector<UInt>::iterator it = vecWFill.begin();
              it != vecWFill.end(); it++ ) {
          entryL_.push_back( vecW[ *it ] / entryU_[ rptrU_[k] ] );
          levelL_.push_back( vecWLevel_[ *it ] );
          ridxL_.push_back( *it );
          vecW[ *it ] = 0.0;
          vecWLevel_[ *it ] = sysMatDim_+1;
          memIndex++;
        }
        cptrL_.push_back( cptrL_[k] + memIndex );

      }

      // Empty the index sets for dense auxilliary vectors
      vecZFill.clear();
      vecWFill.clear();


#ifdef DEBUG_CROUTLU_ILU
      PrintVector( rptrU_ , "rptrU"  );
      PrintVector( cidxU_ , "cidxU"  );
      PrintVector( entryU_, "entryU" );
      PrintVector( cptrL_ , "cptrL"  );
      PrintVector( ridxL_ , "ridxL"  );
      PrintVector( entryL_, "entryL" );
#endif


      // ************************************
      //   UPDATE AUXILLARY DATA STRUCTURES
      // ************************************

      // We can empty the lists for the k-th column of U / row of L, since
      // they are no longer required.
      // listU[k].empty();

      for ( i = 0; i < k; i++ ) {

        // std::cerr << " k = " << k << ", i = " << i << std::endl;

        // Update the index vectors for U (old rows)
        // and the incomplete linked lists too. Note
        // that the last row of U for which we will
        // ever need firstU is (n-1)

        // If old column index is equal to k we must advance it,
        // since it will be too small for the next iteration (k+1)
        // Guard against out-of-bounds access: firstU[i] may equal rptrU_[i+1]
        // (past the last entry of row i) if it was fully advanced in a prior
        // iteration. In that case no entry with column index k can exist.
        if ( firstU[i] < static_cast<Integer>(rptrU_[i+1]) && cidxU_[ firstU[i] ] == k ) {

            // advance pointer
            firstU[i]++;

            // update incomplete linked list
            if ( firstU[i] < static_cast<Integer>(rptrU_[i+1]) ) {
	      listU[ cidxU_[ firstU[i] ] ].push_back( i );

#ifdef DEBUG_CROUTLU_ILU
	      (*debug) << " Added row " << i << " to list for col "
		       << cidxU_[ firstU[i] ] << " of U" << std::endl;
#endif

            }
        }

        // Update the index vectors for L (old columns)
        // and the incomplete linked lists too. Note
        // that the last column of L for which we will
        // ever need firstL is (n-1)

        // If old row index is equal to k we must advance it,
        // since it will be too small for the next iteration (k+1)
        // Guard against out-of-bounds access: firstL[i] may equal cptrL_[i+1]
        // (past the last entry of column i) if it was fully advanced in a prior
        // iteration. In that case no entry with row index k can exist.
        if ( firstL[i] < static_cast<Integer>(cptrL_[i+1]) && ridxL_[ firstL[i] ] == k ) {

            // advance pointer
            firstL[i]++;

            // update incomplete linked list
            if ( firstL[i] < static_cast<Integer>(cptrL_[i+1]) ) {
              listL[ ridxL_[ firstL[i] ] ].push_back( i );

#ifdef DEBUG_CROUTLU_ILU
              (*debug) << " Added col " << i << " to list for row "
                       << ridxL_[ firstL[i] ] << " of L" << std::endl;
#endif

            }
        }
      }

      // Add pointer to first entry with column index >= (k+1) for
      // the new k-th row of U, also update linked list. We also must
      // allow for the case that there is no such entry in the row.
      firstU.push_back( -1 );
      for ( i = rptrU_[k]; i < rptrU_[k+1]; i++ ) {
        if ( cidxU_[i] > k ) {
          firstU[k] = i;
          listU[ cidxU_[i] ].push_back( k );

#ifdef DEBUG_CROUTLU_ILU
          (*debug) << " Added row " << k << " to list for col "
                   << cidxU_[i] << " of U" << std::endl;
#endif

          break;
        }
      }
      if ( firstU[k] == -1 ) {
        firstU[k] = rptrU_[k];
      }

      // Add pointer to first entry with row index >= (k+1) for the
      // new k-th column of L, also update linked list. We also must
      // allow for the case that there is no such entry in the column.
      firstL.push_back( -1 );
      for ( i = cptrL_[k]; i < cptrL_[k+1]; i++ ) {
        if ( ridxL_[i] > k ) {
          firstL[k] = i;
          listL[ ridxL_[i] ].push_back( k );

#ifdef DEBUG_CROUTLU_ILU
          (*debug) << " Added col " << k << " to list for row "
                   << ridxL_[i] << " of L" << std::endl;
#endif

          break;
        }
      }
      if ( firstL[k] == -1 ) {
         firstL[k] = cptrL_[k];
       }

    }

    // *******************************************
    //   Prepare L/U for forward/backward solves
    // *******************************************

    // We replace the diagonal entry of U by its inverse
    // so we don't need a division in the back substitutions
    for ( i = 0; i < sysMatDim_; i++ ) {
      entryU_[ rptrU_[i] ] = 1.0 / entryU_[ rptrU_[i] ];
    }

    // ********************************
    //   Do some administrative stuff
    // ********************************

    // Set status variable
    storingFactors_ = true;

    // Free dynamically allocated memory
    delete [] ( listL );  listL  = NULL;
    delete [] ( listU );  listU  = NULL;

  }


  // =========
  //   Solve
  // =========
  template <typename T>
  void CroutLU<T>::Solve( const Vector<T> &y, Vector<T> &x ) const {


#ifdef DEBUG_CROUTLU_SOLVE
    PrintVector( rptrU_ , "rptrU"  );
    PrintVector( cidxU_ , "cidxU"  );
    PrintVector( entryU_, "entryU" );
    PrintVector( cptrL_ , "cptrL"  );
    PrintVector( ridxL_ , "ridxL"  );
    PrintVector( entryL_, "entryL" );
#endif

    // Check that we have a factorisation at hand and that the
    // dimensions agree
    if ( storingFactors_ == false ) {
      EXCEPTION( "CroutLU::Solve cannot perform a solve, since there are no "
               << "factors L and U available. Call Factorise() first" );
    }
    if ( x.GetSize() != sysMatDim_ || y.GetSize() != sysMatDim_ ) {
      EXCEPTION( "CroutLU::Solve: The size of the input vectors is "
               << "x = " << x.GetSize() << " and y = " << y.GetSize()
               << ", while matrix is " << sysMatDim_ << " x " << sysMatDim_ );
    }

    // We solve LU x = y by a pair of forward/backward substitutions
    // (1) L z = y
    // (2) U x = z
    // The intermediate vector z is omitted in the implemenation, since
    // the step (2) can be done in place.

    UInt i, k;


    // *******************************
    //   Forward substitution with L
    // *******************************

    // L is stored in CCS format, where the diagonal entries are omitted,
    // since they are all equal to one. To solve L z = y we make use of
    // the representation y = z1 * l1 + z2 * l2 + ... + zn  * ln, where
    // zk is the k-th component of z and lk the k-th column vector of L
    // So after zk was computed we replace y by y - zk * lk

#ifdef DEBUG_CROUTLU_SOLVE
    (*debug) << "Starting forward substitution" << std::endl;
#endif

    for ( k = 0; k < sysMatDim_; k++ ) {
      x[k] = y[k];
    }

    for ( k = 0; k < sysMatDim_; k++ ) {
      for ( i = cptrL_[k]; i < cptrL_[k+1]; i++ ) {
        x[ ridxL_[i] ] -= entryL_[i] * x[k];
      }
    }


    // ********************************
    //   Backward substitution with U
    // ********************************

    // U is stored in CRS format, where the entries of a row are in ascending
    // lexicographical ordering w.r.t. their column indices. Thus, the
    // diagonal entry is stored first. The Factorise() method also already
    // inverted it.

#ifdef DEBUG_CROUTLU_SOLVE
    (*debug) << "Starting backward substitution" << std::endl;
#endif

    T aux;
    k = sysMatDim_;
    while ( k >0 ) {
      //    for ( k = sysMatDim_; k >= 1; k-- ) {
      k--;
      aux = x[k];
      for ( i = rptrU_[k] + 1; i < rptrU_[k+1]; i++ ) {
        aux -= entryU_[i] * x[ cidxU_[i] ];
      }
      x[k] = aux * entryU_[ rptrU_[k] ];
    }
  }


  // =============
  //   GetRowOfA
  // =============
  template <typename T>
  inline UInt CroutLU<T>::GetRowOfA( const CRS_Matrix<T> &sysMat,
                                     std::vector<T> &vecZ, UInt k ) {


    // Copy the non-zero entries of the k-th row of the system matrix
    // into the respective positions in the auxilliary full vector.
    //
    // We assume here that the output vectors vecZ and vecZFill are
    // empty on input

    // column index of right-most non-zero entry in row k
    UInt rightMost = 0;

    // Taunt the C++ concept of data encapsulation and obtain direct
    // access to the internal data structures of the system matrix
    const UInt *rptrA  = sysMat.GetRowPointer();
    const UInt *cidxA  = sysMat.GetColPointer();
    const UInt *diagA  = sysMat.GetDiagPointer();
    const T       *entryA = sysMat.GetDataPointer();

    // --------------------------------------------------
    // CASE 1: CRS_Matrix is in LEX_DIAG_FIRST sub-format
    // --------------------------------------------------
    if ( sysMat.GetCurrentLayout() == CRS_Matrix<T>::LEX_DIAG_FIRST ) {

      if ( compFillLevels_ == false ) {

        // Copy the diagonal entry
        vecZ[ k ] = entryA[ rptrA[k] ];

        // Copy the off-diagonal entries on the right
        UInt i;
        for ( i = rptrA[k] + 1; i < (UInt)rptrA[k+1]; i++ ) {
          if ( cidxA[i] >= k ) {
            vecZ[ cidxA[i] ] = entryA[i];
          }
        }

        // Assign right-most column index
        rightMost = cidxA[i-1];
      }

      else {

        // Copy the diagonal entry
        vecZ[ k ] = entryA[ rptrA[k] ];
        vecZLevel_[ k ] = 0;

        // Copy the off-diagonal entries on the right
        UInt i;
        for ( i = rptrA[k] + 1; i < (UInt)rptrA[k+1]; i++ ) {
          if ( cidxA[i] >=  k ) {
            vecZ[ cidxA[i] ] = entryA[i];
            vecZLevel_[ cidxA[i] ] = 0;
          }
        }

        // Assign right-most column index
        rightMost = cidxA[i-1];
      }
    }

    // ---------------------------------------
    // CASE 2: CRS_Matrix is in LEX sub-format
    // ---------------------------------------
    else if ( sysMat.GetCurrentLayout() == CRS_Matrix<T>::LEX ) {

      // Copy entries starting from diagonal
      if ( compFillLevels_ == false ) {
        for ( UInt i = (UInt)diagA[k]; i < (UInt)rptrA[k+1]; i++ ) {
          vecZ[ cidxA[i] ] = entryA[i];
        }
      }

      // Copy entries starting from diagonal and set fill-level to zero
      else {
        for ( UInt i = (UInt)diagA[k]; i < (UInt)rptrA[k+1]; i++ ) {
          vecZ[ cidxA[i] ] = entryA[i];
          vecZLevel_[ cidxA[i] ] = 0;
        }
      }

      // Assign right-most column index
      rightMost = cidxA[rptrA[k+1]-1];
    }

    // ----------------------------------------------
    // CASE 3: CRS_Matrix is in some other sub-format
    // ----------------------------------------------
    else {

      EXCEPTION( "Got a problem in CroutLU<T>::GetRowOfA. "
               << "The system matrix is neither in LEX nor LEX_DIAG_FIRST "
               << "sub-format of CRS storage. Cannot re-order myself, "
               << "however, due to 'const' qualifier!@$%&" );

      // In this case we sort the matrix to LEX and solve
      // problem by recursion
      // sysMat.ChangeLayout( CRS_Matrix<T>::LEX );
      // rightMost = GetRowOfA( sysMat, vecZ, k );
    }

    // Return right-most column index
    return rightMost;

  }


  // =============
  //   GetColOfA
  // =============
  template <typename T>
  inline UInt CroutLU<T>::GetColOfA( const CRS_Matrix<T> &sysMat,
                                     std::vector<T> &vecW, UInt k ) {


    // Our implementation of the Crout algorithm requires that we successively
    // copy the lower part of each column of the system matrix into a full
    // vector in order to compute the respective column of L. This operation
    // is a little bit tricky for a system matrix stored in CRS format.
    //
    // However, if the CRS_Matrix is stored either in LEX_DIAG_FIRST or
    // LEX_DIAG format, we can improve performance by introducing an
    // auxilliary index vector.
    //
    // This vector points to the first entry in a row that is a possible
    // candidate to be copied. If the entry is copied, then the pointer is
    // increased. Thus, we only need one check of a column index per row for
    // each column copy.
    static std::vector<UInt> firstUncheckedAij;

    // Largest row index of non-zero entry in k-th column
    UInt lowest = k;

    // Taunt the C++ concept of data encapsulation and obtain a direct access
    // to the internal data structures of the system matrix
    const UInt *rptrA  = sysMat.GetRowPointer();
    const UInt *cidxA  = sysMat.GetColPointer();
    const UInt *diagA  = sysMat.GetDiagPointer();
    const T       *entryA = sysMat.GetDataPointer();


    // --------------------------------------------------
    // CASE 1: CRS_Matrix is in LEX_DIAG_FIRST sub-format
    // --------------------------------------------------
    if ( sysMat.GetCurrentLayout() == CRS_Matrix<T>::LEX_DIAG_FIRST ) {

      // Initialise auxilliary index vector
      if ( k == 0 ) {
        firstUncheckedAij.clear();
        firstUncheckedAij.reserve( sysMatDim_ + 1 );
        firstUncheckedAij.push_back( 0 );
        firstUncheckedAij.push_back( 0 );
        for ( UInt i = 1; i < sysMatDim_; i++ ) {
          firstUncheckedAij.push_back( rptrA[i] + 1 );
        }
      }

      if ( compFillLevels_ == false ) {

        // Copy the non-zero entries of the k-th column (below the diagonal)
        // of the system Matrix into the respective positions in the
        // auxilliary full vector.
        UInt i, cidx;
        for ( i = k + 1; i < sysMatDim_; i++ ) {

          // Avoid "overflow" into next row
          if ( firstUncheckedAij[i] < (UInt)rptrA[i+1] ) {

            cidx = cidxA[ firstUncheckedAij[i] ];

            if ( cidx == k ) {

              // This entry must be copied into the full auxilliary vector
              vecW[ i ] = entryA[ firstUncheckedAij[i] ];

              // This entry is a candidate for the largest row index
              lowest = i;

              // Increment pointer to next smallest unchecked row entry.
              firstUncheckedAij[i]++;
            }
          }
        }
      }

      else {

        // Copy the non-zero entries of the k-th column (below the diagonal)
        // of the system Matrix into the respective positions in the
        // auxilliary full vector.
        UInt i, cidx;
        for ( i = k + 1; i < sysMatDim_; i++ ) {

          // Avoid "overflow" into next row
          if ( firstUncheckedAij[i] < (UInt)rptrA[i+1] ) {

            cidx = cidxA[ firstUncheckedAij[i] ];

            if ( cidx == k ) {

              // This entry must be copied into the full auxilliary vector
              vecW[ i ] = entryA[ firstUncheckedAij[i] ];
              vecWLevel_[ i ] = 0;

              // This entry is a candidate for the largest row index
              lowest = i;

              // Increment pointer to next smallest unchecked row entry.
              firstUncheckedAij[i]++;
            }
          }
        }
      }
    }

    // ---------------------------------------
    // CASE 2: CRS_Matrix is in LEX sub-format
    // ---------------------------------------
    else if ( sysMat.GetCurrentLayout() == CRS_Matrix<T>::LEX ) {

      // Initialise auxilliary index vector
      if ( k == 0 ) {
        firstUncheckedAij.clear();
        firstUncheckedAij.reserve( sysMatDim_ + 1 );
        firstUncheckedAij.push_back( 0 );
        //        firstUncheckedAij.push_back( 0 );
        for ( UInt i = 1; i < sysMatDim_; i++ ) {
          firstUncheckedAij.push_back( rptrA[i] );
        }
      }

      if ( compFillLevels_ == false ) {

        // Copy the non-zero entries of the k-th column (below the diagonal)
        // of the system Matrix into the respective positions in the
        // auxilliary full vector.
        UInt i, cidx;
        for ( i = k + 1; i < sysMatDim_; i++ ) {

          // No need to search right of diagonal
          if ( firstUncheckedAij[i] < (UInt)diagA[i] ) {

            cidx = cidxA[ firstUncheckedAij[i] ];

            if ( cidx == k ) {

              // This entry must be copied into the full auxilliary vector
              vecW[ i ] = entryA[ firstUncheckedAij[i] ];              

              // This entry is a candidate for the largest row index
              lowest = i;

              // Increment pointer to next smallest unchecked row entry.
              firstUncheckedAij[i]++;
            }
          }
        }
      }

      else {

        // Copy the non-zero entries of the k-th column (below the diagonal)
        // of the system Matrix into the respective positions in the
        // auxilliary full vector.
        UInt i, cidx;
        for ( i = k + 1; i < sysMatDim_; i++ ) {

          // No need to search right of diagonal
          if ( firstUncheckedAij[i] < (UInt)diagA[i] ) {

            cidx = cidxA[ firstUncheckedAij[i] ];

            if ( cidx == k ) {

              // This entry must be copied into the full auxilliary vector
              vecW[ i ] = entryA[ firstUncheckedAij[i] ];
              vecWLevel_[ i ] = 0;

              // This entry is a candidate for the largest row index
              lowest = i;

              // Increment pointer to next smallest unchecked row entry.
              firstUncheckedAij[i]++;
            }
          }
        }
      }
    }

    // ----------------------------------------------
    // CASE 3: CRS_Matrix is in some other sub-format
    // ----------------------------------------------
    else {

      EXCEPTION( "Got a problem in CroutLU<T>::GetColOfA. "
               << "The system matrix is neither in LEX nor LEX_DIAG_FIRST "
               << "sub-format of CRS storage. Cannot re-order myself, "
               << "however, due to 'const' qualifier!@$%&" );

      // In this case we sort the matrix to LEX and solve
      // problem by recursion
      // sysMat.ChangeLayout( CRS_Matrix<T>::LEX );
      // lowest = GetColOfA( sysMat, vecW, k );
    }

    // return largest row index
    return lowest;
  }


  // ============================
  //   Export ILU factorisation
  // ============================
  template <typename T>
  void CroutLU<T>::ExportILUFactorisation( const char *fname ) {


    UInt i, j;

    // ********************
    //   Open output file
    // ********************
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      EXCEPTION( "CroutLU::ExportILUFactorisation: Cannot open file "
               << fname << " for writing!" );
    }

    // *********************
    //   Write file header
    // *********************

    // Matrix Market Format Specification
    if ( EntryType<T>::M_EntryType == BaseMatrix::DOUBLE ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
    }
    else if ( EntryType<T>::M_EntryType == BaseMatrix::COMPLEX ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate complex general\n" );
    }
    else {
      EXCEPTION( "CroutLU::ExportILUFactorisation: Can't identify "
               << "template parameter" );
    }

    // Comment
    fprintf( fp, "%%\n%% (I)LU factorisation matrix F = L + U - I " );
    fprintf( fp, "computed by CroutLU\n%%\n" );

    // Information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", sysMatDim_, sysMatDim_,
             rptrU_[sysMatDim_] + cptrL_[sysMatDim_] -2);

    // **********************
    //   Write entries of U
    // **********************

    // loop over all rows
    for ( i = 0; i < sysMatDim_; i++ ) {

      // First (diagonal) entry must be inverted

      // store row and column index
      fprintf( fp, "%6d\t%6d\t", i+1, i+1 );

      // store non-zero entry
      OpType<T>::ExportEntry( 1.0 / entryU_[ rptrU_[i] ], 0, 0, fp );
      fprintf( fp, "\n" );

      // loop over remaining entries in this row
      for ( j = rptrU_[i] + 1; j < rptrU_[i+1]; j++ ) {

        // store row and column index
        fprintf( fp, "%6d\t%6d\t", i+1, cidxU_[j]+1 );

        // store non-zero entry
        OpType<T>::ExportEntry( entryU_[j], 0, 0, fp );
        fprintf( fp, "\n" );
      }
    }

    // **********************
    //   Write entries of L
    // **********************

    // loop over all columns
    for ( i = 0; i < sysMatDim_; i++ ) {

      // loop over all entries in this column
      for ( j = cptrL_[i]; j < cptrL_[i+1]; j++ ) {

        // store row and column index
        fprintf( fp, "%6d\t%6d\t", ridxL_[j]+1, i+1 );

        // store non-zero entry
        OpType<T>::ExportEntry( entryL_[j], 0, 0, fp );
        fprintf( fp, "\n" );
      }
    }


    // *********************
    //   Close output file
    // *********************
    if ( fclose( fp ) == EOF ) {
      WARN("CroutLU::ExportILUFactorisation: Could not close file "
           << fname << " after writing!");
    }
  }

  // Explicit template instantiation
  template class CroutLU<Double>;
  template class CroutLU<Complex>;

}

#endif
