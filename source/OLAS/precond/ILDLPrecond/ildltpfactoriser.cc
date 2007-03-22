// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include <vector>
#include <algorithm>

#include "precond/ILDLPrecond/ildltpfactoriser.hh"
#include "precond/ILDLPrecond/ildlprecond.hh"

// Used during development phase
// #define DEBUG_ILDLTPFACTORISER
#define ILDLTPFACTORISER_SAFEGUARD

namespace OLAS {


  // ***********************
  //   Default constructor
  // ***********************
  template <class T>
  ILDLTPFactoriser<T>::ILDLTPFactoriser() {

    ENTER_FCN( "ILDLTPFactoriser::ILDLTPFactoriser" );

    (*error) << "Default constructor of ILDLTPFactoriser call was called! "
             << "This constructor is forbidden!";
    Error( __FILE__, __LINE__ );

  }


  // ************************
  //   Standard constructor
  // ************************
  template <class T>
  ILDLTPFactoriser<T>::ILDLTPFactoriser( OLAS_Params *myParams,
                                         OLAS_Report *myReport ) {

    ENTER_FCN( "ILDLTPFactoriser::ILDLTPFactoriser" );

    // Currently this approach is not in a functional state
    (*warning) << "The ILDLTPFactoriser is still in an experimental state. "
               << "Do not use is for production runs!";
    Warning( __FILE__, __LINE__ );

    // Set pointers to communication objects
    this->myParams_ = myParams;
    this->myReport_ = myReport;

    // Initialise remaining attributes
    this->amFactorised_ = false;
    this->sysMatDim_    = 0;

  }


  // **********************
  //   Default Destructor
  // **********************
  template <class T>
  ILDLTPFactoriser<T>::~ILDLTPFactoriser() {

    ENTER_FCN( "ILDLTPFactoriser::~ILDLTPFactoriser" );

  }


  // *************
  //   Factorise
  // *************
  template <class T>
  void ILDLTPFactoriser<T>::Factorise( SCRS_Matrix<T> &sysMat,
                                       std::vector<T> &dataD,
                                       std::vector<UInt> &rptrU,
                                       std::vector<UInt> &cidxU,
                                       std::vector<T> &dataU,
                                       bool newPattern ) {

    ENTER_FCN( "ILDLTPFactoriser::Factorise" );

    UInt i, j, k, current, listPrevElem, listElem, numOffD;
    T elim, aux;


    // ======================
    //  Determine parameters
    // ======================

    // Problem size
    this->sysMatDim_ = sysMat.GetNcols();

    // Dropping threshold
    Double tau = this->myParams_->GetDoubleValue( "ILDLPRECOND_tau" );
    if ( tau <= 0.0 || tau > 1.0 ) {
      (*error) << "ILDLTPFactoriser::Factorise: Dropping threshold tau = "
               << tau << ", but should be in (0,1]. Check setting of "
               << "ILDLPRECOND_tau parameter";
      Error( __FILE__, __LINE__ );
    }

    // Maximal number of entries in a row
    UInt fillVal = this->myParams_->GetIntValue( "ILDLPRECOND_fillVal" );
    if ( fillVal <= 0 ) {
      (*error) << "ILDLTPFactoriser::Factorise: fill value factor fillVal = "
               << fillVal << ", but should be positive integer. Check setting "
               << "of ILDLPRECOND_fillVal parameter";
      Error( __FILE__, __LINE__ );
    }

    UInt entriesPerRow = (sysMat.GetNnz() - this->sysMatDim_);
    entriesPerRow = (UInt) std::floor( (double)(entriesPerRow/(2*this->sysMatDim_))
                                       + 0.5 );
    UInt maxFill = entriesPerRow * fillVal;

    // Shall we be verbose?
    bool logging = this->myParams_->GetIntValue( "ILDLPRECOND_logging" ) > 0;


    // =================
    //  Report start-up
    // =================
    if ( logging ) {
      (*cla) << " + Using ILDL(tau,p) variant with (tau,p) = (" << tau << ", "
             << fillVal << ")\n"
             << " + average number of entries per row (in upper triangle): "
             << entriesPerRow << '\n'
             << " + bound on number of entries per row: " << maxFill
             << '\n' << std::endl;
    }

#ifdef DEBUG_ILDLTPFACTORISER
    (*debug) << " ILDLTPFACTORISER:\n FACTORISE\n" << std::endl;
#endif


    // ===================================
    //  Allocate memory for factorisation
    // ===================================

    // Make sure structure is clear
    dataD.clear();
    rptrU.clear();
    cidxU.clear();
    dataU.clear();

    // Determine maximal number of non-zero entries in factor U = L^T.
    // Avoid the unlikely case of over-allocation and make room for one
    // additional entry (one-based indexing). 
    UInt memSize = maxFill * this->sysMatDim_;
    if ( memSize > (sysMat.GetNnz() - this->sysMatDim_) / 2 ) {
      memSize = (sysMat.GetNnz() - this->sysMatDim_) / 2;
    }
    memSize++;

    // Using STL vectors for the factors, we do a pre-allocation here,
    // (adding one entry for one-based indexing)
    dataD.reserve( this->sysMatDim_ + 1 );
    rptrU.reserve( this->sysMatDim_ + 2 );
    cidxU.reserve( memSize );
    dataU.reserve( memSize );

    // OLAS uses one-base indexing, so we push_back a zero up front
    dataD.push_back( 0 );
    rptrU.push_back( 0 );
    cidxU.push_back( 0 );
    dataU.push_back( 0 );

#ifdef DEBUG_ILDLTPFACTORISER
    (*debug) << " (Pre-)Allocated memory for storing factorisation"
             << " using maxFill = " << maxFill << '\n' << std::endl;
#endif


    // =======================
    //  Prepare status report
    // =======================

    // Needed for writing progress report of factorisation
    UInt percentDone = 0;
    Double actDone = 0.0;
    (*cla) << " Factorisation done:\n" << " 0%" << std::flush;


    // =================
    //  Get matrix data
    // =================
    const Integer *cidxA = sysMat.GetColPointer();
    const Integer *rptrA = sysMat.GetRowPointer();
    const T *dataA = sysMat.GetDataPointer();


    // =======================
    //  Init diagonal factor
    // =======================
    for ( i = 1; i <= this->sysMatDim_; i++ ) {
      // dataD[i] = dataA[ rptrA[i] ];
      dataD.push_back( dataA[ rptrA[i] ] );
    }


    // ======================
    //  Prepare linked lists
    // ======================

    // Allocate memory for linked lists
    // (in future releases these will be attributes and (de-)allocation
    // will be handled by constructor/destructor)
    UInt* scanList_   = New UInt[this->sysMatDim_ + 1];
    UInt* activeList_ = New UInt[this->sysMatDim_ + 1];
    UInt* listIDX_    = New UInt[this->sysMatDim_ + 1];
    T*    listVAL_    = New   T [this->sysMatDim_ + 1];
    AssertMem( scanList_  , sizeof( UInt ) * (this->sysMatDim_ + 1) );
    AssertMem( activeList_, sizeof( UInt ) * (this->sysMatDim_ + 1) );
    AssertMem( listIDX_   , sizeof( UInt ) * (this->sysMatDim_ + 1) );
    AssertMem( listVAL_   , sizeof(   T  ) * (this->sysMatDim_ + 1) );

    // Generate markers for linked lists
    UInt scanListElem, scanListPrevElem;
    UInt activeListElem, activeListPrevElem;
    UInt listEnd = this->sysMatDim_ + 1;
    UInt listNoElem = 0;

    // Initialise linked lists to be empty
    for ( UInt j = 1; j <= this->sysMatDim_; j++ ) {
      scanList_[j]   = listNoElem;
      activeList_[j] = listNoElem;
      listIDX_[j]    = listNoElem;
      listVAL_[j]    = listNoElem;
    }
    scanList_[0]   = listEnd;
    activeList_[0] = listEnd;
    listIDX_[0]    = listEnd;
    listVAL_[0]    = listEnd;


    // =========
    //  Startup
    // =========

    // Generate array for keeping track of first entry in a rwo that
    // has a column index larger or equal to the current row index
    UInt *firstU_;
    NewArray( firstU_, UInt, this->sysMatDim_ );

    // First row starts with index 1 in CRS data array
    rptrU.push_back(1);


    // ============
    //  Main Loop
    // ============
    for ( k = 1; k <= this->sysMatDim_; k++ ) {

#ifdef DEBUG_ILDLTPFACTORISER
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
      actDone = (double)(k*100) / (double)this->sysMatDim_;
      actDone = (UInt)(actDone/10.0)*10;
      if ( actDone > percentDone ) {
        percentDone = (UInt)actDone;
        (*cla) << " .. " << percentDone << "%" << std::flush;
      }

      // Insert row k of A into linked list, but omit the diagonal entry
      if ( rptrA[k+1] - rptrA[k] > 1 ) {

        listIDX_[ 0 ]           = cidxA[ rptrA[k]+1 ];
        listVAL_[ listIDX_[0] ] = dataA[ rptrA[k]+1 ];

        for ( i = (UInt)rptrA[k] + 2; i < (UInt)rptrA[k+1]; i++ ) {
          listIDX_[ cidxA[i-1] ] = cidxA[i];
          listVAL_[ cidxA[ i ] ] = dataA[i];
        }
        listIDX_[ cidxA[i-1] ] = listEnd;
      }

#ifdef DEBUG_ILDLTPFACTORISER
      (*debug) << " + pattern of row in A contains entries:";
      listElem = listIDX_[0];
      while ( listElem != listEnd ) {
        (*debug) << " " << listElem;
        listElem = listIDX_[listElem];
      }
      (*debug) << std::endl;
#endif

      // All rows in active list contribute to current row
      activeListElem = activeList_[0];
      while ( activeListElem != listEnd ) {

        // Get data array index of row entry for column k
        i = firstU_[ activeListElem ];

        // Get value of row entry for column k
        elim = dataU[i];

        // Scale value with corresponding entry in diagonal factor
        elim *= dataD[ activeListElem ];


        // -------------------------------------------
        //  Add multiple of active row to current row
        // -------------------------------------------
        listPrevElem = listNoElem;
        listElem = listIDX_[0];
        for ( j = i + 1; j < rptrU[activeListElem+1]; j++ ) {

          // get current column index
          current = cidxU[j];

          // walk through linked list until column index becomes
          // larger or equal to current column index in active row
          while ( listElem < current ) {
            listPrevElem = listElem;
            listElem = listIDX_[listElem];
          }

          // both rows have an entry with current column index
          if ( listElem == current ) {

            // do numerical operation to update entry
            listVAL_[current] -= elim * dataU[j];

          }

          // the update row has an extra entry, so fill-in occurs
          else {

            // do numerical operation and generate fill-in
            listVAL_[current] = - elim * dataU[j];

            // insert fill-in into list management structure
            listIDX_[listPrevElem] = current;
            listPrevElem = current;
            listIDX_[current] = listElem;
          }
        }

#ifdef DEBUG_ILDLTPFACTORISER
        (*debug) << " + updated row by row " << activeListElem
                 << std::endl;
#endif

        // Proceed to next row
        activeListElem = activeList_[activeListElem];
      }

#ifdef DEBUG_ILDLTPFACTORISER
      (*debug) << " + row pattern contains entries:";
      listElem = listIDX_[0];
      while ( listElem != listEnd ) {
        (*debug) << " " << listElem;
        listElem = listIDX_[listElem];
      }
      (*debug) << std::endl;
#endif

      // ----------------
      //  DO DA DROPPING
      // ----------------
      DropEntries( listIDX_, listVAL_, tau, maxFill );

      // ---------------------
      //  Convert list to CRS
      // ---------------------

      // Now we must insert the column index and value information
      // from the linked list into the data structure for U = L^T

      numOffD = 0;
      listElem = listIDX_[0];
      listIDX_[0] = listEnd;

      // Scan linked list
      while( listElem != listEnd ) {

        // Insert column index into matrix structure
        cidxU.push_back( listElem );

        // Insert numerical value into matrix structure
        // and scale them at the same time by d_kk
        aux = listVAL_[listElem] / dataD[k];
        dataU.push_back( aux );

        // Update entry in diagonal factor
        dataD[ listElem ] -= aux * dataD[k] * aux;

#ifdef ILDLTPFACTORISER_SAFEGUARD
        // Remove old entry from list and proceed to next list element
        listPrevElem           = listElem;
        listElem               = listIDX_[listElem];
        listIDX_[listPrevElem] = 0;
        listVAL_[listPrevElem] = 0.0;
#else
        // Proceed to next list element
        listElem = listIDX_[listElem];
#endif


        // Increment number of row entries
        numOffD++;
      }

#ifdef DEBUG_ILDLTPFACTORISER
      (*debug) << " + converted row " << k << " to CRS: " << numOffD
               << " off-diagonal entries" << std::endl;
#endif

      // Generate pointer to start of next row
      rptrU.push_back( rptrU[k] + numOffD );


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
      while ( scanListElem != listEnd ) {

#ifdef DEBUG_ILDLTPFACTORISER
        (*debug) << " + scanListElem = " << scanListElem;
#endif

        // Check if column index of first entry is equal to k
        // (in this case the row must also be in the active list
        // and activeListElem = scanListElem)
        if ( cidxU[ firstU_[scanListElem] ] == k ) {

#ifdef DEBUG_ILDLTPFACTORISER
        (*debug) << ": incr. firstU_ from " << firstU_[scanListElem]
                 << " / " << cidxU[firstU_[scanListElem]]
                 << " --> ";
#endif

          // If this was the last row entry, delete row from scan list
          // and also from active list
          if ( firstU_[scanListElem] + 1 == rptrU[scanListElem+1] ) {

#ifdef DEBUG_ILDLTPFACTORISER
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
            firstU_[scanListElem]++;

#ifdef DEBUG_ILDLTPFACTORISER
            (*debug) << firstU_[scanListElem]
                     << " / " << cidxU[firstU_[scanListElem]];
#endif

            // Check whether we must delete row from active list now
            if ( cidxU[ firstU_[scanListElem] ] > k + 1 ) {
              activeList_[activeListPrevElem] = activeList_[activeListElem];
              activeList_[activeListElem] = listNoElem;
              activeListElem = activeList_[activeListPrevElem];

#ifdef DEBUG_ILDLTPFACTORISER
              (*debug) << " (active list del)";
#endif
            }


            // If there is no other off-diagonal entry, then we can delete
            // the row from both lists
            else if ( firstU_[scanListElem] + 1 == rptrU[scanListElem+1] ) {

#ifdef DEBUG_ILDLTPFACTORISER
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

#ifdef DEBUG_ILDLTPFACTORISER
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
        else if ( cidxU[ firstU_[scanListElem] ] == k + 1 ) {

          // If there is still another off-diagonal entry with larger column
          // index available add row to active list
          if ( firstU_[scanListElem] + 1 < rptrU[scanListElem+1] ) {

#ifdef DEBUG_ILDLTPFACTORISER
            if ( activeListPrevElem > scanListElem ) {
              (*error) << "Error in active list add: (activeListPrevElem = "
                       << activeListPrevElem << ") < (scanListElem = "
                       << scanListElem << ")";
              Error( __FILE__, __LINE__ );
            }
#endif

            activeList_[ scanListElem       ] = activeListElem;
            activeList_[ activeListPrevElem ] = scanListElem;
            activeListPrevElem = scanListElem;

#ifdef DEBUG_ILDLTPFACTORISER
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

#ifdef DEBUG_ILDLTPFACTORISER
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

#ifdef DEBUG_ILDLTPFACTORISER
        (*debug) << std::endl;
#endif

      }


#ifdef DEBUG_ILDLTPFACTORISER
      // We claim that activeListElem >= scanListElem always holds, so
      // now we should have activeListElem = listEnd
      if ( activeListElem != listEnd ) {
        (*error) << "ILDLTPFACTORISER::Factorise: activeListElem = "
                 << activeListElem << ", but should be " << listEnd;
        Error( __FILE__, __LINE__ );
      }
#endif
      

      // --------------------------------------------
      //  Add k-th row to auxilliary data structures
      // --------------------------------------------

      // Determine number off off-diagonal entries in row
      // numOffD = rptrU[k+1] - rptrU[k];

      // Make sure that row has any off-diagonal elements at all
      if ( numOffD > 0 ) {

        // Set pointer to first off-diagonal entry
        firstU_[k] = rptrU[k];

        // If there are two off-diagonal entries, at least one must
        // have a column index > (k+1), so we must add row to scan list
        if ( numOffD > 1 ) {

          // Addition is performed at the end (so the list will
          // always be ordered with increasing row index)
          // scanList_[k] = scanList_[0];
          // scanList_[0] = k;
          scanList_[scanListPrevElem] = k;
          scanList_[k] = listEnd;

#ifdef DEBUG_ILDLTPFACTORISER
          (*debug) << " + added row " << k << " to scan list" << std::endl;
#endif

          // If the column index of the first off-diagonal entry is (k+1)
          // we must add it to the active list
          if ( cidxU[ firstU_[k] ] == (k+1) ) {

#ifdef DEBUG_ILDLTPFACTORISER
            (*debug) << " + added row " << k << " to active list" << std::endl;
#endif
            // Addition is performed at the end
            // activeList_[k] = activeList_[0];
            // activeList_[0] = k;
            activeList_[activeListPrevElem] = k;
            activeList_[k] = listEnd;
          }
        }
      }

#ifdef DEBUG_ILDLTPFACTORISER
      (*debug) << " --> scanList_:";
      for ( UInt ii = 0; ii <= this->sysMatDim_; ii++ ) {
        (*debug) << " " << scanList_[ii];
      }
      (*debug) << std::endl;
      (*debug) << " --> activeList_:";
      for ( UInt ii = 0; ii <= this->sysMatDim_; ii++ ) {
        (*debug) << " " << activeList_[ii];
      }
      (*debug) << std::endl;
#endif
    }


    // ============
    //   Clean-up
    // ============

    // Finish logging to las-file
    if ( percentDone < 100 ) {
      (*cla) << " .. 100%";
    }
    (*cla) << std::endl;

    delete[] scanList_   ;
    delete[] activeList_ ;
    delete[] listIDX_    ; 
    delete[] listVAL_    ;
    DeleteArray( firstU_ );

#ifdef DEBUG_ILDLTPFACTORISER
    (*debug) << std::endl;
#endif

  }


  // ****************
  //   Drop Entries
  // ****************
  template <class T>
  void ILDLTPFactoriser<T>::DropEntries( UInt *listIDX, T *listVAL, Double tau,
                                         UInt maxFill ) {

    ENTER_FCN( "ILDLTPFactoriser::DropEntries" );

    UInt j;
    UInt listEnd = this->sysMatDim_ + 1;
    UInt listNoElem = 0;
    UInt listElem, listPrevElem;
    UInt listLength = 0;
    Double norm = 0.0;
    Double threshold = 0.0;

    std::vector<double> filter;

    // =========================
    //  Determine 1-norm of row
    // =========================
    listElem = listIDX[0];
    while ( listElem != listEnd ) {
      norm += Abs( listVAL[listElem] );
      listElem = listIDX[listElem];
    }

    // =====================================
    //  Drop entries smaller than threshold
    //  and copy other ones into STL vector
    // =====================================
    threshold = tau * norm;

#ifdef DEBUG_ILDLTPFACTORISER
    (*debug) << " + 1-norm of row = " << norm << std::endl;
    (*debug) << " + threshold = " << threshold << std::endl;
    UInt numDropped = 0;
#endif

    listElem = listIDX[0];
    listPrevElem = listNoElem;

    while ( listElem != listEnd ) {

      // Entry too small, remove from linked list
      if ( Abs(listVAL[listElem]) < threshold ) {
        listIDX[listPrevElem] = listIDX[listElem];
        listIDX[listElem] = listNoElem;
        listVAL[listElem] = 0.0;
        listElem = listIDX[listPrevElem];

#ifdef DEBUG_ILDLTPFACTORISER
        numDropped++;
#endif

      }

      // Determine modulus of entry entry and insert
      // it into STL vector for later sorting
      else {
        filter.push_back( Abs( listVAL[listElem] ) );
        listPrevElem = listElem;
        listElem = listIDX[listElem];
        listLength++;
      }
    }

#ifdef DEBUG_ILDLTPFACTORISER
    (*debug) << " + dropped " << numDropped << " entries\n";
    (*debug) << " + row pattern after dropping:";
    listElem = listIDX[0];
    while ( listElem != listEnd ) {
      (*debug) << " " << listElem;
      listElem = listIDX[listElem];
    }
    (*debug) << std::endl;
#endif

    // Check number of entries left in row
    if ( listLength > maxFill ) {

#ifdef DEBUG_ILDLTPFACTORISER
      (*debug) << " + searching " << maxFill << " largest entries"
               << std::endl;
#endif

      // ==================================
      //  Determine threshold for keeping
      //  only the maxFill largest entries
      // ==================================
      std::vector<double>::iterator cutOff = filter.begin() + maxFill - 1;
      nth_element( filter.begin(), cutOff, filter.end(),
                   std::greater<double>());
      // partial_sort( filter.begin(), cutOff, filter.end(),
      //              std::greater<double>());
      threshold = *cutOff;

      // =======================================================
      //  Walk through list and drop entries that are too small
      // =======================================================
      UInt numTie = 0;
      listElem = listIDX[0];
      listPrevElem = listNoElem;

      while ( listElem != listEnd ) {
        if ( Abs(listVAL[listElem]) < threshold ) {
          listIDX[listPrevElem] = listIDX[listElem];
          listIDX[listElem] = listNoElem;
          listVAL[listElem] = 0.0;
          listElem = listIDX[listPrevElem];
          listLength--;
        }
        else {

          // Count number of items with threshold value
          if ( Abs(listVAL[listElem]) == threshold ) {
            numTie++;
          }
          listPrevElem = listElem;
          listElem = listIDX[listElem];
        }
      }

      // ===============
      //  Test for ties
      // ===============
      Integer excess = listLength - maxFill;
      UInt dropTie;
      if ( excess > 0 ) {

        // Since there is a tie, there must be several elements with the
        // value threshold in the list. Remove the last excess occurences
        if ( excess > numTie ) {
          (*error) << "Internal error. My master cannot program :(";
          Error( __FILE__, __LINE__ );
        }
        dropTie = numTie - excess;

//         std::cerr << "\n\n TIE detected: All men to battle stations!\n"
//                   << " numTie = " << numTie
//                   << " , listLength = " << listLength << " , maxFill = "
//                   << maxFill << " , excess = " << excess << " , dropTie = "
//                   << dropTie
//                   << std::endl;

        numTie = 0;
        listElem = listIDX[0];
        listPrevElem = listNoElem;
        while ( listElem != listEnd ) {
          if ( Abs(listVAL[listElem]) == threshold ) {
            numTie++;

//            std::cerr << " Detected potential dropper (numTie = "
//                      << numTie << " )" << std::endl;

            if ( numTie > dropTie ) {
              listIDX[listPrevElem] = listIDX[listElem];
              listIDX[listElem] = listNoElem;
              listVAL[listElem] = 0.0;
              listElem = listIDX[listPrevElem];
              listLength--;
            }
            else {
              listPrevElem = listElem;
              listElem = listIDX[listElem];
            }
          }
          else {
            listPrevElem = listElem;
            listElem = listIDX[listElem];
          }
        }

        excess = listLength - maxFill;
        if ( excess > 0 ) {

          std::cerr << "\n\n TIE still there!\n"
                  << "listLength = " << listLength << std::endl;
          (*error) << "FUBAR!!!";
          Error( __FILE__, __LINE__ );
        }
      }


#ifdef DEBUG_ILDLTPFACTORISER
      (*debug) << " + row pattern after final dropping:";
      listElem = listIDX[0];
      while ( listElem != listEnd ) {
        (*debug) << " " << listElem;
        listElem = listIDX[listElem];
      }
      (*debug) << std::endl;
#endif

    }
  }
}

