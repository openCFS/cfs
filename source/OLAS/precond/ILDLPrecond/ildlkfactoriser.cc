// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <vector>
#include <cstdio>

#include "MatVec/scrs_matrix.hh"

#include "OLAS/precond/ILDLPrecond/ildlkfactoriser.hh"
#include "OLAS/precond/ILDLPrecond/ildlprecond.hh"


namespace CoupledField {


  // ***********************
  //   Default constructor
  // ***********************
  template <class T>
  ILDLKFactoriser<T>::ILDLKFactoriser() {
    EXCEPTION( "Default constructor of ILDLKFactoriser call was called! "
             << "This constructor is forbidden!" );
  }


  // ************************
  //   Standard constructor
  // ************************
  template <class T>
  ILDLKFactoriser<T>::ILDLKFactoriser( PtrParamNode solverNode,
                                       PtrParamNode olasInfo ) {


    // Set pointers to communication objects
    this->xml_ = solverNode;
    this->olasInfo_ = olasInfo;

    // Initialise remaining attributes
    this->amFactorised_ = false;
    this->sysMatDim_    = 0;
    level_        = 0;
  }


  // **********************
  //   Default Destructor
  // **********************
  template <class T>
  ILDLKFactoriser<T>::~ILDLKFactoriser() {


  }


  // *************
  //   Factorise
  // *************
  template <class T>
  void ILDLKFactoriser<T>::Factorise( SCRS_Matrix<T> &sysMat,
                                      std::vector<T> &dataD,
                                      std::vector<UInt> &rptrU,
                                      std::vector<UInt> &cidxU,
                                      std::vector<T> &dataU,
                                      bool newPattern ) {


    // Get new problem size and perform consistency check
    if ( this->amFactorised_ == false ) {
      this->sysMatDim_ = sysMat.GetNumCols();
    }
    else {
      if ( newPattern == false &&
           this->sysMatDim_ != sysMat.GetNumCols() ) {
        EXCEPTION( "ILDLKFactoriser::Factorise: Now I got you. "
                 << "newPattern = false, but matrix dimension changed "
                 << "from " << this->sysMatDim_ << " to "
                 << sysMat.GetNumCols() );
      }
      else {
        this->sysMatDim_ = sysMat.GetNumCols();
        this->amFactorised_ = false;
      }
    }

    // What fill level is allowed
    level_ = 1;
    this->xml_->Get("ILDK",ParamNode::INSERT)->GetValue( "level", level_, ParamNode::INSERT);

    // Start pattern analysis combined with factorisation
    if ( this->amFactorised_ == false || newPattern == true ) {
      FactoriseAnalytic( sysMat, dataD, rptrU, cidxU, dataU );
    }

    // Start numerical factorisation
    else {
      // FactoriseNumerics( sysMat, dataD, rptrU, cidxU, dataU );
      FactoriseAnalytic( sysMat, dataD, rptrU, cidxU, dataU );
    }

  }


  // *********************
  //   FactoriseAnalytic
  // *********************
  template<typename T>
  void ILDLKFactoriser<T>::FactoriseAnalytic( SCRS_Matrix<T> &sysMat,
                                              std::vector<T> &dataD,
                                              std::vector<UInt> &rptrU,
                                              std::vector<UInt> &cidxU,
                                              std::vector<T> &dataU ) {


    UInt i, j, k, current, listPrevElem, listElem, numOffD;
    UInt lvlParent1, lvlParent2, auxLevel;
    T elim, aux;

    // Shall we be verbose?
    bool logging = false;


    // =================
    //  Report start-up
    // =================
    if ( logging ) {
      (*cla) << " + Using ILDL(k) variant with k = " << level_ << "\n\n"
	     << " Phase: Combined ANALYSE + FACTORISE" << std::endl;
    }

#ifdef DEBUG_ILDLKFACTORISER
    (*debug) << " ILDLKFACTORISER:\n Phase: ANALYSE + FACTORISE\n k = "
             << level_ << '\n' << std::endl;
#endif


    // =================
    //  Get matrix data
    // =================
    const UInt *cidxA = sysMat.GetColPointer();
    const UInt *rptrA = sysMat.GetRowPointer();
    const T *dataA = sysMat.GetDataPointer();


    // ======================
    //  Prepare linked lists
    // ======================

    // Allocate memory for linked lists
    // (in future releases these will be attributes and (de-)allocation
    // will be handled by constructor/destructor)
    UInt* scanList_   = new UInt[this->sysMatDim_ + 1];
    UInt* activeList_ = new UInt[this->sysMatDim_ + 1];
    UInt* listIDX_    = new UInt[this->sysMatDim_ + 1];
    UInt* listLVL_    = new UInt[this->sysMatDim_ + 1];
    T*    listVAL_    = new   T [this->sysMatDim_ + 1];
    ASSERTMEM( scanList_  , sizeof( UInt ) * (this->sysMatDim_ + 1) );
    ASSERTMEM( activeList_, sizeof( UInt ) * (this->sysMatDim_ + 1) );
    ASSERTMEM( listIDX_   , sizeof( UInt ) * (this->sysMatDim_ + 1) );
    ASSERTMEM( listLVL_   , sizeof( UInt ) * (this->sysMatDim_ + 1) );
    ASSERTMEM( listVAL_   , sizeof(   T  ) * (this->sysMatDim_ + 1) );

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
      listVAL_[j]    = listNoElem;
    }
    scanList_[0]   = listEnd;
    activeList_[0] = listEnd;
    listIDX_[0]    = listEnd;
    listLVL_[0]    = listEnd;
    listVAL_[0]    = listEnd;

#ifdef DEBUG_ILDLKFACTORISER
    (*debug) << " Initialised linked lists" << std::endl;
#endif


    // ==================================
    //  Reserve memory for factorisation
    // ==================================

    // Make sure structure is clear
    dataD.clear();
    rptrU.clear();
    cidxU.clear();
    dataU.clear();

    // We need an additional vector for storing the fill levels
    std::vector<UInt> fillU;

    // Estimate memory requirements
    UInt memSize = EstimateFactorMemory( sysMat, listIDX_ );

    // Using STL vectors for the factors, we do a pre-allocation here,
    // (adding one entry for one-based indexing)
    dataD.reserve( this->sysMatDim_ + 1 );
    rptrU.reserve( this->sysMatDim_ + 2 );
    cidxU.reserve( memSize );
    dataU.reserve( memSize );
    fillU.reserve( memSize );

    // OLAS uses one-base indexing, so we push_back a zero up front
    dataD.push_back( 0 );
    rptrU.push_back( 0 );
    cidxU.push_back( 0 );
    dataU.push_back( 0 );
    fillU.push_back( 0 );

#ifdef DEBUG_ILDLKFACTORISER
    (*debug) << " Allocated/Reserved memory for storing factorisation"
             << " using memSize = " << memSize << '\n' << std::endl;
#endif


    // =======================
    //  Init diagonal factor
    // =======================
    for ( i = 1; i <= this->sysMatDim_; i++ ) {
      dataD.push_back( dataA[ rptrA[i] ] );
    }


    // =========
    //  Startup
    // =========

    // Generate array for keeping track of first entry in a row that
    // has a column index larger or equal to the current row index
    UInt *firstU_;
    NEWARRAY( firstU_, UInt, this->sysMatDim_ );

    // First row starts with index 1 in CRS data array
    rptrU.push_back(1);


    // =======================
    //  Prepare status report
    // =======================

    // Needed for writing progress report of factorisation
    UInt percentDone = 0;
    Double actDone = 0.0;
    if ( logging == true ) {
      (*cla) << '\n';
    }
    (*cla) << " Factorisation done:\n" << " 0%" << std::flush;


    // ============
    //  Main Loop
    // ============
    for ( k = 1; k <= this->sysMatDim_; k++ ) {

#ifdef DEBUG_ILDLKFACTORISER
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

      // Insert row k of A into linked list, but omit the diagonal entry,
      // and set the initial fill-levels for its pattern
      if ( rptrA[k+1] - rptrA[k] > 1 ) {

        listIDX_[ 0 ]           = cidxA[ rptrA[k]+1 ];
        listVAL_[ listIDX_[0] ] = dataA[ rptrA[k]+1 ];
        listLVL_[ listIDX_[0] ] = 0;

        for ( i = (UInt)rptrA[k] + 2; i < (UInt)rptrA[k+1]; i++ ) {
          listIDX_[ cidxA[i-1] ] = cidxA[i];
          listVAL_[ cidxA[ i ] ] = dataA[i];
          listLVL_[ cidxA[ i ] ] = 0;
        }
        listIDX_[ cidxA[i-1] ] = listEnd;
      }

#ifdef DEBUG_ILDLKFACTORISER
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

        // Get fill level
        lvlParent1 = fillU[i];


        // -------------------------------------------
        //  Add multiple of active row to current row
        // -------------------------------------------
        listPrevElem = listNoElem;
        listElem = listIDX_[0];
        for ( j = i + 1; j < rptrU[activeListElem+1]; j++ ) {

          // get current column index
          current = cidxU[j];

          // get current fill level
          lvlParent2 = fillU[j];

          // compute new fill level for the two matrix entries causing
          // current update
          auxLevel = lvlParent1 + lvlParent2 + 1;

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

            // check, if fill level must be updated
            if ( auxLevel < listLVL_[current] ) {
              listLVL_[current] = auxLevel;
            }
          }

          // the update row has an extra entry, so fill-in occurs
          else {

            // do numerical operation and generate fill-in
            listVAL_[current] = - elim * dataU[j];

            // set fill level of fill-in entry
            listLVL_[current] = auxLevel;

            // insert fill-in into list management structure
            listIDX_[listPrevElem] = current;
            listPrevElem = current;
            listIDX_[current] = listElem;
          }
        }

#ifdef DEBUG_ILDLKFACTORISER
        (*debug) << " + updated row by row " << activeListElem
                 << std::endl;
#endif

        // Proceed to next row
        activeListElem = activeList_[activeListElem];
      }

#ifdef DEBUG_ILDLKFACTORISER
      (*debug) << " + row pattern contains entries:";
      listElem = listIDX_[0];
      while ( listElem != listEnd ) {
        (*debug) << " " << listElem;
        listElem = listIDX_[listElem];
      }
      (*debug) << std::endl;
      (*debug) << " + with fill levels:";
      listElem = listIDX_[0];
      while ( listElem != listEnd ) {
        (*debug) << " " << listLVL_[listElem];
        listElem = listIDX_[listElem];
      }
      (*debug) << std::endl;
#endif

      // ----------------
      //  DO DA DROPPING
      // ----------------
      DropEntries( listIDX_, listVAL_, listLVL_ );

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

        // Insert fill level into matrix structure
        fillU.push_back( listLVL_[listElem] );

        // Insert numerical value into matrix structure
        // and scale them at the same time by d_kk
        aux = listVAL_[listElem] / dataD[k];
        dataU.push_back( aux );

        // Update entry in diagonal factor
        dataD[ listElem ] -= aux * dataD[k] * aux;

#ifdef ILDLKFACTORISER_SAFEGUARD
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

#ifdef DEBUG_ILDLKFACTORISER
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

#ifdef DEBUG_ILDLKFACTORISER
        (*debug) << " + scanListElem = " << scanListElem;
#endif

        // Check if column index of first entry is equal to k
        // (in this case the row must also be in the active list
        // and activeListElem = scanListElem)
        if ( cidxU[ firstU_[scanListElem] ] == k ) {

#ifdef DEBUG_ILDLKFACTORISER
        (*debug) << ": incr. firstU_ from " << firstU_[scanListElem]
                 << " / " << cidxU[firstU_[scanListElem]]
                 << " --> ";
#endif

          // If this was the last row entry, delete row from scan list
          // and also from active list
          if ( firstU_[scanListElem] + 1 == rptrU[scanListElem+1] ) {

#ifdef DEBUG_ILDLKFACTORISER
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

#ifdef DEBUG_ILDLKFACTORISER
            (*debug) << firstU_[scanListElem]
                     << " / " << cidxU[firstU_[scanListElem]];
#endif

            // Check whether we must delete row from active list now
            if ( cidxU[ firstU_[scanListElem] ] > k + 1 ) {
              activeList_[activeListPrevElem] = activeList_[activeListElem];
              activeList_[activeListElem] = listNoElem;
              activeListElem = activeList_[activeListPrevElem];

#ifdef DEBUG_ILDLKFACTORISER
              (*debug) << " (active list del)";
#endif
            }


            // If there is no other off-diagonal entry, then we can delete
            // the row from both lists
            else if ( firstU_[scanListElem] + 1 == rptrU[scanListElem+1] ) {

#ifdef DEBUG_ILDLKFACTORISER
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

#ifdef DEBUG_ILDLKFACTORISER
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

#ifdef DEBUG_ILDLKFACTORISER
            if ( activeListPrevElem > scanListElem ) {
              EXCEPTION( "Error in active list add: (activeListPrevElem = "
                       << activeListPrevElem << ") < (scanListElem = "
                       << scanListElem << ")" );
            }
#endif

            activeList_[ scanListElem       ] = activeListElem;
            activeList_[ activeListPrevElem ] = scanListElem;
            activeListPrevElem = scanListElem;

#ifdef DEBUG_ILDLKFACTORISER
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

#ifdef DEBUG_ILDLKFACTORISER
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

#ifdef DEBUG_ILDLKFACTORISER
        (*debug) << std::endl;
#endif

      }


#ifdef DEBUG_ILDLKFACTORISER
      // We claim that activeListElem >= scanListElem always holds, so
      // now we should have activeListElem = listEnd
      if ( activeListElem != listEnd ) {
        EXCEPTION( "ILDLKFACTORISER::Factorise: activeListElem = "
                 << activeListElem << ", but should be " << listEnd );
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

#ifdef DEBUG_ILDLKFACTORISER
          (*debug) << " + added row " << k << " to scan list" << std::endl;
#endif

          // If the column index of the first off-diagonal entry is (k+1)
          // we must add it to the active list
          if ( cidxU[ firstU_[k] ] == (k+1) ) {

#ifdef DEBUG_ILDLKFACTORISER
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

#ifdef DEBUG_ILDLKFACTORISER
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


    // =====================
    //   Export Fill Levels
    // =====================
    std::string saveLevelsFile = "";
    this->xml_->Get("ILDLK", ParamNode::INSERT)
      ->GetValue( "saveLevelsFile", saveLevelsFile, ParamNode::INSERT);
    if( saveLevelsFile != "" ) {
      ExportFillLevels( saveLevelsFile.c_str(), rptrU, cidxU, fillU );
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
    delete [] ( firstU_ );  firstU_  = NULL;

#ifdef DEBUG_ILDLKFACTORISER
    (*debug) << std::endl;
#endif

  }


  // *********************
  //   FactoriseNumerics
  // *********************
  template<typename T>
  void ILDLKFactoriser<T>::FactoriseNumerics( SCRS_Matrix<T> &sysMat,
                                              std::vector<T> &dataD,
                                              std::vector<UInt> &rptrU,
                                              std::vector<UInt> &cidxU,
                                              std::vector<T> &dataU ) {


    // Shall we be verbose?
    // TODO: Check if this is still needed
    // bool logging = false;
  }


  // ************************
  //   EstimateFactorMemory
  // ************************
  template<typename T>
  UInt ILDLKFactoriser<T>::EstimateFactorMemory( SCRS_Matrix<T> &sysMat,
                                                 UInt* auxVec ) {


    // For large profiles the estimates will be way too large, or
    // maybe even negative, since unsigned int will be too short.
    // For the Gamm Seminar we used the value below instead. We must
    // develop a way to efficiently do the dynamic phase without an
    // a priori estimate!
    // return level_ * sysMat.GetNnz() - sysMatDim_;

    // NOTE: Changed profile computation. We now do it as in Sloan class,
    //       i.e. we use UInt and represent the large number by modulo
    //       operation w.r.t. UINT_MAX
    UInt i, j = 1;
    UInt k;
    Integer bw, bwLocal;
    UInt profile, profAux, profileMult;

    // Shall we be verbose?
    bool logging = false;

    // Get hold of column index array
    const UInt *cidxA = sysMat.GetColPointer();

    // Get hold of row pointer index array
    const UInt *rptrA = sysMat.GetRowPointer();

    // =====================================
    //  Determine the bandwidth and profile
    // =====================================

    bw = 0;

    // We use a bottom up approach in order to determine for each matrix
    // column the smallest row index of a non-zero entry in that column.
    for ( i = this->sysMatDim_; i > 0; i-- ) {

      // For profile
      for ( k = rptrA[i]; k < rptrA[i+1]; k++ ) {
        auxVec[ cidxA[k] ] = i;
      }

      // For bandwidth
      bwLocal = cidxA[ j - 1 ] - cidxA[ rptrA[i] ];
      if ( bwLocal > bw ) {
        bw = bwLocal;
      }
    }

    // Maximal fill-in will generate a column with all non-zero entries
    // between the diagonal and the row determined above. Summing these
    // differences up gives the profile
    // Also: clear auxilliary vector
    profile = 1;
    profileMult = 0;
    for ( i = 2; i <= this->sysMatDim_; i++ ) {
      profAux = profile;
      profile += i - auxVec[i] + 1;
      if ( profile < profAux ) {
        profileMult++;
      }
      auxVec[i] = 0;
    }

    // Represent profile as floating point number
    Double profileFP = profile + (double) std::numeric_limits<unsigned int>::max() * profileMult;

    // Report
    if ( logging ) {
      (*cla) << " Pattern analysis\n"
             << " + Matrix has a bandwidth of bw = " << bw
             << "\n + Full factor (L + D) would contain " << profileFP
             << " entries (at most)" << std::endl;
    }

    // ======================================================
    //  Compute estimation for memory requirements of factor
    // ======================================================

    Double fillIn, chunk, memSize;

    // Determine additional entries in profile
    fillIn = profile + (double) std::numeric_limits<unsigned int>::max() * profileMult - sysMat.GetNnz();

    // Divide into equivalent chunks
    chunk = fillIn / bw;

    // Assume that for each k we need one additional chunk
    // and add one chunk for safety
    memSize = chunk * ( level_ + 1 ) + sysMat.GetNnz() - this->sysMatDim_;

    // Check that we are not larger than profile
    memSize = memSize > profileFP ? profileFP : memSize;

    // Report
    if ( logging ) {
      (*cla) << " + Estimate of no. of entries in L: "
             << (UInt) memSize << std::endl;
    }

    // That's it
    return (UInt)memSize;
  }


  // ********************
  //   ExportFillLevels
  // ********************
  template <typename T>
  void ILDLKFactoriser<T>::ExportFillLevels( const char *fname,
                                             std::vector<UInt> rptrU,
                                             std::vector<UInt> cidxU,
                                             std::vector<UInt> fillU ) {


    UInt i, j;

    // ====================
    //   Open output file
    // ====================
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      EXCEPTION( "ILDLPrecond::ExportFillLevels: Cannot open file "
               << fname << " for writing!" );
    }

    // =====================
    //   Write file header
    // =====================

    // Matrix Market Format Specification
    fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );

    // Comment
    fprintf( fp, "%%\n%% Fill-levels of entries in matrix factor U = L^T\n" );

    // Information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", this->sysMatDim_, this->sysMatDim_,
             rptrU[this->sysMatDim_+1] - 1 );

    // =========================
    //   Write fill levels of U
    // =========================
    for ( i = 1; i <= this->sysMatDim_; i++ ) {
      for ( j = rptrU[i]; j < rptrU[i+1]; j++ ) {
        fprintf( fp, "%6d\t%6d\t%6d\n", i, cidxU[j], fillU[j] + 1 );
      }
    }

    // =====================
    //   Close output file
    // =====================
    if ( fclose( fp ) == EOF ) {
      WARN("ILDLPrecond::ExportFillLevels: Could not close file "
           << fname << " after writing!");
    }
  }


  // ****************
  //   Drop Entries
  // ****************
  template <class T>
  void ILDLKFactoriser<T>::DropEntries( UInt *listIDX, T *listVAL,
                                        UInt *listLVL ) {


    UInt listEnd = this->sysMatDim_ + 1;
    UInt listNoElem = 0;
    UInt listElem, listPrevElem;

    // ===============================
    //  Drop entries whose fill level
    //  exceeds the threshold
    // ===============================

#ifdef DEBUG_ILDLKFACTORISER
    UInt numDropped = 0;
#endif

    listElem = listIDX[0];
    listPrevElem = listNoElem;

    while ( listElem != listEnd ) {

      if ( listLVL[listElem] > level_ ) {
        listIDX[listPrevElem] = listIDX[listElem];
        listIDX[listElem] = listNoElem;
        listVAL[listElem] = 0.0;
        listLVL[listElem] = 0;
        listElem = listIDX[listPrevElem];

#ifdef DEBUG_ILDLKFACTORISER
        numDropped++;
#endif

      }

      else {
        listPrevElem = listElem;
        listElem = listIDX[listElem];
      }
    }

#ifdef DEBUG_ILDLKFACTORISER
    (*debug) << " + dropped " << numDropped << " entries\n";
    (*debug) << " + row pattern after dropping:";
    listElem = listIDX[0];
    while ( listElem != listEnd ) {
      (*debug) << " " << listElem;
      listElem = listIDX[listElem];
    }
    (*debug) << std::endl;
#endif

  }

  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class ILDLKFactoriser<Double>;
    template class ILDLKFactoriser<Complex>;
  #endif
  
}
