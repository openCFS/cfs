// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include <vector>

#include "MatVec/scrs_matrix.hh"
#include "MatVec/opdefs.hh"

#include "OLAS/precond/ILDLPrecond/ildlcnfactoriser.hh"
#include "OLAS/precond/ILDLPrecond/ildlprecond.hh"

// Used during development phase
// #define DEBUG_ILDLCNFACTORISER
#define ILDLCNFACTORISER_SAFEGUARD

namespace CoupledField {


  // ***********************
  //   Default constructor
  // ***********************
  template <class T>
  ILDLCNFactoriser<T>::ILDLCNFactoriser() {
    EXCEPTION( "Default constructor of ILDLCNFactoriser call was called! "
             << "This constructor is forbidden!" );
  }


  // ************************
  //   Standard constructor
  // ************************
  template <class T>
  ILDLCNFactoriser<T>::ILDLCNFactoriser( ParamNode *solverNode,
                                         InfoNode *olasInfo ) {


    // Currently this approach is not in a functional state
    WARN("The ILDLCNFactoriser is still in an experimental state. "
         << "Do not use is for production runs!");

    // Set pointers to communication objects
    this->xml_ = solverNode;
    this->olasInfo_ = olasInfo;

    // Initialise remaining attributes
    this->amFactorised_ = false;
    this->sysMatDim_    = 0;

  }


  // **********************
  //   Default Destructor
  // **********************
  template <class T>
  ILDLCNFactoriser<T>::~ILDLCNFactoriser() {


  }


  // *************
  //   Factorise
  // *************
  template <class T>
  void ILDLCNFactoriser<T>::Factorise( SCRS_Matrix<T> &sysMat,
                                       std::vector<T> &dataD,
                                       std::vector<UInt> &rptrU,
                                       std::vector<UInt> &cidxU,
                                       std::vector<T> &dataU,
                                       bool newPattern ) {


    UInt i, j, k, current, listPrevElem, listElem, numOffD;
    T elim, aux;


    // ======================
    //  Determine parameters
    // ======================

    // Problem size
    this->sysMatDim_ = sysMat.GetNumCols();

    // Dropping threshold
    Double tau = 0.01;
    this->xml_->Get("ILDLCN", "threshold", tau, false);

    // if ( tau <= 0.0 || tau > 1.0 ) {
    if ( tau < 0.0 ) {
      EXCEPTION( "ILDLCNFactoriser::Factorise: Dropping threshold tau = "
          << tau << ", but should be in (0,1]. Check setting of "
          << "ILDLPRECOND_tau parameter" );
    }
 
    // Shall we be verbose?
    bool logging = false;

    // =================
    //  Report start-up
    // =================
    if ( logging ) {
      (*cla) << " + Using ILDLCN variant with tau = " << tau << "\n"
             << std::endl;
    }

#ifdef DEBUG_ILDLCNFACTORISER
    (*debug) << " ILDLCNFACTORISER:\n FACTORISE\n" << std::endl;
#endif


    // ===================================
    //  Allocate memory for factorisation
    // ===================================

    // Make sure structure is clear
    dataD.clear();
    rptrU.clear();
    cidxU.clear();
    dataU.clear();

    // We use as initial estimate for the number of entries in U the number
    // of entries in the upper triangle of A; if this is not sufficient, the
    // STL will probably first double the vector size involving one copy step.
    // We hope that this then will suffice.
    UInt memSize;
    memSize = (UInt)floor( (double)(sysMat.GetNnz() - this->sysMatDim_) / 2 + 0.5 );
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

#ifdef DEBUG_ILDLCNFACTORISER
    (*debug) << " (Pre-)Allocated memory for storing factorisation"
             << " using memSize = " << memSize << '\n' << std::endl;
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
    const UInt *cidxA = sysMat.GetColPointer();
    const UInt *rptrA = sysMat.GetRowPointer();
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
    UInt* scanList_   = new UInt[this->sysMatDim_ + 1];
    UInt* activeList_ = new UInt[this->sysMatDim_ + 1];
    UInt* listIDX_    = new UInt[this->sysMatDim_ + 1];
    T*    listVAL_    = new   T [this->sysMatDim_ + 1];
    ASSERTMEM( scanList_  , sizeof( UInt ) * (this->sysMatDim_ + 1) );
    ASSERTMEM( activeList_, sizeof( UInt ) * (this->sysMatDim_ + 1) );
    ASSERTMEM( listIDX_   , sizeof( UInt ) * (this->sysMatDim_ + 1) );
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
    }
    scanList_[0]   = listEnd;
    activeList_[0] = listEnd;
    listIDX_[0]    = listEnd;
    listVAL_[0]    = listEnd;


    // =============================
    //  Prepare row norm estimation
    // =============================
    T *xi, *nu;
    T xiPlus, xiMinus;
    NEWARRAY( xi, T, this->sysMatDim_ );
    NEWARRAY( nu, T, this->sysMatDim_ );
    xi[1] = 1.0;
    for ( i = 1; i <= this->sysMatDim_; i++ ) {
      nu[i] = 0.0;
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


    // ============
    //  Main Loop
    // ============
    for ( k = 1; k <= this->sysMatDim_; k++ ) {

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
        (*debug) << " + updated row by row " << activeListElem
                 << std::endl;
#endif

        // Proceed to next row
        activeListElem = activeList_[activeListElem];
      }

#ifdef DEBUG_ILDLCNFACTORISER
      (*debug) << " + row pattern contains entries:";
      listElem = listIDX_[0];
      while ( listElem != listEnd ) {
        (*debug) << " " << listElem;
        listElem = listIDX_[listElem];
      }
      (*debug) << std::endl;
#endif


      // -------------------
      //  Estimate row norm
      // -------------------
      xiPlus  = +1.0 - nu[k];
      xiMinus = -1.0 - nu[k];
      xi[k] = Abs( xiPlus ) > Abs( xiMinus ) ? xiPlus : xiMinus;
      listElem = listIDX_[0];
      while ( listElem != listEnd ) {
        nu[listElem] += xi[k] * listVAL_[listElem];
        listElem = listIDX_[listElem];
      }

#ifdef DEBUG_ILDLCNFACTORISER
        (*debug) << " + norm estimate = " << Abs(xi[k]) << std::endl;
#endif


      // ----------------
      //  DO DA DROPPING
      // ----------------
      DropEntries( listIDX_, listVAL_, tau, Abs(xi[k]) );


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

#ifdef ILDLCNFACTORISER_SAFEGUARD
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

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
        (*debug) << " + scanListElem = " << scanListElem;
#endif

        // Check if column index of first entry is equal to k
        // (in this case the row must also be in the active list
        // and activeListElem = scanListElem)
        if ( cidxU[ firstU_[scanListElem] ] == k ) {

#ifdef DEBUG_ILDLCNFACTORISER
        (*debug) << ": incr. firstU_ from " << firstU_[scanListElem]
                 << " / " << cidxU[firstU_[scanListElem]]
                 << " --> ";
#endif

          // If this was the last row entry, delete row from scan list
          // and also from active list
          if ( firstU_[scanListElem] + 1 == rptrU[scanListElem+1] ) {

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
            (*debug) << firstU_[scanListElem]
                     << " / " << cidxU[firstU_[scanListElem]];
#endif

            // Check whether we must delete row from active list now
            if ( cidxU[ firstU_[scanListElem] ] > k + 1 ) {
              activeList_[activeListPrevElem] = activeList_[activeListElem];
              activeList_[activeListElem] = listNoElem;
              activeListElem = activeList_[activeListPrevElem];

#ifdef DEBUG_ILDLCNFACTORISER
              (*debug) << " (active list del)";
#endif
            }


            // If there is no other off-diagonal entry, then we can delete
            // the row from both lists
            else if ( firstU_[scanListElem] + 1 == rptrU[scanListElem+1] ) {

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
            if ( activeListPrevElem > scanListElem ) {
              EXCEPTION( "Error in active list add: (activeListPrevElem = "
                       << activeListPrevElem << ") < (scanListElem = "
                       << scanListElem << ")" );
            }
#endif

            activeList_[ scanListElem       ] = activeListElem;
            activeList_[ activeListPrevElem ] = scanListElem;
            activeListPrevElem = scanListElem;

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
        (*debug) << std::endl;
#endif

      }


#ifdef DEBUG_ILDLCNFACTORISER
      // We claim that activeListElem >= scanListElem always holds, so
      // now we should have activeListElem = listEnd
      if ( activeListElem != listEnd ) {
        EXCEPTION( "ILDLCNFACTORISER::Factorise: activeListElem = "
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

#ifdef DEBUG_ILDLCNFACTORISER
          (*debug) << " + added row " << k << " to scan list" << std::endl;
#endif

          // If the column index of the first off-diagonal entry is (k+1)
          // we must add it to the active list
          if ( cidxU[ firstU_[k] ] == (k+1) ) {

#ifdef DEBUG_ILDLCNFACTORISER
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

#ifdef DEBUG_ILDLCNFACTORISER
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
    delete [] ( firstU_ );  firstU_  = NULL;
    delete [] ( xi );  xi  = NULL;
    delete [] ( nu );  nu  = NULL;

#ifdef DEBUG_ILDLCNFACTORISER
    (*debug) << std::endl;
#endif

  }


  // ****************
  //   Drop Entries
  // ****************
  template <class T>
  void ILDLCNFactoriser<T>::DropEntries( UInt *listIDX, T *listVAL, Double tau,
                                         Double norm ) {


    std::cerr << " Row norm estimate = " << norm << std::endl;

    Double threshold = tau / norm;

    UInt listNoElem = 0;
    UInt listEnd = this->sysMatDim_ + 1;

    UInt listElem = listIDX[0];
    UInt listPrevElem = listNoElem;

#ifdef DEBUG_ILDLTPFACTORISER
    UInt numDropped = 0;
#endif

    while ( listElem != listEnd ) {

      if ( Abs(listVAL[listElem]) <= threshold ) {
        listIDX[listPrevElem] = listIDX[listElem];
        listIDX[listElem] = listNoElem;
        listVAL[listElem] = 0.0;
        listElem = listIDX[listPrevElem];

#ifdef DEBUG_ILDLTPFACTORISER
        numDropped++;
#endif

      }
      else {
        listElem = listIDX[listElem];
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

  }

  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class ILDLCNFactoriser<Double>;
    template class ILDLCNFactoriser<Complex>;
  #endif
  
}

