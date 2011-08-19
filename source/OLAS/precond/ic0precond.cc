// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/scrs_matrix.hh"

#include "OLAS/precond/ic0precond.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  template <typename T>
  IC0Precond<T>::IC0Precond( const StdMatrix& mat, PtrParamNode solverNode,
                             PtrParamNode olasInfo )
  {
    this->xml_ = solverNode;
    this->olasInfo_ = olasInfo;
    size_           = mat.GetNumRows();
    amFactorised_   = false;
  }


  // **************
  //   Destructor
  // **************
  template <typename T>
  IC0Precond<T>::~IC0Precond() {

    delete [] ( rptrU_ );
    delete [] ( cidxU_ );
    delete [] ( dataU_ );
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void IC0Precond<T>::Setup( SCRS_Matrix<T> &sysmat ) {

    UInt nnzA = (size_ + sysmat.GetNnz() ) / 2;

    if ( amFactorised_ ) {
      // IC0 has lready beeen performed; now we have to do it again, since
      // the entries of the system matrix have changed )e.g. nonlinear analysis)
      delete [] ( rptrU_ );  rptrU_  = NULL;
      delete [] ( cidxU_ );  cidxU_  = NULL;
      delete [] ( dataU_ );  dataU_  = NULL;
    }

    // Initialise internal data arrays
    NEWARRAY( rptrU_, UInt, size_+1 );
    NEWARRAY( cidxU_, UInt, nnzA );
    NEWARRAY( dataU_, T, nnzA );

    // =================
    //  Get matrix data
    // =================
    const UInt *cidxA = sysmat.GetColPointer();
    const UInt *rptrA = sysmat.GetRowPointer();
    const T *dataA = sysmat.GetDataPointer();

    //copy the structure
    for (UInt k=0; k<=size_; k++) {
      rptrU_[k] = rptrA[k];
    }

    for (UInt k=0; k<nnzA; k++) {
      cidxU_[k] = cidxA[k];
    }


    //compute maximal number of entries per row
    UInt maxRowEntries = 1;
    UInt rowEntries;
    for (UInt k=0; k<size_; k++) {
      rowEntries = rptrA[k+1] - rptrA[k];
      if ( rowEntries > maxRowEntries )
        maxRowEntries = rowEntries;
    }

    // some help variables
    T  Rkk, Rik, sqrRkk;
    UInt idx;

    T * Rk;
    NEWARRAY( Rk, T, size_ );

    UInt   nPD = 1;
    UInt   counterFailed = 1;
    UInt   percentDone = 0;
    Double diagScale = 1.0;
    Double actDone;
    UInt   jstart, jstop;

    while (nPD) {
      //set U = systemMatrix
      for (UInt k=0; k<nnzA; k++) {
        dataU_[k] = diagScale*dataA[k];
      }

      nPD = 0;

      //set zero
      for (UInt i=0; i<size_; i++) {
        Rk[i] = (T) 0;
      }

      // go over all equations
      for ( UInt k=0; k<size_; k++ ) {
        //	std::cout << "perform EQ: " << k << std::endl;

        // Keep user informed on progress
        actDone = (double)(k*100) / (double)(size_-1);
        actDone = (UInt)(actDone/10.0)*10;
        if ( actDone > percentDone ) {
          percentDone = (UInt)actDone;
          (*cla) << " .. " << percentDone << "%" << std::flush;
        }

        Rkk = dataU_[rptrA[k]];

        jstart = rptrA[k];
        jstop = rptrA[k+1];
        //copy the values of row k to Rk
        for ( UInt j=jstart; j<jstop; j++ ) {
          Rk[cidxA[j]] = dataU_[j];
        }

        //compute starting index
        UInt starti = 0;
        if (maxRowEntries - k  > 2) {
          starti = maxRowEntries - k -2;
        }

        for  ( UInt i=starti; i<k; i++ ) {
          Rik = (T) 0;
          Integer found = 0;
          for ( UInt j=rptrA[i]; j<rptrA[i+1]; j++) {
            if ( cidxA[j] == k) {
              Rik = dataU_[j];
              found = 1;
              break;;
            }
          }
          if ( found == 1) {
            Rkk = Rkk - Rik*Rik;
            for ( UInt j=rptrA[i]+1; j<rptrA[i+1]; j++ ) {
              idx = cidxA[j];
              if (idx > k) {
                Rk[idx] -= Rik * dataU_[j];
              }
            }
          }

        }

        if (  CheckPositivity(Rkk) > 0 ) {
          sqrRkk = sqrt(Rkk);
          dataU_[rptrA[k]] = sqrRkk;

          for (UInt i=rptrA[k]+1; i<rptrA[k+1]; i++) {
            dataU_[i] = Rk[cidxA[i]] / sqrRkk;
          }
        }
        else {
          // we have to left the loop, scale the diagonal and try again
          diagScale = 2.0 * counterFailed;
          counterFailed++;
          nPD = 1;
//           (*cla) << "\n Perform scaling of system matrix: " << counterFailed
//          << std::endl;
        }

        if (nPD ==1 )
          break;

        // set Rk to zero
        for ( UInt j=jstart; j<jstop; j++ ) {
          Rk[cidxA[j]] = 0;
        }

      } // end of equation loop

    } // end of while loop (to get postiv definite precond. matrix)

    delete [] ( Rk );  Rk  = NULL;
    amFactorised_ = true;
  }


/* -------------------- new version ------------------------- Fabian -----
  template<typename T>
  void IC0Precond<T>::Setup( SCRS_Matrix<T> &sysmat )
  {

    PROFILE( "IC0Precond::Setup",
             size_ * BlockSize<T>::size * BlockSize<T>::size );


    // store all l-values concurrently in this structure to have a faster access
    // A column is represented via a map where the key is the row index.
    vector<Col_Map>(size_+1) colCM;
    // result of a find of Col_Map
    Col_Map::iterator *searchCM;

    bool   spd = false; //  check for positive definiteness
    UInt   counterFailed = 1;
    UInt   percentDone = 0;
    Double diagScale = 1.0;
    Double actDone;

    while (!spd)
    {
      //set U = systemMatrix
      for (UInt k=1; k<=nnzA; k++)
      {
	    dataU_[k] = diagScale*dataA[k];
      }

      spd = true;

      // go over all equations
      for ( UInt k=1; k<=size_; k++ )
      {
         // Keep user informed on progress
         actDone = (double)(k*100) / (double)size_;
	     actDone = (UInt)(actDone/10.0)*10;
         if ( actDone > percentDone ) {
           percentDone = (UInt)actDone;
           (*cla) << " .. " << percentDone << "%" << std::flush;
         }

         // diagonal element
         T l_kk = dataU_[rptrA[k]];

         // compute t_kk via colum k = a_kk - sum l_jk * l_jk
         Col_Map& k_col = colCM[k];
         for (Col_Map::iterator iter = k_col.begin(); iter != k_col.end(); iter++) {
         {
            cout << "l[k=" << k << "][j=" << (*iter).first << "]=" << (*iter).second << endl;
            T val = (*iter).second;
            l_kk -= val * val;
         }

         // check if we are positive before taking the root
         if(!CheckPositivity(l_kk))
         {
           // we have to leave the loop, scale the diagonal and try again
           diagScale = 2.0 * counterFailed; // killme - do smarter -> Saad
           counterFailed++;
           spd = false;
           (*cla) << "\n Perform scaling of system matrix: " << counterFailed  << endl;
           break; // cancel the k loop
        }

        // now l_kk is correct
        l_kk = sqrt(l_kk);
        // store it in our sparse structure
        colCM[k][k] = l_kk;

        // compute all the l_ki values right from the diagonal element
        for(UInt i = k+1; i <= size_; i++)
        {
            // sum the l_ji * l_jk pairs
            T l_ki;
            sysmat.GetMatrixEntry(k, i, &l_ki);

            // we do this by fetching all i_columns and compare the elements with the j_colum
            Col_Map& i_col = colCM[i];
            for (Col_Map::iterator iter = i_col.begin(); iter != i_col.end(); iter++)
            {
                UInt j_row = (*iter).first;
                searchCM = k_col.find(j_row);
                if(searchCM != k_col.end())
                {
                    cout << "l[j=" << j_row << "][i=" << i << "]=" << (*iter).second << "\t l[j=" << j_row << "][k=" << k << "]" << (*seachCM).second << endl;
                    l_ki -=  (*iter).second * (*seachCM).second;
                }
                else
                {
                    cout << "l[j=" << j_row << "][i=" << i << "]=" << (*iter).second << "\t l[j=" << j_row << "][k=" << k << "] = <empty>" << endl;
                }
            }

            l_ki /= l_kk;
            // store in sparse structure
            colCM[i][k] = l_ik; // note that we store at "row" k in the map of column k
        } // end of l_ki loop
      } // end of k loop (equation look)
      (*cla) << "\n \n " << std::endl;
    } // end of while loop (to get postiv definite precond. matrix)


    UInt nnzA = (size_ + sysmat.GetNnz() ) / 2;

    // Initialise internal data arrays
    NEWARRAY( rptrU_, UInt, size_+1 );
    NEWARRAY( cidxU_, UInt, nnzA );
    NEWARRAY( dataU_, T, nnzA );


    // =================
    //  Get matrix data
    // =================
    const Integer *cidxA = sysmat.GetColPointer();
    const Integer *rptrA = sysmat.GetRowPointer();
    const T *dataA = sysmat.GetDataPointer();

    //copy the structure of A, we will fill this structure with the decomposition
    for (UInt k=1; k<=size_+1; k++) {
      rptrU_[k] = rptrA[k];
    }

    for (UInt k=1; k<=nnzA; k++) {
      cidxU_[k] = cidxA[k];
    }



  }



  // copy the upper triangonal L factorization to a SCRC structure
    UInt nnzA = (size_ + sysmat.GetNnz() ) / 2;

    // Initialise internal data arrays
    NEWARRAY( rptrU_, UInt, size_+1 );
    NEWARRAY( cidxU_, UInt, nnzA );
    NEWARRAY( dataU_, T, nnzA );


    // =================
    //  Get matrix data
    // =================
    const Integer *cidxA = sysmat.GetColPointer();
    const Integer *rptrA = sysmat.GetRowPointer();
    const T *dataA = sysmat.GetDataPointer();

    //copy the structure of A, we will fill this structure with the decomposition
    for (UInt k=1; k<=size_+1; k++) {
      rptrU_[k] = rptrA[k];
    }

    for (UInt k=1; k<=nnzA; k++) {
      cidxU_[k] = cidxA[k];
    }


*/
  // *******************************
  //   Apply: apply preconditioner
  // *******************************
  template <typename T>
  void IC0Precond<T>::Apply( const SCRS_Matrix<T> &sysmat,
                             const Vector<T> &r, Vector<T> &z ) const {

    //set values of solution to RHS
    for (UInt k=0; k<size_; k++) {
      z[k] = r[k];
    }

    // do forward substitution
    z[0] /= dataU_[0];

    UInt idx;
    T  zprev;

    for ( UInt k=1; k<size_; k++) {
      zprev = z[k-1];

      for ( UInt i=rptrU_[k-1]+1; i<rptrU_[k]; i++ ) {
        idx = cidxU_[i];
        z[idx] -= zprev * dataU_[i];
      }

      z[k] /= dataU_[rptrU_[k]];

    }

    //solve backward
    z[size_-1] /= dataU_[rptrU_[size_-1]];

    UInt k = size_-1;
    while ( k > 0 ) {
      k--;
      for ( UInt i=rptrU_[k]+1; i<rptrU_[k+1]; i++ ) {
	z[k] -= z[cidxU_[i]] * dataU_[i];
      }

      z[k] /= dataU_[rptrU_[k]];
    }
  }


  template <typename T>
  UInt IC0Precond<T>::CheckPositivity(Double val) {

    if (val > 0)
      return 1;
    else
      return 0;
  }

  template <typename T>
  UInt IC0Precond<T>::CheckPositivity(Complex val) {

    if (std::norm(val) > 0)
      return 1;
    else
      return 0;
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class IC0Precond<Double>;
  template class IC0Precond<Complex>;
#endif

} // namespace
