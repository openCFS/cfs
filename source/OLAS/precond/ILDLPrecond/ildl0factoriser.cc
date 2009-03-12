// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <vector>

#include "General/exception.hh"

#include "MatVec/scrs_matrix.hh"

#include "ildl0factoriser.hh"

namespace CoupledField {


  // ***********************
  //   Default constructor
  // ***********************
  template <class T>
  ILDL0Factoriser<T>::ILDL0Factoriser() {


    EXCEPTION( "Default constructor of ILDL0Factoriser call was called! "
             << "This constructor is forbidden!" );

  }


  // ************************
  //   Standard constructor
  // ************************
  template <class T>
  ILDL0Factoriser<T>::ILDL0Factoriser( OLAS_Params *myParams,
                                       OLAS_Report *myReport ) {


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
  ILDL0Factoriser<T>::~ILDL0Factoriser() {


  }


  // *************
  //   Factorise
  // *************
  template <class T>
  void ILDL0Factoriser<T>::Factorise( SCRS_Matrix<T> &sysMat,
                                      std::vector<T> &dataD,
                                      std::vector<UInt> &rptrU,
                                      std::vector<UInt> &cidxU,
                                      std::vector<T> &dataU,
                                      bool newPattern ) {


    this->sysMatDim_ = sysMat.GetNumCols();
    Integer nnzA = (sysMat.GetNnz() + this->sysMatDim_ ) / 2;
    //std::cout << "Dim=" << this->sysMatDim_ << " NNZ=" << nnzA << std::endl;

    dataD.reserve( this->sysMatDim_ + 1 );
    rptrU.reserve( this->sysMatDim_ + 2 );
    cidxU.reserve( nnzA - this->sysMatDim_ +1 ); // we do not store diagonal elements
    dataU.reserve( nnzA - this->sysMatDim_ +1 ); // since these are 1

    // OLAS uses one-base indexing, so we push_back a zero up front
    dataD.push_back( 0 );
    rptrU.push_back( 0 );
    cidxU.push_back( 0 );
    dataU.push_back( 0 );

    // do the factorization
    FactoriseNumerics( sysMat, dataD, rptrU, cidxU, dataU );

  }

 
  // *********************
  //   FactoriseNumerics
  // *********************
  template<typename T>
  void ILDL0Factoriser<T>::FactoriseNumerics( SCRS_Matrix<T> &sysMat,
                                              std::vector<T> &dataD,
                                              std::vector<UInt> &rptrU,
                                              std::vector<UInt> &cidxU,
                                              std::vector<T> &dataU ) {


    // Shall we be verbose?
    // COMPWARNING: unused variable bool logging = this->myParams_->GetIntValue( "ILDLPRECOND_logging" ) > 0;

    // =================
    //  Get matrix data
    // =================
    const UInt *cidxA = sysMat.GetColPointer();
    const UInt *rptrA = sysMat.GetRowPointer();
    const T *dataA = sysMat.GetDataPointer();

    //write out the system matrix
    for (UInt k=1; k<=this->sysMatDim_; k++) {
      Integer numCol = rptrA[k+1] - rptrA[k]; 
      // Integer startIdx = rptrA[k]; TODO: Check if this is still needed
      for (Integer j=0; j<numCol; j++) {
        //std::cout << dataA[startIdx+j] << "  ";
      }
      //std::cout << std::endl;
    }

    //std::cout << "Do factorise" << std::endl;

    //help vector; same structure as dataU
    Integer nnzA = (sysMat.GetNnz() + this->sysMatDim_ ) / 2;
    std::vector<T> dataH;
    dataH.resize(nnzA - this->sysMatDim_+1);
    dataH.push_back( 0 );

		// COMPWARNING: unused variable Integer idx 
    Integer idxK1, idxK2, idxI1, idxI2;

    //help vectors
    std::vector<Integer> colIdxRowK;
    std::vector<Integer> colIdxRowI;
    colIdxRowK.resize(this->sysMatDim_+1); 
    colIdxRowI.resize(this->sysMatDim_+1); 

    colIdxRowK.push_back( 0 );
    colIdxRowI.push_back( 0 );

    T val;
    //loop over all rows
    for (UInt j=1; j<=this->sysMatDim_; j++) {
      //std::cout << "EQj: " << j << std::endl;
      for (UInt i=1; i<=j-1; i++) {
        //std::cout << "EQi: " << i << std::endl;
        //get entry positions if row i
        idxI1 = rptrA[i] + 1;
        idxI2 = rptrA[i+1] - 1;
        //std::cout << "idxI1=" << idxI1 << " idxI2=" << idxI2 << std::endl;
        for (Integer l=idxI1; l<=idxI2; l++) {
          colIdxRowI[cidxA[l]] = l;
        }

        //for (UInt pos=1; pos<=this->sysMatDim_; pos++) {
        //  std::cout <<  colIdxRowI[pos] << "  ";
        //}
        std::cout << std::endl;

        if (colIdxRowI[j] > 0 ) {
          dataH[idxI1-1] = dataA[idxI1];

          //std::cout << "dataA=" <<  dataA[idxI1] << std::endl;

          val = 0.0;
          for (UInt k=1; k<=i-1; k++) {
            //get entry positions if row i
            idxK1 = rptrA[k] + 1;
            idxK2 = rptrA[k+1] - 1;
            //std::cout << "idxK1=" << idxK1 << " idxK2=" << idxK2 << std::endl;

            for (Integer l=idxK1; l<=idxK2; l++) {
              colIdxRowK[cidxA[l]] = l;
            }

            if (colIdxRowK[i] > 0 && colIdxRowK[j] ) {
              val = val + dataU[idxK1]; // * dataH[
            }

            for (Integer l=idxK1; l<=idxK2; l++) {
              colIdxRowK[cidxA[l]] = 0;
            }
          }

        }

        for (Integer l=idxI1; l<=idxI2; l++) {
          colIdxRowI[cidxA[l]] = 0; 
        }

      } // for i
    }  // for j


    //std::cout << "OK, that's it" << std::endl;
  }

  // Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class ILDL0Factoriser<Double>;
    template class ILDL0Factoriser<Complex>;
  #endif

}
