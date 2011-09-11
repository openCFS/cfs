// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/sbmmatrix.hh"
#include "MatVec/sbmvector.hh"

#include "OLAS/algsys/algebraicSys.hh"
#include "DataInOut/Logging/cfslog.hh"

#include "idbchandlerpenalty.hh"


DECLARE_LOG(idbcPenalty)
DEFINE_LOG(idbcPenalty, "idbcPenalty")

namespace CoupledField {

  template <typename T>
  IDBC_HandlerPenalty<T>::IDBC_HandlerPenalty( StdVector<UInt> numIDBC ) {

    // Log some information
    (*cla) << "\n IDBC_HandlerPenalty: Administrating "
           << numIDBC << " inhom. Dirichlet BCs.\n\n";

    // Initialise penalty term
    penaltyTerm_ = 0.0;

    // Determine entry type from template parameter
    eType_ = BaseMatrix::NOENTRYTYPE;
    if ( AssocType<T>::tagM == AssocType<Double>::tagM ) {
      eType_ = BaseMatrix::DOUBLE;
    }
    else if ( AssocType<T>::tagM == AssocType<Complex>::tagM ) {
      eType_ = BaseMatrix::COMPLEX;
    }
    else {
      EXCEPTION( "Internal template error! No swearing please!" );
    }

    // Generate vector for storing Dirchlet values
    SingleVector *stdVec = NULL;
    BaseVector *bVec = NULL;
    SBM_Vector *sbmVec = new SBM_Vector( numIDBC.GetSize() );
    for ( UInt i = 0; i < numIDBC.GetSize(); i++ ) {

      // If there are IDBCs generate and insert sub-vector
      //
      // NOTE: We use SPARSE_NONSYM as matrix storage type in order to
      //       obtain a Vector<T> object
      if ( numIDBC[i] > 0 ) {
        bVec = GenerateSingleVectorObject( BaseMatrix::SPARSE_NONSYM, eType_, numIDBC[i]);
        stdVec = dynamic_cast<SingleVector*>(bVec);
        sbmVec->SetSubVector( stdVec, i );
      }
    }
    dirichletValue_ = sbmVec;

  }


  // **************
  //   Destructor
  // **************
  template <typename T>
  IDBC_HandlerPenalty<T>::~IDBC_HandlerPenalty() {

    // Delete vector of Dirichlet values
    delete dirichletValue_;
  }


  // ***********************
  //   InitDirichletValues
  // ***********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::InitDirichletValues() {
    dirichletValue_->Init();
  }


  // *********************
  //   AdaptSystemMatrix
  // *********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::AdaptSystemMatrix( SBM_Matrix &sysMat ) {

    // Compute penalty term
    penaltyTerm_ = sysMat.GetMaxDiag() * 1e12;
      
    StdMatrix *stdMat = NULL;
    const unsigned int srows(sysMat.GetNumRows());
    for ( UInt i = 0; i < srows; i++ ) {
      stdMat = sysMat.GetPointer(i,i);
      if ( dirichletEqns_[i].GetSize() > 0 && stdMat != NULL ) {
        StdVector<UInt> & myEqns = dirichletEqns_[i];
        UInt numIDBC = myEqns.GetSize();
        for( UInt j = 0; j < numIDBC; ++j ) {
          stdMat->SetDiagEntry( myEqns[j]-1, penaltyTerm_ );  
        }
      }
    }

    // Be verbose
    (*cla) << "\n IDBC_HandlerPenalty:\n"
           << " -> Adapted system matrix\n"
           << " -> Penalty term = "
           << std::scientific << penaltyTerm_ << '\n'
           << std::endl;
    cla->unsetf( std::ios::scientific );
  }


  // ****************
  //   AddIDBCToRHS
  // ****************
  template <typename T>
  void IDBC_HandlerPenalty<T>::AddIDBCToRHS( SBM_Vector *rhs ) {
    

    // We need pointers to sub-vectors
    SingleVector *stdVec = NULL;
    SingleVector *stdVal = NULL;

    // Now loop over all blocks / sub-vectors and
    // set the Dirichlet values
    T entry;
    for ( UInt i = 0; i < rhs->GetSize() ; i++ ) {

      // obtain sub-vectors
      stdVec = rhs->GetPointer( i );
      stdVal = dirichletValue_->GetPointer( i );

      // Now adapt the sub-vector
      if ( dirichletEqns_[i].GetSize() > 0 
          && stdVec != NULL && stdVal != NULL ) {
        StdVector<UInt> & myEqns = dirichletEqns_[i];
        UInt numIDBC = myEqns.GetSize();

        for( UInt j = 0; j < numIDBC; j++ ) {
          stdVal->GetEntry(j, entry);
          entry *=  penaltyTerm_;
          stdVec->SetEntry( myEqns[j] - 1, entry );
        }
      }
    }
  }


  // ***********
  //   SetIDBC
  // ***********
  template <typename T>
  void IDBC_HandlerPenalty<T>::SetIDBC( UInt rowBlock, UInt rowNum,
                                        const T &val ) {
    
    // get hold of block specific data
    std::map<Integer, UInt> & myIndices = bcIndices_[rowBlock];
    StdVector<UInt> & myEqns = dirichletEqns_[rowBlock];
    SingleVector *stdVec = dirichletValue_->GetPointer( rowBlock );
    
    // Check, if equation was already set and determine index
    UInt index = 0;
    std::map<Integer, UInt>::iterator it = myIndices.find(rowNum); 
    if( it != myIndices.end() ) {
      index = it->second; 
    } else {
      myIndices[rowNum] = myEqns.GetSize();
      index = myEqns.GetSize();
      myEqns.Push_back(rowNum);
    }
    stdVec->SetEntry(index,val);
  }
  
  // ***********
   //   GetIDBC
   // ***********
   template <typename T>
   void IDBC_HandlerPenalty<T>::GetIDBC( UInt rowBlock, UInt rowNum,
                                         T &val ) {
     
     // get hold of block specific data
     std::map<Integer, UInt> & myIndices = bcIndices_[rowBlock];
     SingleVector *stdVec = dirichletValue_->GetPointer( rowBlock );
     
     // Check, if equation was already set and determine index
     UInt index = 0;
     std::map<Integer, UInt>::iterator it = myIndices.find(rowNum); 
     if( it != myIndices.end() ) {
       index = it->second; 
     } else {
       EXCEPTION("Index " << rowNum << " in sbmBlock #"
                <<  rowBlock << " is no Dirichlet entry"); 
     }
     stdVec->GetEntry(index,val);
   }


  // *****************
  //   SetDofsToIDBC
  // *****************
  template <typename T>
  void IDBC_HandlerPenalty<T>::SetDofsToIDBC( SBM_Vector *vec ) {

    // Loop over all sub-vectors
    SingleVector *stdVec = NULL;
    SingleVector *stdVal = NULL;
    T aux;

    for ( UInt i = 0; i < vec->GetSize(); i++ ) {
      stdVec = vec->GetPointer( i );
      stdVal = dirichletValue_->GetPointer( i );

      // leave, if no Dirichlet values are defined for this block
      if( stdVal == NULL ) 
        continue;

      StdVector<UInt> & myEqns = dirichletEqns_[i];

      // Insert values
      for ( UInt j = 0; j < stdVal->GetSize(); j++ ) {
        stdVal->GetEntry( j, aux );
        stdVec->SetEntry( myEqns[j]-1, aux );
      }
    }
  }

  // **********************
  //   PrintInternalState
  // **********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::PrintInternalState( std::ostream &os ) {

    REFACTOR;

//      os << " IDBC_HandlerPenalty --- Internal State:\n"
//         << " --------------------------------------\n"
//         << " numIDBC_ = " << numIDBC_ << '\n' << std::endl;
//
//      SingleVector *stdVal = NULL;
//      SBM_Vector* sbmVec = dynamic_cast<SBM_Vector*>(dirichletValue_);
//      if(sbmVec == NULL) EXCEPTION("Invalid cast attempt!");
//      UInt aux = 0;
//      UInt idx = 0;
//      T val;
//
//      // Compute column widths for pretty-printing
//      aux = numIDBC_ + 1;
//      UInt cw1 = aux > 0 ? (UInt)std::log10( (float)aux ) + 1 : 1;
//      UInt cw2 = 5;
//      // UInt cw3 = 4; // TODO: Check if this is still needed
//
//      for ( UInt i = 0; i < dirichletValue_->GetSize(); i++ ) {
//
//        // Compute number of IDBCs
//        aux = offset_[i+1] - offset_[i];
//
//        // Get hold of sub-vector
//        stdVal = sbmVec->GetPointer(i);
//
//        // Print IDBCs for this PDE
//        os << " PDE # " << i << " / numIDBCs = " << aux << std::endl;
//        if ( aux > 0 ) {
//          os << " # | eqnNo | comp | value" << std::endl;
//          for ( UInt j = 0; j < aux; j++ ) {
//
//            stdVal->GetEntry( j, val );
//            idx = offset_[i] + j;
//
//            os << std::setfill( ' ' )
//               << std::setw(cw1) << idx << " | "
//               << std::setw(cw2) << dirichletEQN_[idx] << " | "
//               << val << std::endl;
//          }
//        }
//        os << std::endl;
//      }
//
//    os << " --------------------------------------\n" << std::endl;

  }

  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class IDBC_HandlerPenalty<Double>;
    template class IDBC_HandlerPenalty<Complex>;
  #endif  
}
