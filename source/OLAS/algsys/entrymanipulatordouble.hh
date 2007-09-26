// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_ENTRYMANIPULATORDOUBLE_HH
#define OLAS_ENTRYMANIPULATORDOUBLE_HH


#include <vector>
#include "utils/utils.hh"
#include "matvec/vector.hh"
#include "algsys/baseentrymanipulator.hh"



namespace OLAS {


  // Forward declaration of classes
  class StandardSystem;
  class StdMatrix;
  class SparseVector;
  class BaseIDBC_Handler;


  //! Auxilliary class for manipulating matrix and vector entries

  //! This is the base class for a collection of child classes that act as
  //! auxilliaries for the children of the BaseSystem class. The child classes
  //! fix a problem of the CFS++/OLAS interface. CFS++ does not know about
  //! tiny matrices or tiny vectors and scarcely supports use of the Complex
  //! data type. Thus floating point values are passed between CFS++ and OLAS
  //! mostly in the form of plain Double arrays. The children of this class
  //! allow conversion of the double arrays into Complex numbers or tiny
  //! matrices/vectors and relieve the children of BaseSystem from this
  //! task. The classes also provide some specialised methods used in the
  //! generation of the linear system.
  struct EntryManipulatorDouble : public BaseEntryManipulator {

  public:

    // ======================================================================
    // CONSTRUCTION and DESTRUCTION
    // ======================================================================

    //@{ \name Methods for construction and destruction

    //! Constructor
    EntryManipulatorDouble() : BaseEntryManipulator() {
    }
 
    //! Destructor
    ~EntryManipulatorDouble() {
    }

    //@}

    // ======================================================================
    // REMAINING DOC IS INHERITED FROM BASE CLASS!
    // ======================================================================

    void SetElementMatrix( FEMatrixType matrix_id,
                           const PdeIdType pdeID1,
                           const PdeIdType pdeID2,
                           StdMatrix *stdMat,
                           StdMatrix *counterPart,
                           BaseIDBC_Handler *idbcHandler,
                           Double *elemMat,
                           Integer *connect1, UInt length1,
                           Integer *connect2, UInt length2,
                           UInt limit1, UInt limit2,
                           bool setCounterPart );

    void SetCounterPartOnly( FEMatrixType matrix_id,
                             const PdeIdType pdeID1,
                             const PdeIdType pdeID2,
                             StdMatrix *stdMat,
                             StdMatrix *counterPart,
                             BaseIDBC_Handler *idbcHandler,
                             Double *elemMat,
                             Integer *connect1, UInt length1,
                             Integer *connect2, UInt length2,
                             UInt limit1, UInt limit2 );

    void SetElementRHS( SparseVector *rhs, Double *elemRHS,
                        Integer *connect, UInt elemSize, UInt limit );

    //! Get an entry of a matrix (real)
    void GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                         const PdeIdType pdeID2, StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         Double & entry, Integer eqnNr1, 
                         Integer eqnNr2, 
                         UInt limit1, UInt limit2 );
    
    //! Get an entry of a matrix (complex)
    void GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                         const PdeIdType pdeID2, StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         Complex & entry, Integer eqnNr1, 
                         Integer eqnNr2, 
                         UInt limit1, UInt limit2 );

    //! Set directly an entry of a  matrix (real)
    void SetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                         const PdeIdType pdeID2, StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         Double entry,
                         Integer eqnNr1, 
                         Integer eqnNr2, 
                         UInt limit1, UInt limit2, bool setCounterPart );
    
    //! Set directly an entry of a matrix (complex)
    void SetMatrixEntry( FEMatrixType matrix_id,
                         const PdeIdType pdeID1,const PdeIdType pdeID2,
                         StdMatrix *stdMat, BaseIDBC_Handler *idbcHandler,
                         Complex entry,
                         Integer eqnNr1, 
                         Integer eqnNr2, 
                         UInt limit1, UInt limit2, bool setCounterPart );

    void SetNodeRHS( SparseVector *rhs, Double val, 
                     Integer node );

    void SetNodeRHS( SparseVector *rhs, Complex val, 
                     Integer node);
                          
    void InitRHS( SparseVector *rhs, const Double *newRHS);
                          
    void UpdateRHS( SparseVector *rhs, StdMatrix *stdMat,
                    Double *fup );

    void SetVectorEntry( SparseVector *vec, UInt index,
                         Double &newVal );

    void SetVectorEntry( SparseVector *vec, UInt index,
                         Complex &newVal );

    void GetSolutionVal( SparseVector *sol, Double* &ptSol, 
                         const PdeIdType identifierPDE = NO_PDE_ID );

    void GetSolutionVal( SparseVector *sol, Complex* &ptSol,
                         const PdeIdType identifierPDE = NO_PDE_ID );

    void InvalidateSolBuffer() {
      solBufferIsValid_ = false;
    }

    void GetRHSVal( SparseVector *rhs, Double* &ptRhs, 
                    const PdeIdType identifierPDE = NO_PDE_ID);

    void GetRHSVal( SparseVector *rhs, Complex* &ptRhs, 
                    const PdeIdType identifierPDE = NO_PDE_ID);

    void InvalidateRhsBuffer() {
      rhsBufferIsValid_ = false;
    }

    void AddToDiagMatrixEntry( StdMatrix *stdMat, Integer eqnNum,
                               Double *val );

    void AdaptSystemMatrix( StdMatrix &stdMat, UInt *dirichletEQN,
                            UInt numIDBC,
                            Double &penaltyTerm );

    void AdaptRHSForIDBC( SparseVector &rhs,
                          UInt *dirichletEQN,
                          SparseVector &dirichletValue,
                          Double &penaltyTerm,
                          UInt numIDBC );

  private:

    //! return buffer for solution values
    Vector<Double> solBuffer_;

    //! return buffer for rhs values
    Vector<Double> rhsBuffer_;

    //@{
    //! Auxilliary vector
    std::vector<UInt> rowIndList_ ;
    std::vector<UInt> colIndList1_;
    std::vector<UInt> colIndList2_;
    //@}

  };

}

#endif
