// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEENTRYMANIPULATOR_HH
#define OLAS_BASEENTRYMANIPULATOR_HH


#include "utils/utils.hh"
#include "matvec/vector.hh"



namespace OLAS {


  // Forward declaration of classes
  class StandardSystem;
  class StdMatrix;
  class SparseVector;
  class BaseIDBC_Handler;


  //! Auxilliary class for manipulating matrix and vector entries

  //! This is the base class for a collection of child classes that act as
  //! auxilliaries for the children of the BaseSystem class. The child classes
  //! fix a problem of the CFS++/%OLAS interface. CFS++ does not know about
  //! tiny matrices or tiny vectors and scarcely supports use of the Complex
  //! data type. Thus floating point values are passed between CFS++ and OLAS
  //! mostly in the form of plain Double arrays. The children of this class
  //! allow conversion of the double arrays into Complex numbers or tiny
  //! matrices/vectors and relieve the children of BaseSystem from this
  //! task. The classes also provide some specialised methods used in the
  //! generation of the linear system.
  struct BaseEntryManipulator {

  public:

    // ======================================================================
    // CONSTRUCTION and DESTRUCTION
    // ======================================================================

    //@{ \name Methods for construction and destruction

    //! Constructor
    BaseEntryManipulator() {
      ENTER_FCN( "BaseEntryManipulator::BaseEntryManipulator" );
      solBufferIsValid_ = false;
      rhsBufferIsValid_ = false;
    }
 
    //! Destructor
    virtual ~BaseEntryManipulator() {
      ENTER_FCN( "BaseEntryManipulator::~BaseEntryManipulator" );
    }

    //@}

    // ======================================================================
    // ASSEMBLY of LINEAR SYSTEM
    // ======================================================================

    //@{ \name Methods used in the assembly of the linear system

    //! Insert local element matrix into global finite element matrix
    virtual void SetElementMatrix( FEMatrixType matrix_id,
                                   const PdeIdType pdeID1,
                                   const PdeIdType pdeID2,
                                   StdMatrix *stdMat,
                                   StdMatrix *counterPart,
                                   BaseIDBC_Handler *idbcHandler,
                                   Double *elemMat,
                                   Integer *connect1, UInt length1,
                                   Integer *connect2, UInt length2,
                                   UInt limit1, UInt limit2,
                                   bool setCounterPart ) = 0;

    //! Only insert counter part of local element matrix into global one
    virtual void SetCounterPartOnly( FEMatrixType matrix_id,
                                     const PdeIdType pdeID1,
                                     const PdeIdType pdeID2,
                                     StdMatrix *stdMat,
                                     StdMatrix *counterPart,
                                     BaseIDBC_Handler *idbcHandler,
                                     Double *elemMat,
                                     Integer *connect1, UInt length1,
                                     Integer *connect2, UInt length2,
                                     UInt limit1, UInt limit2 ) = 0;

    //! Method for insertion of local element vector into right-hand side

    //! This auxilliary method performs the insertion of an element-local
    //! right-hand side vector into the global right-hand side vector of
    //! the linear system. In doing so only values for free degrees of
    //! freedom (equation numbers in [1,limit]) are inserted. Values of
    //! degrees of freedom fixed by (in)homogeneous Dirichlet boundary
    //! conditions are discarded, while values of slaves (marked by a
    //! negative equation number) are added to the corresponding master
    //! in the constraint case.
    //! \param rhs      pointer to right-hand side vector of linear system
    //! \param elemRHS  array containing right-hand side vector computed
    //!                 for the finite element; in the complex case we
    //!                 assume that the array contains all the real parts
    //!                 followed in positions [elemSize+1,2*elemSize] by
    //!                 all the imaginary parts
    //! \param connect  array containing the equation numbers assigned
    //!                 to the degrees of freedom of the finite element
    //! \param elemSize length of the connect array
    //! \param limit    the interval [1,limit] contains equation numbers
    //!                 of free degrees of freedom
    //! \note For using the penalty approach of handling inhomogeneous
    //!       Dirichlet boundary conditions simply set limit = size of
    //!       the linear system
    virtual void SetElementRHS( SparseVector *rhs, Double *elemRHS,
                                Integer *connect, UInt elemSize,
                                UInt limit ) = 0;
    
    //! Get an entry of a matrix (real)
    virtual void GetMatrixEntry( FEMatrixType matrix_id,
                                 const PdeIdType pdeID1,
                                 const PdeIdType pdeID2,
                                 StdMatrix *stdMat,
                                 BaseIDBC_Handler *idbcHandler,
                                 Double & entry,
                                 Integer eqnNr1, 
                                 Integer eqnNr2, 
                                 UInt limit1, UInt limit2 ) = 0;
    
    //! Get an entry of a matrix (complex)
    virtual void GetMatrixEntry( FEMatrixType matrix_id,
                                 const PdeIdType pdeID1,
                                 const PdeIdType pdeID2,
                                 StdMatrix *stdMat,
                                 BaseIDBC_Handler *idbcHandler,
                                 Complex & entry,
                                 Integer eqnNr1, 
                                 Integer eqnNr2, 
                                 UInt limit1, UInt limit2 ) = 0;

    //! Set directly an entry of a  matrix (real)
    virtual void SetMatrixEntry( FEMatrixType matrix_id,
                                 const PdeIdType pdeID1,
                                 const PdeIdType pdeID2,
                                 StdMatrix *stdMat,
                                 BaseIDBC_Handler *idbcHandler,
                                 Double entry,
                                 Integer eqnNr1, 
                                 Integer eqnNr2, 
                                 UInt limit1, UInt limit2,
                                 bool setCounterPart ) = 0;
    
    //! Set directly an entry of a matrix (complex)
    virtual void SetMatrixEntry( FEMatrixType matrix_id,
                                 const PdeIdType pdeID1,
                                 const PdeIdType pdeID2,
                                 StdMatrix *stdMat,
                                 BaseIDBC_Handler *idbcHandler,
                                 Complex entry,
                                 Integer eqnNr1, 
                                 Integer eqnNr2,
                                 UInt limit1, UInt limit2,
                                 bool setCounterPart ) = 0;
    //@}
    

    // ======================================================================
    // GETTER METHODS
    // ======================================================================

    //! \name Methods related to passing vectors back to CFS++
    //@{

    //@{
    //! Return the pointer to the current solution of a PDE

    //! This method passes the current solution of one given PDE
    //! as a pointer to a buffer containing the solution entries. 
    //! If no identifier is given, the pointer to the beginning of the
    //! complete solution array is passed.
    //!
    //! \param sol   pointer to solution vector of linear system
    //! \param ptSol pointer (0-based) to the solution buffer
    //! \param identifierPDE identifier for PDE related to sub-graph
    //! \return number of array entries
    //! 
    //! \note The return buffer is guaranteed to retain the current solution
    //! until the next call of this method (after solving the next step)!
    virtual void GetSolutionVal( SparseVector *sol, Double* &ptSol, 
                                 const PdeIdType identifierPDE
                                 = NO_PDE_ID ) = 0;

    virtual void GetSolutionVal( SparseVector *sol, Complex* &ptSol, 
                                 const PdeIdType identifierPDE 
                                 = NO_PDE_ID ) = 0;
    //@}

    //! Invalidates the buffer containing the solution values
    
    //! This method invalidates the buffer containing the solution array.
    //! This causes an update of the buffer the next time 
    //! GetSolutionVal()-method is called.
    void InvalidateSolBuffer() {
     solBufferIsValid_ = false;
    }

    //@{
    //! Return the pointer to the current rhs value of a PDE
    
    //! This method passes the current rhs as a pointer to buffer with
    //! right hand side entries.
    //! \param rhs   pointer to right-hand side vector of linear system
    //! \param ptRhs pointer (0-based) to the buffer with rhs values
    //! \param identifierPDE identifier for PDE related to sub-graph
    //! \return number of array entries
    //! \note The return buffer is guaranteed to retain the current rhs
    //! until the next call of this method!
    virtual void GetRHSVal( SparseVector *rhs, Double* &ptRhs, 
                            const PdeIdType identifierPDE = NO_PDE_ID) = 0;

    virtual void GetRHSVal( SparseVector *rhs, Complex* &ptRhs, 
                            const PdeIdType identifierPDE = NO_PDE_ID) = 0;
    //@}

    //! Invalidates the buffer containing the rhs values

    //! This method invalidates the buffer containing the right hand side 
    //! array. This causes an update of the buffer the next time 
    //! GetRhsVal()-method is called.
    void InvalidateRhsBuffer() {
      rhsBufferIsValid_ = false;
    }
    //@}


    // ======================================================================
    // SPECIAL RHS METHODS
    // ======================================================================

    //! \name Special methods for manipulation of right-hand side vector
    //@{

    //@{
    //! Set single entry of right-hand side vector
    virtual void SetNodeRHS( SparseVector *rhs, Double val, 
                             Integer node ) = 0;

    virtual void SetNodeRHS( SparseVector *rhs, Complex val, 
                             Integer node ) = 0;
    //@}

    //! Initialise right-hand side with given vector
    virtual void InitRHS( SparseVector *rhs, const Double *newRHS) = 0;
                          
    //! ...
    virtual void UpdateRHS( SparseVector *rhs, StdMatrix *stdMat,
                            Double *fup ) = 0;

    //@}


    // ======================================================================
    // GENERAL VECTOR METHODS
    // ======================================================================

    //! \name General methods for setting matrix/vector entries
    //@{

    //! Set entry in a vector to specified value
    virtual void SetVectorEntry( SparseVector *vec, UInt index,
                                 Double &newVal ) = 0;

    //! Set entry in a vector to specified value
    virtual void SetVectorEntry( SparseVector *vec, UInt index,
                                 Complex &newVal ) = 0;

    //! Add a value to a diagonal matrix entry

    //! This method allows to directly add a value to the value of an entry
    //! on the main diagonal of a certain of the different problem matrices.
    //! \param stdMat   matrix to be manipulated
    //! \param eqnNum   determines the diagonal entry to be manipulated, i.e.
    //!                 the parameter gives the equation number of the unknown
    //!                 coupled to the row/column of the entry; e.g. for
    //!                 \f$\mbox{eqnNum}=k\f$ it is \f$a_{kk}\f$ that will be
    //!                 altered
    //! \param val      zero-based array containing the numerical value that
    //!                 is to be added to the matrix entry; in the case of
    //!                 real entries, the first array value is used, in the
    //!                 case of complex entries the first entry is used as
    //!                 real and the second as imaginary part
    virtual void AddToDiagMatrixEntry( StdMatrix *stdMat, Integer eqnNum,
                                       Double *val ) = 0;
    //@}


    // ======================================================================
    // METHODS FOR HANDLING IDBCs
    // ======================================================================

    //! \name Methods for handling IDBCs
    //! The following methods are used by the IDBC_HandlerPenalty class to
    //! incorporate inhomogeneous Dirichlet boundary conditions into the
    //! linear system with the help of the penalty approach.

    //@{
    //! Modify system matrix for penalty approach
    virtual void AdaptSystemMatrix( StdMatrix &stdMat, UInt *dirichletEQN,
                                    UInt numIDBC,
                                    Double &penaltyTerm ) = 0;

    //! Modify right-hand side vector following penalty approach
    virtual void AdaptRHSForIDBC( SparseVector &rhs,
                                  UInt *dirichletEQN,
                                  SparseVector &dirichletValue,
                                  Double &penaltyTerm,
                                  UInt numIDBC ) = 0;
    //@}

  protected:
  
    //! flag for state of solution buffer
    bool solBufferIsValid_;

    //! flag for state of rhs buffer
    bool rhsBufferIsValid_;
    
  };

}

#endif
