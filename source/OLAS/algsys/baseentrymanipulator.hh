// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEENTRYMANIPULATOR_HH
#define OLAS_BASEENTRYMANIPULATOR_HH


#include "MatVec/vector.hh"

#include "General/exception.hh"

namespace CoupledField {


  // Forward declaration of classes
  class StandardSystem;
  class StdMatrix;
  class SingleVector;
  class BaseIDBC_Handler;


  //! Auxilliary class for manipulating matrix and vector entries

  //! This is the base class for a collection of child classes that act as
  //! auxilliaries for the children of the AlgebraicSys class. The child classes
  //! fix a problem of the CFS++/%OLAS interface. CFS++ does not know about
  //! tiny matrices or tiny vectors and scarcely supports use of the Complex
  //! data type. Thus floating point values are passed between CFS++ and OLAS
  //! mostly in the form of plain Double arrays. The children of this class
  //! allow conversion of the double arrays into Complex numbers or tiny
  //! matrices/vectors and relieve the children of AlgebraicSys from this
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
      solBufferIsValid_ = false;
      rhsBufferIsValid_ = false;
    }
 
    //! Destructor
    ~BaseEntryManipulator() {
    }

    //@}

    // ======================================================================
    // ASSEMBLY of LINEAR SYSTEM
    // ======================================================================

    //@{ \name Methods used in the assembly of the linear system

    //! Insert local element matrix into global finite element matrix
    template<typename T>
    void SetElementMatrix( FEMatrixType matrix_id,
                           const FeFctIdType pdeID1,
                           const FeFctIdType pdeID2,
                           StdMatrix *stdMat,
                           StdMatrix *counterPart,
                           BaseIDBC_Handler *idbcHandler,
                           const Matrix<T>& elemMat,
                           const StdVector<Integer>& connect1,
                           const StdVector<Integer>& connect2,
                           UInt limit1, UInt limit2,
                           bool setCounterPart);

    //! Only insert counter part of local element matrix into global one
    template<typename T>
    void SetCounterPartOnly( FEMatrixType matrix_id,
                             const FeFctIdType pdeID1,
                             const FeFctIdType pdeID2,
                             StdMatrix *stdMat,
                             StdMatrix *counterPart,
                             BaseIDBC_Handler *idbcHandler,
                             const Matrix<T>& elemMat,
                             const StdVector<Integer>& connect1,
                             const StdVector<Integer>& connect2,
                             UInt limit1, UInt limit2 );

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
    template<typename T>
    void SetElementRHS( SingleVector *rhs, const Vector<T>& elemRhs,
                        const StdVector<Integer>& connect,
                        UInt limit );
    
    //! Get an entry of a matrix
    template<typename T>
    void GetMatrixEntry( FEMatrixType matrix_id,
                         const FeFctIdType pdeID1,
                         const FeFctIdType pdeID2,
                         StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         T & entry,
                         Integer eqnNr1, 
                         Integer eqnNr2, 
                         UInt limit1, UInt limit2 );
                         
    //! Set directly an entry of a  matrix
    template<typename T>
    void SetMatrixEntry( FEMatrixType matrix_id,
                         const FeFctIdType pdeID1,
                         const FeFctIdType pdeID2,
                         StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         const T& entry,
                         Integer eqnNr1, 
                         Integer eqnNr2, 
                         UInt limit1, UInt limit2,
                         bool setCounterPart );
    
    //@}
    

    // ======================================================================
    // SPECIAL RHS METHODS
    // ======================================================================

    //! \name Special methods for manipulation of right-hand side vector
    //@{

    //! Set single entry of right-hand side vector
    template<typename T>
    void SetNodeRHS( SingleVector *rhs, const T& val, Integer node );
   
    //! Initialize RHS side with given entries
    template<typename T>
    void InitRHS( SingleVector *rhs, const Vector<T>& newRHS );
    
    //! Initialise sol vector  with given vector
    template<typename T>
    void InitSol( SingleVector *sol, const Vector<T>& newSol);

    //! ...
    template<typename T>
    void UpdateRHS( SingleVector *rhs, StdMatrix *stdMat,
                    const Vector<T>& fup );

    //@}


    // ======================================================================
    // GENERAL VECTOR METHODS
    // ======================================================================

    //! \name General methods for setting matrix/vector entries
    //@{

    //! Set entry in a vector to specified value
    template<typename T>
    void SetVectorEntry( SingleVector *vec, UInt index,
                         const T& newVal );

    
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
    template<typename T>
    void AddToDiagMatrixEntry( StdMatrix *stdMat, Integer eqnNum,
                               const T& val );
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
    void AdaptSystemMatrix( StdMatrix &stdMat, 
                            StdVector<UInt>& dirichletEQN,
                            Double& penaltyTerm );
                            
    //! Modify right-hand side vector following penalty approach
    void AdaptRHSForIDBC( SingleVector &rhs,
                          StdVector<UInt>& dirichletEQN,
                          SingleVector& dirichletValues,
                          Double& penaltyTerm,
                          UInt numIDBC );
    //@}

  protected:
  
    //! flag for state of solution buffer
    bool solBufferIsValid_;

    //! flag for state of rhs buffer
    bool rhsBufferIsValid_;
    
    //@{
    //! Auxilliary vector
    std::vector<UInt> rowIndList_ ;
    std::vector<UInt> colIndList1_;
    std::vector<UInt> colIndList2_;
    //@}

  };

}

#endif
