// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_ENTRYMANIPULATORCOMPLEX_HH
#define OLAS_ENTRYMANIPULATORCOMPLEX_HH


#include <vector>
#include "utils/utils.hh"
#include "matvec/vector.hh"
#include "algsys/baseentrymanipulator.hh"



namespace OLAS {


  // Forward declaration of classes
  class StandardSystem;
  class StdMatrix;
  class StdVector;
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
  struct EntryManipulatorComplex : public BaseEntryManipulator {

  public:

    // ======================================================================
    // CONSTRUCTION and DESTRUCTION
    // ======================================================================

    //@{ \name Methods for construction and destruction

    //! Constructor
    EntryManipulatorComplex() : BaseEntryManipulator() {
      ENTER_FCN( "EntryManipulatorComplex::EntryManipulatorComplex" );
    }
 
    //! Destructor
    ~EntryManipulatorComplex() {
      ENTER_FCN( "EntryManipulatorComplex::~EntryManipulatorComplex" );
    }

    //@}

    // ======================================================================
    // ASSEMBLY of LINEAR SYSTEM
    // ======================================================================

    //@{ \name Methods used in the assembly of the linear system

    //!
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

    //!
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

    //!
    void SetElementRHS( StdVector *rhs, Double *elemRHS,
                        Integer *connect, UInt elemSize, UInt limit );


    //! Get an entry of a matrix (real)
    void GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                         const PdeIdType pdeID2, StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         Double & entry, Integer eqnNr1, Integer eqnDof1,
                         Integer eqnNr2, Integer eqnDof2,
                         UInt limit1, UInt limit2 );
    
    //! Get an entry of a matrix (complex)
    void GetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                         const PdeIdType pdeID2, StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         Complex & entry, Integer eqnNr1, Integer eqnDof1,
                         Integer eqnNr2, Integer eqnDof2,
                         UInt limit1, UInt limit2 );

    //! Set directly an entry of a  matrix (real)
    void SetMatrixEntry( FEMatrixType matrix_id, const PdeIdType pdeID1,
                         const PdeIdType pdeID2, StdMatrix *stdMat,
                         BaseIDBC_Handler *idbcHandler,
                         Double entry,
                         Integer eqnNr1, Integer eqnDof1,
                         Integer eqnNr2, Integer eqnDof2,
                         UInt limit1, UInt limit2, bool setCounterPart );
    
    //! Set directly an entry of a matrix (complex)
    void SetMatrixEntry( FEMatrixType matrix_id,
                         const PdeIdType pdeID1,const PdeIdType pdeID2,
                         StdMatrix *stdMat, BaseIDBC_Handler *idbcHandler,
                         Complex entry,
                         Integer eqnNr1, Integer eqnDof1,
                         Integer eqnNr2, Integer eqnDof2,
                         UInt limit1, UInt limit2, bool setCounterPart );
    //@}


    // ======================================================================
    // ASSEMBLY of LINEAR SYSTEM
    // ======================================================================

    //!
    void SetNodeRHS( StdVector *rhs, Double val, 
                     Integer node, Integer dof );
    //!
    void SetNodeRHS( StdVector *rhs, Complex val, 
                     Integer node, Integer dof );
                          
    //!
    void InitRHS( StdVector *rhs, const Double *newRHS);
                          
    //!
    void UpdateRHS( StdVector *rhs, StdMatrix *stdMat,
                    Double *fup );

    //@{
    //! Set entry in a vector to specified value
    void SetVectorEntry( StdVector *vec, UInt index, UInt component,
                         Double &newVal );
    void SetVectorEntry( StdVector *vec, UInt index, UInt component,
                         Complex &newVal );
    //@}

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
    void GetSolutionVal( StdVector *sol, Double* &ptSol, 
                         const PdeIdType identifierPDE = NO_PDE_ID );

    void GetSolutionVal( StdVector *sol, Complex* &ptSol, 
                         const PdeIdType identifierPDE = NO_PDE_ID );
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
    //!
    //! \param rhs   pointer to right-hand side vector of linear system
    //! \param ptRhs pointer (0-based) to the buffer with rhs values
    //! \param identifierPDE identifier for PDE related to sub-graph
    //! \return number of array entries
    //!
    //! \note The return buffer is guaranteed to retain the current rhs
    //! until the next call of this method!
    void GetRHSVal( StdVector *rhs, Double* &ptRhs, 
                    const PdeIdType identifierPDE = NO_PDE_ID);

    void GetRHSVal( StdVector *rhs, Complex* &ptRhs, 
                    const PdeIdType identifierPDE = NO_PDE_ID);
    //@}

    //! Invalidates the buffer containng the rhs values

    //! This method invalidates the buffer containing the right hand side 
    //! array. This causes an update of the buffer the next time 
    //! GetRhsVal()-method is called.
    void InvalidateRhsBuffer() {
      rhsBufferIsValid_ = false;
    }

    //! Add a value to a diagonal matrix entry

    //! This method allows to directly add a value to the value of an entry
    //! on the main diagonal of a certain of the different problem matrices.
    //! \param stdMat   matrix to be manipulated
    //! \param eqnNum   determines the diagonal entry to be manipulated, i.e.
    //!                 the parameter gives the equation number of the unknown
    //!                 coupled to the row/column of the entry; e.g. for
    //!                 \f$\mbox{eqnNum}=k\f$ it is \f$a_{kk}\f$ that will be
    //!                 altered
    //! \param dof      gives the degree of freedom in the diagonal entry
    //!                 that is to be changed; in the case of scalar entries
    //!                 this should be 1, only in the case of block entries
    //!                 may this be larger than 1
    //! \param val      zero-based array containing the numerical value that
    //!                 is to be added to the matrix entry; in the case of
    //!                 real entries, the first array value is used, in the
    //!                 case of complex entries the first entry is used as
    //!                 real and the second as imaginary part
    void AddToDiagMatrixEntry( StdMatrix *stdMat, Integer eqnNum,
                               Integer dof, Double *val );


    // ======================================================================
    // METHODS FOR HANDLING IDBCs
    // ======================================================================

    //! \name Methods for handling IDBCs
    //! The following methods are used by the IDBC_HandlerPenalty class to
    //! incorporate inhomogeneous Dirichlet boundary conditions into the
    //! linear system with the help of the penalty approach.

    //@{
    //! Modify system matrix for penalty approach
    void AdaptSystemMatrix( StdMatrix &stdMat, UInt *dirichletEQN,
                            UInt *dirichletComponent,
                            UInt numIDBC,
                            Double &penaltyTerm );

    //! Modify right-hand side vector following penalty approach
    void AdaptRHSForIDBC( StdVector &rhs,
                          UInt *dirichletEQN,
                          UInt *dirichletComponent,
                          StdVector &dirichletValue,
                          Double &penaltyTerm,
                          UInt numIDBC );
    //@}

  private:

    //! return buffer for solution values
    Vector<Complex> solBuffer_;

    //! return buffer for rhs values
    Vector<Complex> rhsBuffer_;

    //@{
    //! Auxilliary vector
    std::vector<UInt> rowIndList_ ;
    std::vector<UInt> colIndList1_;
    std::vector<UInt> colIndList2_;
    //@}

  };

}

#endif
