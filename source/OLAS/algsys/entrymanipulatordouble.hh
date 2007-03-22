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
  struct EntryManipulatorDouble : public BaseEntryManipulator {

  public:

    // ======================================================================
    // CONSTRUCTION and DESTRUCTION
    // ======================================================================

    //@{ \name Methods for construction and destruction

    //! Constructor
    EntryManipulatorDouble() : BaseEntryManipulator() {
      ENTER_FCN("EntryManipulatorDouble::EntryManipulatorDouble");
    }
 
    //! Destructor
    ~EntryManipulatorDouble() {
      ENTER_FCN( "EntryManipulatorDouble::~EntryManipulatorDouble" );
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

    void SetNodeRHS( StdVector *rhs, Double val, 
                     Integer node, Integer dof );

    void SetNodeRHS( StdVector *rhs, Complex val, 
                     Integer node, Integer dof );
                          
    void InitRHS( StdVector *rhs, const Double *newRHS);
                          
    void UpdateRHS( StdVector *rhs, StdMatrix *stdMat,
                    Double *fup );

    void SetVectorEntry( StdVector *vec, UInt index, UInt component,
                         Double &newVal );

    void SetVectorEntry( StdVector *vec, UInt index, UInt component,
                         Complex &newVal );

    void GetSolutionVal( StdVector *sol, Double* &ptSol, 
                         const PdeIdType identifierPDE = NO_PDE_ID );

    void GetSolutionVal( StdVector *sol, Complex* &ptSol,
                         const PdeIdType identifierPDE = NO_PDE_ID );

    void InvalidateSolBuffer() {
      solBufferIsValid_ = false;
    }

    void GetRHSVal( StdVector *rhs, Double* &ptRhs, 
                    const PdeIdType identifierPDE = NO_PDE_ID);

    void GetRHSVal( StdVector *rhs, Complex* &ptRhs, 
                    const PdeIdType identifierPDE = NO_PDE_ID);

    void InvalidateRhsBuffer() {
      rhsBufferIsValid_ = false;
    }

    void AddToDiagMatrixEntry( StdMatrix *stdMat, Integer eqnNum,
                               Integer dof, Double *val );

    void AdaptSystemMatrix( StdMatrix &stdMat, UInt *dirichletEQN,
                            UInt *dirichletComponent,
                            UInt numIDBC,
                            Double &penaltyTerm );

    void AdaptRHSForIDBC( StdVector &rhs,
                          UInt *dirichletEQN,
                          UInt *dirichletComponent,
                          StdVector &dirichletValue,
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
