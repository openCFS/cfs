#ifndef BASE_IDBC_HANDLER_HH
#define BASE_IDBC_HANDLER_HH

#include <set>
#include <map>

#include "General/Environment.hh"
#include "General/defs.hh"

#include "MatVec/SBM_Vector.hh"
#include "MatVec/SBM_Matrix.hh"

#include "OLAS/algsys/BaseIDBC_Handler.hh"

namespace CoupledField {


  //! Base class for all IDBC_Handler classes

  //! This is the base class for all IDBC_Handler classes. The only current
  //! raison d'etre for this base class is that we need it for being able to
  //! manage Double and Complex instances of the IDBC_Handler class, which
  //! is the only child, in the algebraic system transparently.
  class BaseIDBC_Handler {

  public:

    //! Destructor
    virtual ~BaseIDBC_Handler() {
    }

    //! Combine different FE matrices into a single system matrix

    //! This method must be called to combine the different internal FE
    //! matrices into a single system matrix. The system matrix is constructed
    //! by computing a weighted sum of the internal FE matrices using the
    //! weights supplied in the factors map input parameter.
    //! The operation only works on a subset of the free equations, denoted
    //! by the indicesPerBlock (i.e. for each SBM matrix (first map) there
    //! is a set of row-indices (second index ) ).
    virtual void
    BuildSystemMatrix( const std::map<FEMatrixType, Double> &factors,
                       std::map<UInt, std::set<UInt> >& indicesPerBlock ) = 0;
    //! OVERLOAD: BuildSystemMatrix with complex factors
    virtual void
    BuildSystemMatrix( const std::map<FEMatrixType, Complex> &factors, std::map<UInt, std::set<UInt> >& indicesPerBlock ) = 0;
    //! Adapt system matrix

    //! This method can be used by those approaches that require a
    //! modification of the system matrix in order to incorporate
    //! inhomogeneous Dirichlet boundary conditions into the linear system,
    //! like e.g. the penalty approach.
    virtual void AdaptSystemMatrix( SBM_Matrix &sysMat ) = 0;

    //! Incorporate inhomogeneous Dirichlet BCs into right hand side

    //! Call this method in order to update an exsisting right hand side
    //! such that it incorporates the modification resulting from eliminating
    //! the dofs for inhomogeneous Dirichlet boundary conditions from the
    //! linear system.
    //! \note It is the caller's responsibility to ensure that the initial
    //!       right-hand side passed to the method does not contain Dirichlet
    //!       values from previous calls, since the method will not eliminate
    //!       them.
    //! \param rhs vector with right-hand side entries
    //!         deltaIDBC: if true, add idbc_values - oldIdbc_values
    virtual void AddIDBCToRHS( SBM_Vector *rhs, bool deltaIDBC = false ) = 0;

    virtual void SetOldDirichletValues(){};

    virtual void ToString(){};

    //! Remove inhomogeneous Dirichlet BCs from right hand side

    //! Call this method in order to update an existing right hand side
    //! such that it incorporates the modification resulting from eliminating
    //! the dofs for inhomogeneous Dirichlet boundary conditions from the
    //! linear system.
    //! \note It is the caller's responsibility to ensure that the initial
    //!       right-hand side passed to the method does not contain Dirichlet
    //!       values from previous calls, since the method cannot eliminate
    //!       them.
    //! \param rhs vector with right-hand side entries
    virtual void RemoveIDBCFromRHS( SBM_Vector *rhs, bool deltaIDBC = false ) = 0;
    //@{
    //! Set value for a Dirichlet boundary condition

    //! This method can be used to set the value of a degree of freedom that
    //! is fixed by an inhomogeneous Dirichlet boundary condition.
    //! \param rowBlock sbm row Number
    //! \param rowNum number of the row for the degree of freedom whose value
    //!               should be set
    //! \param val    inhomogeneous Dirichlet value
    virtual void SetIDBC( UInt rowBlock, UInt rowNum, const Double &val ) {
      EXCEPTION("BaseIDBC_Handler::SetIDBC: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");
    }

    virtual void SetIDBC( UInt rowBlock, UInt rowNum, const Complex &val )  {
      EXCEPTION("BaseIDBC_Handler::SetIDBC: The derived class does " \
                << "obviously not support the Complex version of this " \
                << "interface! So it is probably a Double instance!");
    }
    //@}
    
    //@{
    //! Get value for a Dirichlet boundary condition

    //! This method can be used to set the value of a degree of freedom that
    //! is fixed by an inhomogeneous Dirichlet boundary condition.
    //! \param rowBlock sbm row Number
    //! \param rowNum number of the row for the degree of freedom whose value
    //!               should be set
    //! \param val    inhomogeneous Dirichlet value
    virtual void GetIDBC( UInt rowBlock, UInt rowNum, Double &val, bool deltaIDBC=false ) {
      EXCEPTION("BaseIDBC_Handler::GetIDBC: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");
    }

    virtual void GetIDBC( UInt rowBlock, UInt rowNum, Complex &val, bool deltaIDBC=false )  {
      EXCEPTION("BaseIDBC_Handler::GetIDBC: The derived class does " \
                << "obviously not support the Complex version of this " \
                << "interface! So it is probably a Double instance!");
    }
    //@}

    //@{
    //! Add Dof Values fixed by Dirichlet boundary condition to a RHS

    //! In case of elemenination IDBC handling in the time domain we need to 
    //! add the DOFs fixed by an IDBC to the RHS by multiplication with the 
    //! corresponding matrix entries.
    //! We do this for each row block and for each rowblock in each row of 
    //! the SBM matrix.
    //! \param matID matrix type to be multiplied with the value
    //! \param colBlock The column block the value is associated with
    //! \param colInd The index inside the sbm subvector
    //! \param rhs  the right hand side we are operating on
    //! \param val  The value we need to multiply
    virtual void AddFixedToFreeRHS( FEMatrixType matID, UInt colBlock,
                                        UInt colInd, SBM_Vector *rhs, 
                                        const Double& val ) {
      EXCEPTION("BaseIDBC_Handler::AddFixedToFreeRHS: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");
    }

    virtual void AddFixedToFreeRHS( FEMatrixType matID, UInt colBlock,
                            UInt colInd, SBM_Vector *rhs, const Complex& val )  {
      EXCEPTION("BaseIDBC_Handler::AddFixedToFreeRHS: The derived class does " \
                << "obviously not support the Complex version of this " \
                << "interface! So it is probably a Double instance!");
    }
    //@}

    //@{
    //! Queue a fixed dof value for a batched fixed-to-free RHS update

    //! Batched variant of AddFixedToFreeRHS(): instead of performing one
    //! matrix-vector product per fixed dof, the values of all fixed dofs are
    //! first gathered with this method and then applied with a single
    //! matrix-vector product in FinishFixedToFreeRHS().
    //! \param colBlock The column block the value is associated with
    //! \param colInd The index inside the sbm subvector
    //! \param val  The value we need to multiply
    virtual void QueueFixedToFreeRHS( UInt colBlock, UInt colInd,
                                      const Double& val ) {
      EXCEPTION("BaseIDBC_Handler::QueueFixedToFreeRHS: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");
    }

    virtual void QueueFixedToFreeRHS( UInt colBlock, UInt colInd,
                                      const Complex& val ) {
      EXCEPTION("BaseIDBC_Handler::QueueFixedToFreeRHS: The derived class does " \
                << "obviously not support the Complex version of this " \
                << "interface! So it is probably a Double instance!");
    }
    //@}

    //! Apply all queued fixed dof values to the RHS with one multiplication

    //! Performs rhs += auxMat[matID] * queuedValues for all values gathered
    //! via QueueFixedToFreeRHS() since the last call of this method.
    //! The default implementation is a no-op: handlers without fixed dofs
    //! (e.g. IDBC_HandlerVoid) never get values queued, so there is nothing
    //! to apply.
    //! \param matID matrix type to be multiplied with the queued values
    //! \param rhs  the right hand side we are operating on
    virtual void FinishFixedToFreeRHS( FEMatrixType matID, SBM_Vector *rhs ) {
    }

    //@{
    //! Add weight of coupling between a fixed and a free dof into matrix

    //! This method provides an interface to add the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition into the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxilliary matrix
    //! \param rowBlock row index of sub-matrix in the SBM_Matrix case
    //! \param colBlock column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to be set, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to be set, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param val value of the weight of the coupling
    virtual void AddWeightFixedToFree( FEMatrixType matID,
                                       UInt rowBlock,
                                       UInt colBlock,
                                       UInt rowInd,
                                       UInt colInd,
                                       const Double& val ) {
      EXCEPTION("BaseIDBC_Handler::AddWeightFixedToFree: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");
    }
    
    virtual void AddWeightFixedToFree( FEMatrixType matID,
                                       UInt rowBlock,
                                       UInt colBlock,
                                       UInt rowInd,
                                       UInt colInd,
                                       const Complex& val ) {
      EXCEPTION("BaseIDBC_Handler::AddWeightFixedToFree: The derived class does " \
                << "obviously not support the Complex version of this " \
                << "interface! So it is probably a Complex instance!");
    }
    //@}
    
    //@{
    //! Set weight of coupling between a fixed and a free dof into matrix

    //! This method provides an interface to set the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition into the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxilliary matrix
    //! \param rowBlock row index of sub-matrix in the SBM_Matrix case
    //! \param colBlock column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to be set, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to be set, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param val      value of the weight of the coupling
    virtual void SetWeightFixedToFree( FEMatrixType matID,
                                       UInt rowBlock,
                                       UInt colBlock,
                                       UInt rowInd,
                                       UInt colInd,
                                       const Double& val ) {
      EXCEPTION("BaseIDBC_Handler::SetWeightFixedToFree: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");
    }
    virtual void SetWeightFixedToFree( FEMatrixType matID,
                                       UInt rowBlock,
                                       UInt colBlock,
                                       UInt rowInd,
                                       UInt colInd,
                                       const Complex& val ) {
      EXCEPTION("BaseIDBC_Handler::SetWeightFixedToFree: The derived class does " \
                << "obviously not support the Complex version of this " \
                << "interface! So it is probably a Complex instance!");
    }
    //@}

    //@{
    //! Get weight of coupling between a fixed and a free dof from matrix

    //! This method provides an interface to get the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition from the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxiliary matrix
    //! \param rowBlock row index of sub-matrix in the SBM_Matrix case
    //! \param colBlock column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to get, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to get, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param val      value of the weight of the coupling
    virtual void GetWeightFixedToFree( FEMatrixType matID,
                                       UInt rowBlock,
                                       UInt colBlock,
                                       UInt rowInd,
                                       UInt colInd,
                                       Double & val ) {
      EXCEPTION("BaseIDBC_Handler::GetWeightFixedToFree: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");

    }
    virtual void GetWeightFixedToFree( FEMatrixType matID,
                                       UInt rowBlock,
                                       UInt colBlock,
                                       UInt rowInd,
                                       UInt colInd,
                                       Complex & val ) {
      EXCEPTION("BaseIDBC_Handler::GetWeightFixedToFree: The derived class does " \
                << "obviously not support the Double version of this " \
                << "interface! So it is probably a Complex instance!");

    }
    //@}

    //! Re-set specified internal matrix to zero

    //! Calling this method re-sets the specified internal matrix of the
    //! object to its initial state, i.e. Init() is called on that matrix
    //! object.
    virtual void InitMatrix( FEMatrixType matrixType ) = 0;

    //! Re-set vector of Dirichlet values

    //! Calling this method deletes all information stored by the object
    //! internally on the values and degrees of freedom that are fixed
    //! by inhomogeneous Dirichlet boundary conditions.
    virtual void InitDirichletValues() = 0;

    //! Set fixed dofs in given vector to specified Dirichlet boundary values

    //! This method replaces the values of all fixed degrees of freedom in the
    //! specified input vector by new values. These new values are taken to
    //! be the values specified via the inhomogeneous Dirichlet boundary
    //! condition that fixes the respective degrre of freedom.
    virtual void SetDofsToIDBC( SBM_Vector *vec, bool deltaIDBC = false ) = 0;

    //Just to print the IDBCvec
    virtual void PrintIDBCvec(){};

  };

}

#endif
