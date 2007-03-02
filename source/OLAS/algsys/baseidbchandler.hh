#ifndef BASE_IDBC_HANDLER_HH
#define BASE_IDBC_HANDLER_HH


#include <set>
#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "algsys/baseidbchandler.hh"


namespace OLAS {


  //! Base class for all IDBC_Handler classes

  //! This is the base class for all IDBC_Handler classes. The only current
  //! raison d'etre for this base class is that we need it for being able to
  //! manage Double and Complex instances of the IDBC_Handler class, which
  //! is the only child, in the algebraic system transparently.
  class BaseIDBC_Handler {

  public:

    //! Destructor
    virtual ~BaseIDBC_Handler() {
      ENTER_FCN( "BaseIDBC_Handler::~BaseIDBC_Handler" );
    }

    //! Combine different FE matrices into a single system matrix

    //! This method must be called to combine the different internal FE
    //! matrices into a single system matrix. The system matrix is constructed
    //! by computing a weighted sum of the internal FE matrices using the
    //! weights supplied in the factors map input parameter.
    virtual void
    BuiltSystemMatrix( const std::map<FEMatrixType, Double> &factors ) = 0;

    //! Adapt system matrix

    //! This method can be used by those approaches that require a
    //! modification of the system matrix in order to incorporate
    //! inhomogeneous Dirichlet boundary conditions into the linear system,
    //! like e.g. the penalty approach.
    virtual void AdaptSystemMatrix( BaseMatrix &sysMat ) = 0;

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
    virtual void AddIDBCToRHS( BaseVector *rhs ) = 0;

    //! Remove inhomogeneous Dirichlet BCs from right hand side

    //! Call this method in order to update an exisiting right hand side
    //! such that it incorporates the modification resulting from eliminating
    //! the dofs for inhomogeneous Dirichlet boundary conditions from the
    //! linear system.
    //! \note It is the caller's responsibility to ensure that the initial
    //!       right-hand side passed to the method does not contain Dirichlet
    //!       values from previous calls, since the method cannot eliminate
    //!       them.
    //! \param rhs vector with right-hand side entries
    virtual void RemoveIDBCFromRHS( BaseVector *rhs ) = 0;

    //@{
    //! Set value for a Dirichlet boundary condition

    //! This method can be used to set the value of a degree of freedom that
    //! is fixed by a homogeneous Dirichlet boundary condition.
    //! \param pdeID identifier of PDE; this is only used in the case of an
    //!              SBM_System in order to identify the sub-vector in which
    //!              the value must be stored
    //! \param eqnNo equation number for the degree of freedom whose value
    //!              should be set
    //! \param comp  identifies the block component of the fixed degree of
    //!              freedom in the case of dof blocking
    //! \param bcNum specifies the number of this Dirichlet value in the
    //!              range [1,..., no. inhom. Dirchlet BCs]
    //! \param val   inhomogeneous Dirichlet value
    virtual void SetIDBC( PdeIdType pdeID, UInt eqnNo, UInt comp, UInt bcNum,
                          const Double &val ) {
      (*error) << "BaseIDBC_Handler::SetIDBC: The derived class does "
               << "obviously not support the Double version of this "
               << "interface! So it is probably a Complex instance!";
      Error( __FILE__, __LINE__ );
    }

    virtual void SetIDBC( PdeIdType pdeID, UInt eqnNo, UInt comp, UInt bcNum,
                          const Complex &val )  {
      (*error) << "BaseIDBC_Handler::SetIDBC: The derived class does "
               << "obviously not support the Complex version of this "
               << "interface! So it is probably a Double instance!";
      Error( __FILE__, __LINE__ );
    }
    //@}

    //! Add weight of coupling between a fixed and a free dof into matrix

    //! This method provides an interface to add the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition into the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxilliary matrix
    //! \param pdeID1   row index of sub-matrix in the SBM_Matrix case
    //! \param pdeID2   column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to be set, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to be set, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param realPart real valued part of the weight of the coupling
    //! \param imagPart imaginary part of the weight of the coupling
    virtual void AddWeightFixedToFree( FEMatrixType matID,
                                       PdeIdType pdeID1,
                                       PdeIdType pdeID2,
                                       UInt rowInd,
                                       UInt colInd,
                                       Double realPart,
                                       Double imagPart = 0.0 ) = 0;

    //! Set weight of coupling between a fixed and a free dof into matrix

    //! This method provides an interface to set the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition into the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxilliary matrix
    //! \param pdeID1   row index of sub-matrix in the SBM_Matrix case
    //! \param pdeID2   column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to be set, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to be set, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param realPart real valued part of the weight of the coupling
    //! \param imagPart imaginary part of the weight of the coupling
    virtual void SetWeightFixedToFree( FEMatrixType matID,
                                       PdeIdType pdeID1,
                                       PdeIdType pdeID2,
                                       UInt rowInd,
                                       UInt colInd,
                                       Double realPart,
                                       Double imagPart = 0.0 ) = 0;

    //! Get weight of coupling between a fixed and a free dof from matrix

    //! This method provides an interface to get the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition from the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxiliary matrix
    //! \param pdeID1   row index of sub-matrix in the SBM_Matrix case
    //! \param pdeID2   column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to get, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to get, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param realPart real valued part of the weight of the coupling
    //! \param imagPart imaginary part of the weight of the coupling
    virtual void GetWeightFixedToFree( FEMatrixType matID,
                                       PdeIdType pdeID1,
                                       PdeIdType pdeID2,
                                       UInt rowInd,
                                       UInt colInd,
                                       Double & realPart,
                                       Double & imagPart ) const = 0;

    //! Re-set specified internal matrix to zero

    //! Calling this method re-sets the specified internal matrix of the
    //! object to its initial state, i.e. Init() is called on that matrix
    //! object.
    virtual void InitMatrix( FEMatrixType matrixID ) = 0;

    //! Re-set vector of Dirichlet values

    //! Calling this method deletes all information stored by the object
    //! internally on the values and degrees of freedom that are fixed
    //! by inhomogeneous Dirichlet boundary conditions.
    virtual void InitDirichletValues() = 0;

    //! Set fixed dofs to specified Dirichlet boundary values

    //! This method replaces the values of all fixed degrees of freedom in the
    //! specified input vector by new values. These new values are taken to
    //! be the values specified via the inhomogeneous Dirichlet boundary
    //! condition that fixes the respective degrre of freedom.
    virtual void SetDofsToIDBC( BaseVector *vec ) = 0;

  };

}

#endif
