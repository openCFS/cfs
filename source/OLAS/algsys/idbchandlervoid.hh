#ifndef IDBC_HANDLER_VOID_HH
#define IDBC_HANDLER_VOID_HH


#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "algsys/baseidbchandler.hh"


namespace OLAS {


  //! Null version of the IDBC_Handler concept

  //! This is class is a sort of null version of the %IDBC_Handler concept.
  //! It implements all methods required by the BaseIDBC_Handler interface
  //! in the form of stubs. It is intended for those cases in which no
  //! inhomogeneous Dirichlet boundary conditions are present in the linear
  //! system.
  class IDBC_HandlerVoid : public BaseIDBC_Handler {

  public:

    IDBC_HandlerVoid() {};

    //! Destructor
    virtual ~IDBC_HandlerVoid() {};

    //! Combine different FE matrices into a single system matrix
    void BuiltSystemMatrix( const std::map<FEMatrixType,
                            Double> &factors ) {};

    //! Adapt system matrix
    void AdaptSystemMatrix( BaseMatrix &sysMat ) {};

    //! Incorporate inhomogeneous Dirichlet BCs into right hand side
    void AddIDBCToRHS( BaseVector *rhs ) {};

    //! Remove inhomogeneous Dirichlet BCs from right hand side
    void RemoveIDBCFromRHS( BaseVector *rhs ) {};

    //! Set value for a Dirichlet boundary condition
    void SetIDBC( PdeIdType pdeID, UInt eqnNo, UInt comp, UInt bcNum,
                  const Double &val ) {};

    //! Set value for a Dirichlet boundary condition
    void SetIDBC( PdeIdType pdeID, UInt eqnNo, UInt comp, UInt bcNum,
                  const Complex &val ) {};

    //! Re-set specified internal matrix to zero
    void InitMatrix( FEMatrixType matrixID ) {};

    //! Re-set vector of Dirichlet values
    void InitDirichletValues() {};

    //! Set fixed dofs to specified Dirichlet boundary values
    void SetDofsToIDBC( BaseVector *vec ) {};

    //! Add weight of coupling between a fixed and a free dof into matrix
    void AddWeightFixedToFree( FEMatrixType matID,
                               PdeIdType pdeID1,
                               PdeIdType pdeID2,
                               UInt rowInd,
                               UInt colInd,
                               Double realPart,
                               Double imagPart = 0.0 ) {};

    //! Set weight of coupling between a fixed and a free dof into matrix
    void SetWeightFixedToFree( FEMatrixType matID,
                               PdeIdType pdeID1,
                               PdeIdType pdeID2,
                               UInt rowInd,
                               UInt colInd,
                               Double realPart,
                               Double imagPart = 0.0 ) {};

     //! Get weight of coupling between a fixed and a free dof from matrix
    void GetWeightFixedToFree( FEMatrixType matID, PdeIdType pdeID1,
                               PdeIdType pdeID2, UInt rowInd, UInt colInd,
                               Double & realPart,Double & imagPar ) const  {};
    
    //! Set the value of all coupling weights of a free dof to its fixed ones
    void SetRowWeights( FEMatrixType matID, PdeIdType pdeID, UInt rowInd,
                        Double realPart, Double imagPart = 0.0 ) {};
    
    
    //! Set the value of all coupling weights of a fixed dof to its free ones
    void SetColWeights( FEMatrixType matID, PdeIdType pdeID,UInt colInd,
                        Double realPart, Double imagPart = 0.0 ) {};
    
  };

}

#endif
