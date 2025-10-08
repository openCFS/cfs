#ifndef IDBC_HANDLER_VOID_HH
#define IDBC_HANDLER_VOID_HH


#include "OLAS/algsys/BaseIDBC_Handler.hh"


namespace CoupledField {


  //! Null version of the IDBC_Handler concept

  //! This class is a sort of null version of the %IDBC_Handler concept.
  //! It implements all methods required by the BaseIDBC_Handler interface
  //! in the form of stubs. It is intended for those cases in which no
  //! inhomogeneous Dirichlet boundary conditions are present in the linear
  //! system.
  class IDBC_HandlerVoid : public BaseIDBC_Handler {

  public:

    //! Default constructor
    IDBC_HandlerVoid() {};

    //! Destructor
    virtual ~IDBC_HandlerVoid() {};

    //! @copydoc BaseIDBC_Handler::BuiltSystemMatrix
    void BuildSystemMatrix( const std::map<FEMatrixType,
                            Double> &factors,
                            std::map<UInt, std::set<UInt> >& indPerBlock) {};
    void BuildSystemMatrix( const std::map<FEMatrixType,
                            Complex> &factors,
                            std::map<UInt, std::set<UInt> >& indPerBlock) {};

    //! @copydoc BaseIDBC_Handler::AdaptSystemMatrix()
    void AdaptSystemMatrix( SBM_Matrix &sysMat ) {};

    //! @copydoc BaseIDBC_Handler::AddIDBCToRHS()
    void AddIDBCToRHS( SBM_Vector *rhs, bool deltaIDBC = false ) {};

    //! @copydoc BaseIDBC_Handler::RemoveIDBCFromRHS()
    void RemoveIDBCFromRHS( SBM_Vector *rhs, bool deltaIDBC = false ) {};

    //! @copydoc BaseIDBC_Handler::SetIDBC()
    void SetIDBC( UInt rowBlock, UInt rowNum, const Double &val ) {};

    //! @copydoc BaseIDBC_Handler::SetIDBC()
    void SetIDBC( UInt rowBlock, UInt rowNum, const Complex &val ) {};
    
    //! @copydoc BaseIDBC_Handler::GetIDBC()
    void GetIDBC( UInt rowBlock, UInt rowNum, Double &val ) {};

    //! @copydoc BaseIDBC_Handler::GetIDBC()
    void GetIDBC( UInt rowBlock, UInt rowNum, Complex &val ) {};

    //! @copydoc BaseIDBC_Handler::InitMatrix()
    void InitMatrix( FEMatrixType matrixType ) {};

    //! @copydoc BaseIDBC_Handler::InitDirichletValues()
    void InitDirichletValues() {};

    void SetOldDirichletValues() {};

    //! @copydoc BaseIDBC_Handler::SetDofsToIDBC()
    void SetDofsToIDBC( SBM_Vector *vec, bool deltaIDBC = false ) {};

    //! @copydoc BaseIDBC_Handler::AddWeightFixedToFree()
    void AddWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const Double& val ) {};

    //! @copydoc BaseIDBC_Handler::AddWeightFixedToFree()
    void AddWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const Complex& val ) {};


    //! @copydoc BaseIDBC_Handler::SetWeightFixedToFree()
    void SetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const Double& val ) {};

    //! @copydoc BaseIDBC_Handler::SetWeightFixedToFree()
    void SetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const Complex& val ) {};

    //! @copydoc BaseIDBC_Handler::GetWeightFixedToFree()
    void GetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               Double & val )  {};

    //! @copydoc BaseIDBC_Handler::GetWeightFixedToFree()
    void GetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               Complex & val )  {};

    void PrintIDBCvec(){
      EXCEPTION("Not implemented here");
    }

  };

}

#endif
