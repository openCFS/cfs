//==============================================================================
/*!
 *       \file     CoefFunctionMeanFlowConvection.hh
 *       \brief    Coefficient function which compute sthe stabilization 
 *                 parameters for SUPG
 *
 *       \date     07/04/2013
 *       \author   Simon Triebenbacher
 */
//==============================================================================

#ifndef FILE_COEFFUNCTION_MEANFLOWCONVECTION_HH
#define FILE_COEFFUNCTION_MEANFLOWCONVECTION_HH

#include "CoefFunction.hh"
#include "CoefFunctionMulti.hh"
#include "FeBasis/FeFunctions.hh"
#include "Forms/Operators/BaseBOperator.hh"

namespace CoupledField {

  // forward class declaration
  class BaseBOperator;
  
  // ============================================================================
  //  Stabilization parameters for SUPG
  // ============================================================================
  //! Computes the stabilization parameters for SUPG


  template<typename T, UInt DOFS = 2>
  class CoefFunctionMeanFlowConvection : public CoefFunction
  {
  public:
    
    //! Constructor for tensor case: scalar * Grad(v0)
    CoefFunctionMeanFlowConvection(PtrCoefFct scalar,
                                   BaseBOperator* optGrad,
                                   shared_ptr<BaseFeFunction> feFncV0);

    //! Constructor for Vector-case: v0 . Grad(v0)
    CoefFunctionMeanFlowConvection(BaseBOperator* optGrad,
                                   shared_ptr<BaseFeFunction> feFncV0);

    //! Destructor
    virtual ~CoefFunctionMeanFlowConvection(){;}

    //! Return real-valued tensor at integration point
    virtual void GetTensor( Matrix<T>& tensor, 
                            const LocPointMapped& lpm );

    //! Return real-valued tensor at integration point
    virtual void GetVector( Vector<T>& vector, const LocPointMapped& lpm );

    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      numRows = dim_;
      numCols = dim_;
    }
    virtual UInt GetVecSize() const {
          return dim_;
        }

  protected:
    
    //! density
    PtrCoefFct factor_;

    BaseBOperator* bOperator_;
    
    //! Depending FeFunction
    shared_ptr<BaseFeFunction>  feFct_;
    
    //! dimension of the tensor/vector
    UInt dim_;
  };
}

#endif
