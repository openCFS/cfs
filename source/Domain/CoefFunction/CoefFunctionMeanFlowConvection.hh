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

#include <boost/tr1/type_traits.hpp>
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


  template<typename T>
  class CoefFunctionMeanFlowConvection : public CoefFunction
  {
  public:
    
    //! Constructor
    CoefFunctionMeanFlowConvection(PtrCoefFct density,
                                   PtrCoefFct viscosity,
                                   BaseBOperator* opt,
                                   shared_ptr<BaseFeFunction> feFnc);
    //! Destructor
    virtual ~CoefFunctionMeanFlowConvection(){;}
    
    //! Return real-valued tensor at integration point
    virtual void GetTensor( Matrix<T>& tensor, 
                            const LocPointMapped& lpm );
    
  protected:
    
    //! density
    PtrCoefFct density_;
    
    //! viscosity
    PtrCoefFct viscosity_;

    BaseBOperator* bOperator_;
    
    //! Depending FeFunction
    shared_ptr<BaseFeFunction>  feFct_;
    
  };
}

#endif
