#ifndef BASE_ILD_FACTORISER_HH
#define BASE_ILD_FACTORISER_HH

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  template<typename> class SCRS_Matrix;
  class ParamNode;

  //! Base class for all incomplete LDL factorisation variants

  //! This class acts as base class for all variants of incomplete \f$LDL^T\f$
  //! factorisation variants.
  template <class T>
  class BaseILDLFactoriser {

  public:

    //! Default constructor
    BaseILDLFactoriser() {};

    //! Default destructor
    virtual ~BaseILDLFactoriser() {};

    //! Incomplete \f$LDL^T\f$ factorisation of a square matrix

    //! This is the central method of the class. It computes an
    //! incomplete \f$LDL^T\f$ factorisation of a square matrix. It must be
    //! over-written by the derived class, which implements its special
    //! version of the factorisation, i.e. a special dropping strategy.
    virtual void Factorise( SCRS_Matrix<T> &sysMat, std::vector<T> &dataD,
                            std::vector<UInt> &rptrU, std::vector<UInt> &cidxU,
                            std::vector<T> &dataU, bool newPattern ) = 0;

  protected:

    //! Pointer to parameter object

    //! This is a pointer to a parameter object containing the steering
    //! parameters for this preconditioner.
    PtrParamNode xml_;

    //! Pointer to report object

    //! This is a pointer to a report object which the preconditioner can use
    //! to store general information about its performance or setup phase.
    PtrParamNode olasInfo_;

    //! Dimension of problem matrix which is factorised
    UInt sysMatDim_;

    //! Status of factorisation
    bool amFactorised_;
    
    //! Flag for logging
    bool logging_;
    
  };

}

#endif
