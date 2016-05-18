#ifndef FILE_RHS_KX_INTEGRATOR_
#define FILE_RHS_KX_INTEGRATOR_

#include "LinearForm.hh"
#include <boost/tr1/type_traits.hpp>
#include "Domain/CoefFunction/CoefFunction.hh"


namespace CoupledField{


  //! Calculates product of element matrix K and element solution X 

  //! This class implements the general integrator for RHS integrators of 
  //! the form
  //!   \[ {\mathbf K} \cdot \vec{x}  \]
  //! So we have a quantity X specified by a FeFunction
  //! passed to the constructor as well as a bilinearform K. 
  template< class VEC_DATA_TYPE=Double>
  class KXIntegrator : public LinearForm{
  public:
    
    //! Constructor
    KXIntegrator( BiLinearForm * biLinForm,
                  VEC_DATA_TYPE factor,
                  shared_ptr<BaseFeFunction> feFct );

    //! Copy constructor
    KXIntegrator(const KXIntegrator& right );

    //! \copydoc LinearForm::Clone
    virtual KXIntegrator* Clone(){
      return new KXIntegrator( *this );
    }

    //! Destructor 
    ~KXIntegrator(){
    }

    //! Calculate element vector
    void CalcElemVector(Vector<VEC_DATA_TYPE> & elemVec,EntityIterator& ent);

    //! Calculate element vector
    void CalcRhsVector(Vector<VEC_DATA_TYPE> & elemVec,shared_ptr<CoefFunction > rhsCoefs,EntityIterator& ent);

    //! \copydoc LinearForm::IsSolDependent
    bool IsSolDependent() {
      return true;
    }
    
    //! Return if linearform is complex
    bool IsComplex(){
      return std::tr1::is_same<VEC_DATA_TYPE,Complex>::value;
    }

  protected:

    //! Bilinearform
    BiLinearForm * form_;

    //! Additional scalar factor
    VEC_DATA_TYPE factor_;

    //! FeFunction of the solution vector X
    shared_ptr<FeFunction<VEC_DATA_TYPE> > feFct_;
  };
}
//Include template definition file
#include "KXInt.cc"
#endif
