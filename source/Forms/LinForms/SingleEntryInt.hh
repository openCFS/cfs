#ifndef CFS_SINGLEENTRY_INT_HH
#define CFS_SINGLEENTRY_INT_HH

#include "LinearForm.hh" 
#include "MatVec/Vector.hh"

namespace CoupledField {

  //! (Bi)Linearform, which sets just one single entry 
  class SingleEntryInt : public LinearForm {

  public:

    //! Constructor
    //! \param val Coefficient function (vector valued)
    SingleEntryInt( shared_ptr<CoefFunction>& val );

    //! Destructor
    virtual ~SingleEntryInt();

    //! Calculation of element 'vector' (real case )
    void CalcElemVector(Vector<Double>& elemVec,EntityIterator& ent);

    //! Calculation of element 'vector' (complex case )
    void CalcElemVector(Vector<Complex>& elemVec,EntityIterator& ent);

    //! \copydoc LinearForm::IsSolDependent
    bool IsSolDependent() {
      return val_->GetDependency() == CoefFunction::SOLUTION;
    }
    
    //! \see LinearForm::IsComplex
    bool IsComplex() {
      return val_->IsComplex();
    }
    
  protected:

    //! Coefficient function
    shared_ptr<CoefFunction> val_;
  };
}

#endif
