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
    SingleEntryInt( PtrCoefFct& val );

    //! Copy constructor
    SingleEntryInt(const SingleEntryInt& right );

    //! \copydoc LinearForm::Clone
    virtual SingleEntryInt* Clone(){
      return new SingleEntryInt( *this );
    }

    //! Destructor
    virtual ~SingleEntryInt();

    //! Calculation of element 'vector' (real case )
    void CalcElemVector(Vector<Double>& elemVec,EntityIterator& ent);

    //! Calculation of element 'vector' (complex case )
    void CalcElemVector(Vector<Complex>& elemVec,EntityIterator& ent);

    //! \see LinearForm::IsComplex
    bool IsComplex() {
      return val_->IsComplex();
    }
    
  protected:

    //! Coefficient function
    PtrCoefFct val_;
  };
}

#endif
