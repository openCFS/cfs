#ifndef CFS_SINGLEENTRY_INT_HH
#define CFS_SINGLEENTRY_INT_HH

#include "baseForm.hh"

namespace CoupledField {

  class SingleEntryInt : public BaseForm {

  public:

    //! Constructor
    SingleEntryInt( Double value, UInt dof, UInt numDofs );

    //! Destructor
    virtual ~SingleEntryInt();

    //! Calculation of element 'matrix'
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
  protected:

    //! factor to be set 
    Double entry_;

    //! dof to be set
    UInt dof_;

    //! maximum number of dofs
    UInt numDofs_;
  };
}

#endif
