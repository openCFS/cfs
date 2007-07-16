// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_SINGLEENTRY_INT_HH
#define CFS_SINGLEENTRY_INT_HH

#include "baseForm.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {

  class SingleEntryInt : public BaseForm {

  public:

    //! Constructor (real case)
    SingleEntryInt( const std::string& real,
                    UInt dof, UInt numDofs );

    //! Constructor (complex case )
    SingleEntryInt( const std::string& real,
                    const std::string& imag,
                    UInt dof, UInt numDofs );

    //! Destructor
    virtual ~SingleEntryInt();

    //! Calculation of element 'matrix' (real case )
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

    //! Calculation of element 'matrix' (complex case )
    void CalcElementMatrix( Matrix<Complex>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
  protected:

    //! factors to be set 
    std::string real, imag;

    //! handles for math parser
    MathParser::HandleType rHandle_, iHandle_;

    //! dof to be set
    UInt dof_;

    //! maximum number of dofs
    UInt numDofs_;
  };
}

#endif
