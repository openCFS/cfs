// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINMASSINT_1
#define FILE_NLINMASSINT_1

#include "baseForm.hh"

namespace CoupledField
{

  //! Class for calculation  element mass matrix
  class NlinMassInt : public BaseForm
  {
  public:

    // Constructor
    NlinMassInt( Double density , const UInt nrDofsPerNode=1, 
                 bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~NlinMassInt();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
  
  protected:
    
  
  private:

    Double density_;          //!< multiplicative value for mass integrator
    UInt nrDofsPerNode_;   //!< degrees of freedom per node
  
  };





}

#endif // FILE_MASSINT
