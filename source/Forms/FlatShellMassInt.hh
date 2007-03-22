// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FLATSHELLMASSINT_1
#define FILE_FLATSHELLMASSINT_1

#include "FlatShellInt.hh"
#include "Domain/Composite.hh"

namespace CoupledField
{

  //! Class for calculation  element mass matrix
  class FlatShellMassInt : public FlatShellInt
  {
  public:
    
    //!Constructor for normal Material
    FlatShellMassInt(BaseMaterial * matData);
        
    //Constructor for composite material
    FlatShellMassInt(Composite * composite);
    
    // Destructor
    virtual ~FlatShellMassInt();

    //! Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Double>& stiffMat,
			    EntityIterator& ent1, 
			    EntityIterator& ent2 );
     
     //!returns nr. of degrees of freedom
    virtual Integer getNrDofs(){return 6;};    
    
  protected:
    
  private:

    Double density_;          //!< multiplicative value for mass integrator
      
    //Parametres of the lamina structures
    Vector<Double> z_;
    Vector<Double> orAngle_;
  
  };

}

#endif // FILE_FLATSHELLMASSINT
