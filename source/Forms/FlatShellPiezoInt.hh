// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FLATSHELLPIEZOINT
#define FILE_FLATSHELLPIEZOINT

#include "FlatShellInt.hh"

namespace CoupledField
{
  class FlatShellStiffInt;
     	
  //! Class for flat shell Integrator
  class FlatShellPiezoInt : public FlatShellInt
  {
  public:
    
    //! Constructor with pointer to composite material
    FlatShellPiezoInt( Composite *composite);
   
    //!  Destructor
    virtual ~FlatShellPiezoInt();

    //!returns nr. of degrees of freedom
    virtual Integer getNrDofs(){return 6;};
    
    //!returns dimension of D matrix of bending stiffness
    virtual Integer getDimD(){return 8;};
   				    
    //!Function for calculation BE matrix 
    void CalcElementMatrix( Matrix<Double>& elemMat,
			    EntityIterator& ent1, 
			    EntityIterator& ent2 );

  protected:
           
    //!returns B - matrix for BtE
    void calcBMat(Matrix<Double> & bMat, Integer ip, Matrix<Double> & ptCoord );
    
    //!returns E - matrix for BE multiplication
    void calcEMat(Matrix<Double> & eMat );  
 
    
    //Parametres of the lamina structures
    Vector<Double> z_;
    Vector<Double> orAngle_;
    
    //Pointer to FlatShellStiffInt
    FlatShellStiffInt * ptr_StiffInt_;   
  };

}//end of namespace

#endif 
