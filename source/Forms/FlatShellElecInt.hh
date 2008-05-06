// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FLATSHELLELECINT
#define FILE_FLATSHELLELECINT

#include "FlatShellInt.hh"

namespace CoupledField
{
     //! Class for flat shell Integrator
    class FlatShellElecInt : public FlatShellInt
 {
    public:
    
   //! Constructor with pointer to composite material
   FlatShellElecInt( Composite *composite,  bool hasDrillDof );
   
   //! Destructor
   virtual ~FlatShellElecInt();
   
   //! Function for calculation BE matrix 
   void CalcElementMatrix( Matrix<Double>& elemMat,
			   EntityIterator& ent1, 
			   EntityIterator& ent2 );
   
  protected:    
       
    //!returns E - matrix for BE multiplication
    void calcEMat(Matrix<Double> & eMat );  
 
    
    //Parametres of the lamina structures
    Vector<Double> z_;
    Vector<Double> orAngle_;
 };

}//end of namespace

#endif 
