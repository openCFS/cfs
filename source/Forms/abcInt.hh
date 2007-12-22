// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef File_ABC_MECH_INT
#define File_ABC_MECH_INT

#include "baseForm.hh"

namespace CoupledField {

  //! Class for calculation surface integrals of absorbing boundary
  //! elements
  class AbsorbingBCsInt : public SurfForm
  {
  public:
    
    //! Default constructor
    //! \param isAxi flag indicating axi-symmetrix geometry
    AbsorbingBCsInt( bool isAxi );

    //! Default destructor
    virtual ~AbsorbingBCsInt();
    
    //! Sets a multiplicative factor for element matrix
    void SetFactor( const std::string& factor );

    //! Calculate elementwise matrix
    //! This method calculates the element matrix for absorbing boundary
    //! conditions of first order. 
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

  private:
    
  };
  
  
  

  //! heat flux term according to Newton's law of cooling
  class HeatFluxInt : public SurfForm
  {
  public:
    //! Default constructor
    //! \param factor multiplicative factor
    //! \param isaxi flag indicating axi-symmetrix geometry
    HeatFluxInt( const std::string& factor, bool isaxi );
    
    virtual ~HeatFluxInt();
    
    void SetFactor( const std::string& factor ) {
      mParser_->SetExpr( mHandle_, factor );
    }
    
    //! Calculate elementwise matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
  private:
    //
    
  };


}

#endif
