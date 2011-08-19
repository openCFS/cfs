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
    //! \param imp_magn math parser expression for magnitude of the impedance
    //! \param imp_phase math parser expression for phase of the impedance
    //! \param isAxi flag indicating axisymmetric geometry
    AbsorbingBCsInt( bool optimalImpedance,
                     const std::string &imp_magn,
                     const std::string &imp_phase,
                     bool isAxi );

    //! Default destructor
    virtual ~AbsorbingBCsInt();
    
    //! Calculate elementwise matrix
    //! This method calculates the element matrix for absorbing boundary
    //! conditions of first order. 
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

    //! Calculate elementwise matrix
    //! This method calculates the element matrix for absorbing boundary
    //! conditions of first order. 
    void CalcElementMatrix( Matrix<Complex>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );

  private:
    
    //! should the impedance be absorbing optimally (i.e. ABC)?
    bool optImp_;
    
    //! stores the magnitude of the impedance
    std::string impedMagn_;
    
    //! stores the phase of the impedance
    std::string impedPhase_;
    
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
