// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOUPOWERDENSITYOP
#define FILE_ACOUPOWERDENSITYOP

#include <complex>

#include "Forms/baseoperator.hh"
#include "General/defs.hh"
#include "MatVec/exprt/xpr2.hh"



namespace CoupledField {

class EntityIterator;
class EqnMap;
  // Forward declaration of classes
  class Grid;
class StdPDE;
  template<class TYPE> class Vector;
  
  
  //! This operator class calculates the acoustic power density per
  //! element for a acoustic PDE formulation in velocity potential

  template<class TYPE>
  class AcouPowerDensityOp : public BaseOperator
  {

  
  public:

    //! Constructor
    /*!
      \param ptGrid (input) Pointer to grid
      \param ptPDE  (input) Pointer to PDE
      \param ptEQN  (input) Pointer to EQN
      \param isaxi  (input) Flag for axi-symmetric geomtetry
    */
    AcouPowerDensityOp(Grid * ptGrid,
                       StdPDE * ptPDE,
                       shared_ptr<EqnMap> eqnMap,
                       bool isaxi=false);
    
    //! Destructor
    virtual ~AcouPowerDensityOp();
  
    //! Calculate acoustic power density for element
    /*!
      \param elemPD (output) Element vector of acoustic power density
      \param ptElem (input) Pointer to element
      \param density (input) multiplicative factor
    */
    virtual void CalcElemPD(Vector<Double> & elemPD,
                            const EntityIterator& it,
                            const Double density);
    

//     //! Calculate acoustic power density for elements in subdomains
//     /*!
//       \param elemPD (output) Vector containing
//       \param SD (input) Identifier of the subdomain
//       \param density (input) multiplicative factor
//     */
//     virtual void CalcSubdomPD(Vector<Double> & elemPD,
//                               const StdVector<RegionIdType> & SD,
//                               const Double density);
                                                       

  protected:
  
    bool isaxi_;
    Double dimensions_;

  private:

    //! Calculate first product for transient(Double) and harmonic(Complex) case
    Double ComputeN1( Vector<TYPE> solGrAtIp,
                      Vector<TYPE> solD1GrAtIp );

    //! Calculate second product for transient(Double) and harmonic(Complex) case
    Vector<Double> ComputeN2( Vector<TYPE> solGrAtIp,
                              TYPE solD1AtIp );
  };


} // end of namespace

#endif
