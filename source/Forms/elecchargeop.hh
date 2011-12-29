// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ELECCHARGEOP_2004
#define FILE_ELECCHARGEOP_2004

#include <complex>

#include "Forms/baseoperator.hh"
#include "General/defs.hh"
#include "MatVec/exprt/xpr2.hh"

namespace CoupledField {

class EntityIterator;
class EqnMap;
  class Grid;
class StdPDE;
  class SurfElemList;
struct ResultInfo;
  template<class TYPE> class Vector;
 
  //! Operator for calculating the electric charge on surface elements
  
  //! This operator class calculates electric charges on surface
  //! elements by evaluating a surface integral of the electric
  //! flux density
   

  template<class TYPE>
  class ElecChargeOp : public BaseOperator {

  
  public:

    //! Constructor
    /*!
      \param ptGrid (input) Pointer to grid
      \param ptPDE (input) Pointer to PDE
      \param ptEQN (input) Pointer to EQN
      \param isaxi (input) Flag for axi-symmetric geomtetry
    
    */
    ElecChargeOp(Grid * ptGrid,
                 StdPDE * ptPDE,
                 shared_ptr<EqnMap> eqnMap,
                 shared_ptr<ResultInfo> result,
                 bool isaxi=false,
                 bool coordUpdate = false );

    //! Destructor
    virtual ~ElecChargeOp();
  
    //! Calculate charge for one surfac element
    /*!
      \param charge (output) Charge of given surface element
      \param ent (input) EntityIterator pointing to current element
      \param lCoord (input) Local coordinates of evaluation point
      \param eFluxDensity (input) Normal component of Flux density in lCoord
    */
    virtual void CalcElemCharge(TYPE & charge,
                                const EntityIterator& ent,
                                const Vector<Double> & lCoord,
                                const TYPE & eNormalFluxDensity);
 

    //! Calculate charges for whole surface
    /*!
      \param charges (output) Vector of element charges
      \param surfElems (input) Entity list with surface elements
      \param lCoord (input) Local coordinate of evaluation
      \param eFluxDensity (input) Vector of normal components of 
      Flux density in lCoord
    */
    virtual void CalcElemCharges(Vector<TYPE> & charges,
                                 const shared_ptr<SurfElemList> surfElems,
                                 const Vector<Double> & lCoord,
                                 const Vector<TYPE> & eNormalFluxDensity);
  protected:

    //! Pointer to resultDof object
    shared_ptr<ResultInfo> result_;
  };

} // end of namespace

#endif
