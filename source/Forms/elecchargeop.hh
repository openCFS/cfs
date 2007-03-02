#ifndef FILE_ELECCHARGEDOP_2004
#define FILE_ELECCHARGEDOP_2004

#include "Forms/baseoperator.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

#include "olas.hh"

namespace CoupledField {

  class Grid;
  class Elem;
  template<class TYPE> class ElemStoreSol;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;
 
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
