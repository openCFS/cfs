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
   


  class ElecChargeOp : public BaseOperator {

  
  public:

    //! Constructor
    /*!
      \param ptGrid (input) Pointer to grid
      \param ptPDE (input) Pointer to PDE
      \param ptEQN (input) Pointer to EQN
      \param level (input) Multigrid level
      \param isaxi (input) Flag for axi-symmetric geomtetry
    
    */
    ElecChargeOp(Grid * ptGrid,
		 StdPDE * ptPDE,
		 NodeEQN * ptEQN,
		 const Integer level,
		 Boolean isaxi=FALSE);

    //! Destructor
    virtual ~ElecChargeOp();
  
    //! Calculate charge for one surfac element
    /*!
      \param charge (output) Charge of given surface element
      \param ptElem (input) Pointer to elemeent
      \param lCoord (input) Local coordinates of evaluation point
      \param eFluxDensity (input) Normal component of Flux density in lCoord
    */
    virtual void CalcElemCharge(Double & charge,
				const Elem * ptElement,
				const Vector<Double> & lCoord,
				const Double & eNormalFluxDensity);
 
    virtual void CalcElemCharge(Complex & charge,
				const Elem * ptElement,
				const Vector<Double> & lCoord,
				const Complex & eNormalFluxDensity);
  

    //! Calculate charges for whole surface
    /*!
      \param charges (output) Vector of element charges
      \param surfElems (input) Vector of surface elements
      \param lCoord (input) Local coordinate of evaluation
      \param eFluxDensity (input) Vector of normal components of 
      Flux density in lCoord
    */
    virtual void CalcElemCharges(CFSVector & charges,
				 const StdVector<Elem*> & surfElems,
				 const Vector<Double> & lCoord,
				 const CFSVector & eNormalFluxDensity);
  protected:
  
  };




} // end of namespace

#endif
