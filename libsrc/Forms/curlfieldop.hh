#ifndef FILE_CURLFIELDOP_2004
#define FILE_CURLFIELDOP_2004

#include "Forms/baseoperator.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

#ifdef USE_OLAS
#include "olas.hh"
#else
#include "multigrid.hh"
#endif

namespace CoupledField
{

  class Grid;
  class Elem;
  template<class TYPE> class ElemStoreSol;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;
 
  //! Curl operator for edge elements
  
class CurlEdgeOp : public BaseOperator
{  
public:

  //! Constructor
  /*!
    \param ptGrid (input) Pointer to grid
    \param EPotential (input) Pointer to vector containing the electric potential for all nodes of domain
    \param level (input) Multigrid level
    \param algsys (input) pointer to algebraic system
  */
  CurlEdgeOp(Grid * ptGrid,
	     BasePDE * ptPDE,
	     NodeEQN * ptEQN,
	     NodeStoreSol<Double> & sol,
	     const Integer level,
	     BaseSystem * algsys);

  //! Destructor
  virtual ~CurlEdgeOp();
  
  //! Calculate element electric field
  /*!
    \param E (output) Element vector of electric field components
    \param ptElem (input) Pointer to element
    \param LCoord (input) Local coordinates of evaluation point
   */
  virtual void CalcElemCurlEdge(Vector<Double> & E,
				const Elem * ptElement,
				const Vector<Double> & LCoord);
  

  void CalcElemMagVec(Vector<Double> & magVecPot, 
		      const Elem * ptElement,
		      const Vector<Double> & LCoord);
  
protected:
  
  NodeStoreSol<Double> * sol_;
  BaseSystem * algsys_;
};




  //! curl operator for nodeal elements


class CurlNodeOp : public BaseOperator
{  
public:

  //! Constructor
  CurlNodeOp(Grid * ptGrid,
	     BasePDE * ptPDE,
	     NodeEQN * ptEQN,
	     NodeStoreSol<Double> & sol,
	     const Integer level);

  //! Destructor
  virtual ~CurlNodeOp();
  
  void Set2DType(Boolean axi) { isaxi_ = axi;};

  //! Calculate element magnetic field
  /*!
    \param E (output) Element vector of electric field components
    \param ptElem (input) Pointer to element
    \param LCoord (input) Local coordinates of evaluation point
   */
  virtual void CalcElemCurlNode(Vector<Double> & E,
				const Elem * ptElement,
				const Vector<Double> & LCoord);

  void CalcElemMagVec(Vector<Double> & magVecPot, 
		      const Elem * ptElement,
		      const Vector<Double> & LCoord)
  {Error("CalcElemMagVec not implemented",__FILE__,__LINE__);};
  
protected:
  
  NodeStoreSol<Double> * sol_;
};

} // end of namespace

#endif
