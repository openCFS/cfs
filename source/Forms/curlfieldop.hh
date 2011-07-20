// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CURLFIELDOP_2004
#define FILE_CURLFIELDOP_2004

#include "Forms/baseoperator.hh"
#include "Utils/tools.hh"
#include "Utils/StdVector.hh"

#include "OLAS/algsys/basesystem.hh"

namespace CoupledField
{

  class Grid;
  struct Elem;
  template<class TYPE> class ElemStoreSol;
  template<class TYPE> class NodeStoreSol;
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
      \param algsys (input) pointer to algebraic system
    */
    CurlEdgeOp(Grid * ptGrid,
               StdPDE * ptPDE,
               shared_ptr<EqnMap> eqnMap,
               NodeStoreSol<Double> & sol,
               BaseSystem * algsys,
               bool coordUpdate = false );

    //! Destructor
    virtual ~CurlEdgeOp();
  
    //! Calculate element electric field
    /*!
      \param E (output) Element vector of electric field components
      \param ptElem (input) Pointer to element
      \param LCoord (input) Local coordinates of evaluation point
    */
    virtual void CalcElemCurlEdge(Vector<Double> & E,
                                  const EntityIterator& ent,
                                  const Vector<Double> & lCoord);
  

    void CalcElemMagVec(Vector<Double> & magVecPot, 
                        const EntityIterator& ent,
                        const Vector<Double> & LCoord);
  
  protected:
  
    NodeStoreSol<Double> * sol_;
    BaseSystem * algsys_;
  };




  //! curl operator for nodeal elements

  template<class TYPE>
  class CurlNodeOp : public BaseOperator
  {  
  public:

    //! Constructor
    CurlNodeOp(Grid * ptGrid,
               StdPDE * ptPDE,
               shared_ptr<EqnMap> eqnMap,
               NodeStoreSol<TYPE> & sol,
               bool coordUpdate = false );

    //! Destructor
    virtual ~CurlNodeOp();
  
    void Set2DType(bool axi) { isaxi_ = axi;};

    //! Calculate element magnetic field
    /*!
      \param E (output) Element vector of electric field components
      \param ptElem (input) Pointer to element
      \param LCoord (input) Local coordinates of evaluation point
    */
    virtual void CalcElemCurlNode(Vector<TYPE> & E,
                                  const EntityIterator& ent,
                                  const Vector<Double> & lCoord );

    void CalcElemMagVec(Vector<Double> & magVecPot, 
                        const Elem * ptElement,
                        const Vector<Double> & LCoord)
    {EXCEPTION("CalcElemMagVec not implemented");};
  
  protected:
  
    NodeStoreSol<TYPE> * sol_;
  };

} // end of namespace

#endif
