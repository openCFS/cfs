#ifndef FILE_ELECFIELDOP_2003
#define FILE_ELECFIELDOP_2003

#include <Forms/baseoperator.hh>
#include <Utils/vector.hh>
#include <Utils/array.hh>
#include <Matrix/matrix.hh>
#include <multigrid.hh>

namespace CoupledField
{

  class Grid;
  class Elem;
 
  //! Operator for calculating the electric field
  /*! This operator class calculates the electric field in an element
      at an arbitrary point
   */


class ElecFieldOp : public BaseOperator
{

  
public:

  //! Constructor
  /*!
    \param ptGrid (input) Pointer to grid
    \param EPotential (input) Pointer to vector containing the electric potential for all nodes of domain
    \param level (input) Multigrid level
  */
  ElecFieldOp(Grid * ptGrid,
	      BasePDE * ptPDE,
	      std::vector<Integer> * ptMesh2PDENode,
	      Vector<Double> * EPotential,
	      const Integer level);

  //! Destructor
  virtual ~ElecFieldOp();
  
  //! Calculate element electric field
  /*!
    \param E (output) Element vector of electric field components
    \param ptElem (input) Pointer to element
    \param LCoord (input) Local coordinates of evaluation point
   */
  virtual void CalcElemElecField(Vector<Double> & E,
				 const Elem * ptElement,
				 const std::vector<Double> & LCoord);
  

  //! Calculate electric field for list of subdomains
  /*!
    \param E (input) Array of vectors containing Electric field components for whole domain
    \f[ E[0] = \left( E_{1,x}, E_{2,x}, \cdots \right),
        E[1] = \left( E_{1,y}, E_{2,y}, \cdots \right), \cdots 
    \f]
    \param SD (input) Name of the subdomain
    \param LCoord (input) Local coordinates of evalutation point
  */
  virtual void CalcSDElecField(Array<Double> & E,
			       const std::vector<std::string> & SD,
			       const std::vector<Double> & LCoord);
			    			       

protected:
  
  Vector<Double> * EPotential_;

};





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
	     std::vector<Integer> * ptMesh2PDENode,
	     Vector<Double> * sol,
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
				const std::vector<Double> & LCoord);
  

  void CalcElemMagVec(Vector<Double> & magVecPot, 
		      const Elem * ptElement,
		      const std::vector<Double> & LCoord);
  
protected:
  
  Vector<Double> * sol_;
  BaseSystem * algsys_;
};


class CurlNodeOp : public BaseOperator
{  
public:

  //! Constructor
  CurlNodeOp(Grid * ptGrid,
	     BasePDE * ptPDE,
	     std::vector<Integer> * ptMesh2PDENode,
	     Vector<Double> * sol,
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
				const std::vector<Double> & LCoord);

  void CalcElemMagVec(Vector<Double> & magVecPot, 
		      const Elem * ptElement,
		      const std::vector<Double> & LCoord)
  {Error("CalcElemMagVec not implemented",__FILE__,__LINE__);};
  
protected:
  
  Vector<Double> * sol_;
  Boolean isaxi_;
};

} // end of namespace

#endif
