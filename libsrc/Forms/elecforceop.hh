#ifndef FILE_ELECFORCEOP_2003
#define FILE_ELECFORCEOP_2003

#include <Forms/baseoperator.hh>
#include <Forms/elecfieldop.hh>


namespace CoupledField
{

  class Grid;
  class Elem;
 
  //! Operator for calculating the electric force 
  /*! This operator class calculates the electric force in an element
      \f$ F_{p,r} =   \epsilon_{0} E \cdot J^{-1} \frac{\delta J}{\delta r} E \left| J \right|  -
      \frac{\epsilon_{0} E \cdot E}{2} \cdot \frac{\delta \left| J \right|}{\delta r} \f$
   */


class ElecForceOp : public BaseOperator
{

  
public:
  
  static const Double  eps0 = 8.854187817e-12;
  //! Constructor
  /*!
    \param ptGrid (input) Pointer to grid
    \param EPotential (input) Pointer to vector containing the electric potential for all nodes of domain
    \param level (input) Multigrid level
  */
  ElecForceOp(Grid * ptGrid, 
	      std::vector<Integer> * ptMesh2PDENode,
	      Vector<Double> * EPotential,
	      const Integer level);
  
  //! Destructor
  virtual ~ElecForceOp();
  
  //! Calculate element electric force pressure
  /*!
    \param F (output) Element vector of electric field pressure components
    \param ptElem (input) Pointer to element
    \param IsBoundaryNode (input) contains 1, if corresponding node is a boundary node,
    otherwise 0
    \param LCoord (input) Local coordinates of evaluation point
  */
  virtual void CalcElemElecForce(Vector<Double> & F,
				     const Elem * ptElement,
				     Double epsilon,
				     const std::vector<ShortInt> & IsBoundaryNode);
  

			    			       

protected:
  
  ElecFieldOp * ElecFieldOp_;

};

} // end of namespace

#endif
