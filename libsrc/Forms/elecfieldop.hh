#ifndef FILE_ELECFIELDOP_2003
#define FILE_ELECFIELDOP_2003

#include <Forms/baseoperator.hh>
#include <Utils/vector.hh>
#include <Matrix/matrix.hh>


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
  virtual void CalcSDElecField(Vector<Double> * E,
			       const std::vector<std::string> & SD,
			       const std::vector<Double> & LCoord);
			    			       

protected:
  
  Vector<Double> * EPotential_;

};

} // end of namespace

#endif
