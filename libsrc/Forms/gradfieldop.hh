#ifndef FILE_GRADIENTFIELDOP_2004
#define FILE_GRADIENTFIELDOP_2004

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
  template<class TYPE> class NodeStoreSol;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;
  
  //! Operator for calculating gradient fields
  
  //! This operator class calculates elemt gradient fields like electric
  //! field or electric flux density at arbitrary points for a given
  //! potential field

class GradientFieldOp : public BaseOperator
{

  
public:

  //! Constructor
  /*!
    \param ptGrid (input) Pointer to grid
    \param ptPDE (input) Pointer to PDE
    \param ptEQN (input) Pointer to EQN
    \param potential (input) NodeStoreSol containing nodal potential
    \param solType (input) SolutionType of the potentialField
    \param level (input) Multigrid level
    \param isaxi (input) Flag for axi-symmetric geomtetry
  */
  GradientFieldOp(Grid * ptGrid,
		  BasePDE * ptPDE,
		  NodeEQN * ptEQN,
		  NodeStoreSol<Double> & potential,
		  const SolutionType solType,
		  const Integer level,
		  Boolean isaxi=FALSE);
  
  //! Destructor
  virtual ~GradientFieldOp();
  
  //! Calculate element gradient field
  /*!
    \param elemField (output) Element vector of gradient field
    \param ptElem (input) Pointer to element
    \param LCoord (input) Local coordinates of evaluation point
    \param factor (input) Scaling factor (e.g. permittivity for E-Field)
   */
  virtual void CalcElemGradField(Vector<Double> & elemField,
				 const Elem * ptElement,
				 const Vector<Double> & LCoord,
				 const Double factor);
  

  //! Calculate electric field for list of subdomains
  /*!
    \param elemField (input) Vector containing
    \f[ elemFIeld[0] = \left( elemField_{1,x}, elemField_{2,x}, \cdots \right),
        elemField[1] = \left( elemField_{1,y}, elemField_{2,y}, \cdots \right), \cdots 
    \f]
    \param SD (input) Name of the subdomain
    \param LCoord (input) Local coordinates of evalutation point
  */
  virtual void CalcSDGradField(Vector<Double> & elemField,
			       const StdVector<std::string> & SD,
			       const Vector<Double> & lCoord,
			       const Vector<Double> & factors);
			    			       

protected:
  
  //! StoreSolution containing potentialfield
  NodeStoreSol<Double> * potential_;

  //! Soltution type of the potential field
  SolutionType solType_;
};

} // end of namespace

#endif
