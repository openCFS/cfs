#ifndef FILE_ELECFORCEOP_2003
#define FILE_ELECFORCEOP_2003

#include "Forms/baseoperator.hh"
#include "Forms/gradfieldop.hh"


namespace CoupledField
{

  // Forward declaration of classes
  class Grid;
  class Elem;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;

  //! Operator for calculating the electric force 
  //! This operator class calculates the electric force in an element
  //! \f$ F_{p,r} = \epsilon_{0} E \cdot J^{-1} \frac{\delta J}{\delta r}
  //! E \left| J \right| - \frac{\epsilon_{0} E \cdot E}{2} \cdot \frac{\delta
  //! \left| J \right|}{\delta r} \f$
  class ElecForceOp : public BaseOperator
  {

  public:

    //! This is a static const Double

    //! Warning: This violates the ISO C++ standard. Only integral types
    //!          can be static and const!
    //! \todo eps0 violates the ISO C++ standard. Only integral types
    //!       can be static and const!
    //static const Double  eps0 = 8.854187817e-12;

    //! Constructor

    //! \param ptGrid     (input) Pointer to grid
    //! \param EPotential (input) Pointer to vector containing the electric
    //!                           potential for all nodes of domain
    //! \param level      (input) Multigrid level
    ElecForceOp(Grid * ptGrid, 
		BasePDE * ptPDE,
		NodeEQN * ptEQN,
		NodeStoreSol<Double> & EPotential,
		const Integer level, Boolean isaxi);

    //! Destructor
    virtual ~ElecForceOp();
  
    //! Calculate element electric force pressure

    //! \param F              (output) Array containing nodal forces
    //!                                (dim x nodes) of each element
    //! \param ptElem         (input)  Pointer to element
    //! \param IsBoundaryNode (input)  contains 1, if corresponding node is a
    //!                                boundary node, otherwise 0
    //! \param LCoord         (input)  Local coordinates of evaluation point
    virtual void CalcElemElecForce(ElemStoreSol<Double> & F,
				   const Elem * ptElement,
				   Double epsilon,
				   const StdVector<ShortInt> & IsBoundaryNode);

  protected:
  
    //! Calculates the expression \f[ \frac{\delta \vert J \vert}{\delta r} /f]

    //! \param J (input) Jacobian matrix
    //! \param J_dr (input) derivative of Jacobian matrix in r-direction
    //! \param dim (input) dimension (= direction) of r
    Double CalcDetJDr(Matrix<Double> &J, Matrix<Double> &Jdr, Integer dim);

    //! I'm a class attribute (please add documentation for me)
    GradientFieldOp * gradFieldOp_;

  };

} // end of namespace

#endif
