#ifndef FILE_MAGFORCEOP_2003
#define FILE_MAGFORCEOP_2003

#include "Forms/baseoperator.hh"
#include "Forms/curlfieldop.hh"


namespace CoupledField
{

  // Forward declaration of classes
  class Grid;
  class Elem;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;

  //! Operator for calculating the magnetic Lorentz force 
  //! J x B
  class MagLorentzForceOp : public BaseOperator
  {

  public:

    //! Constructor

    //! \param ptGrid     (input) Pointer to grid
    //! \param magPotential (input) Pointer to vector containing the magnetic
    //!                           potential for all nodes of domain
    //! \param level      (input) Multigrid level
    MagLorentzForceOp(Grid * ptGrid, 
		BasePDE * ptPDE,
		NodeEQN * ptEQN,
		NodeStoreSol<Double> & magPotential,
		const Integer level, Boolean isaxi);

    //! Destructor
    virtual ~MagLorentzForceOp();
  
    //! Calculate element Lorentz force
    //! \param F              (output) force
    //! \param Jeddy          (input) eddy current density at finite element nodes
    //! \param ptElem         (input)  Pointer to element
    virtual void CalcElemMagLorentzForce(Matrix<Double>& F,
					 Vector<Double>& Jeddy,
					 const Elem * ptElement);


  protected:
  
    //! I'm a class attribute (please add documentation for me)
    CurlNodeOp * curlFieldOp_;

  };

} // end of namespace

#endif
