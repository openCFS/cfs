#ifndef FILE_MAGFORCEOP_2003
#define FILE_MAGFORCEOP_2003

#include "Forms/baseoperator.hh"
#include "Forms/curlfieldop.hh"
#include "Forms/baseForceOp.hh"


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
    MagLorentzForceOp(Grid * ptGrid, 
		      StdPDE * ptPDE,
		      NodeEQN * ptEQN,
		      NodeStoreSol<Double> & magPotential,
		      Boolean isaxi);

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



  //! Operator for calculating the magnetic force 
  //! This operator class calculates the electric force in an element
  //! \f$ F_{p,r} = \frac{1}{\mu_0} B \cdot J^{-1} \frac{\delta J}{\delta r}
  //! B \left| J \right| - \frac{ B \cdot B}{2 \mu_0} \cdot \frac{\delta
  //! \left| J \right|}{\delta r} \f$
  class MagForceOp : public BaseForceOp
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
    //! \param sol        (input) Pointer to vector containing the magnetic
    //!                           vector potential for all nodes of domain
    MagForceOp(Grid * ptGrid, 
	       StdPDE * ptPDE,
	       NodeEQN * ptEQN,
	       NodeStoreSol<Double> & sol,
	       Integer dim,
	       MaterialData* &matData,
	       StdVector<RegionIdType>& allSubdoms,
	       Boolean isaxi);

    //! Destructor
    virtual ~MagForceOp();


  protected:
  
    //! returns the scalar material value, used for force computation
    virtual Double GetMatVal(Integer actSD);

    //! computes the field quantity
    virtual void ComputeField(Vector<Double> & Field, const Elem * ptElement,
			      const Vector<Double> & lCoord);

    //! I'm a class attribute (please add documentation for me)
    CurlNodeOp * curlFieldOp_;


  };

} // end of namespace

#endif
