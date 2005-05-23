#ifndef FILE_BASEFORCEOP_2005
#define FILE_BASEFORCEOP_2005

#include "Forms/baseoperator.hh"


namespace CoupledField
{

  // Forward declaration of classes
  class Grid;
  class Elem;
  template<class TYPE> class Vector;
  template<class TYPE> class Matrix;

  //! Basae operator for calculating the force with VWP 
  class BaseForceOp : public BaseOperator
  {

  public:

    //! Constructor

    //! \param ptGrid     (input) Pointer to grid
    //! \param sol (input) Pointer to vector containing the electric
    //!                           potential for all nodes of domain
    BaseForceOp(Grid * ptGrid, 
		StdPDE * ptPDE,
		NodeEQN * ptEQN,
		NodeStoreSol<Double> & sol,
		Integer dim,
		MaterialData* &matData,
		StdVector<RegionIdType> & allSubdoms,
		Boolean isaxi);

    //! Destructor
    virtual ~BaseForceOp();

    //!
    void Setup(StdVector<RegionIdType> & neighRegions, 
	       StdVector<Integer>& couplingnodes);

    //!  
    void CalcNodeForce(Vector<Double> & force, Vector<Double> & totalForce );

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

    //! returns the scalar material value, used for force computation
    virtual Double GetMatVal(Integer actSD)=0;

    //! computes the field quantity
    virtual void ComputeField(Vector<Double> & Field, const Elem * ptElement,
			      const Vector<Double> & lCoord)=0;

    //! coupling nodes
    StdVector<Integer> couplingNodes_;

    //!neighbor-regions
    StdVector<RegionIdType> neighRegions_;

    //! interface elements
    StdVector<Elem*> interfaceElems_;

    //! vector containing flag array for element boundary nodes
    StdVector<StdVector<ShortInt> > isBoundaryNode_;

    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<Integer> > elemNodeToCouplingNode_; 

    //! dimension
    Integer dim_;

    //! material data
    MaterialData* materialData_;

    //! all subdomains belonging to the PDE
    StdVector<RegionIdType> PDEsubdoms_;

    //! solution type;
    SolutionType solType_;

    // sign of force
    Double sign_;

  };

} // end of namespace

#endif
