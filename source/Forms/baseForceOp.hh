// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEFORCEOP_2005
#define FILE_BASEFORCEOP_2005

#include <map>

#include "Forms/baseoperator.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class EntityIterator;
class EqnMap;
class StdPDE;
template <typename T> class ElemStoreSol;
template <typename T> class NodeStoreSol;
}  // namespace CoupledField

namespace CoupledField
{

  class BaseMaterial;
  // Forward declaration of classes
  class Grid;
  struct Elem;
  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;

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
                shared_ptr<EqnMap> eqnMap,
                NodeStoreSol<Double> & sol,
                UInt dim,
                std::map<RegionIdType, BaseMaterial*>& materials,
                bool isaxi, bool coordUpdate = false);

    //! Destructor
    virtual ~BaseForceOp();

    //!
    void Setup(StdVector<RegionIdType> & neighRegions, 
               StdVector<UInt>& couplingnodes);
    //! returns the number of coupling nodes
    UInt GetNumCouplingNodes() {
      return couplingNodes_.GetSize();
    }

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
    Double CalcDetJDr(Matrix<Double> &J, Matrix<Double> &Jdr, UInt dim);

    //! returns the scalar material value, used for force computation
    virtual Double GetMatVal(RegionIdType regionId)=0;

    //! computes the field quantity
    virtual void ComputeField(Vector<Double> & Field, const EntityIterator& ent,
                              const Vector<Double> & lCoord)=0;

    //! coupling nodes
    StdVector<UInt> couplingNodes_;

    //!neighbor-regions
    StdVector<RegionIdType> neighRegions_;

    //! interface elements
    StdVector<Elem*> interfaceElems_;

    //! vector containing flag array for element boundary nodes
    StdVector<StdVector<ShortInt> > isBoundaryNode_;

    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<UInt> > elemNodeToCouplingNode_; 

    //! dimension
    UInt dim_;

    //! material data of all subdomains
    std::map<RegionIdType, BaseMaterial*> materials_;

    //! solution type;
    SolutionType solType_;

    // sign of force
    Double sign_;

  };

} // end of namespace

#endif
