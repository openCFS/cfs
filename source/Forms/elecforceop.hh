// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ELECFORCEOP_2003
#define FILE_ELECFORCEOP_2003

#include <map>

#include "Forms/baseForceOp.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"

namespace CoupledField {
class BaseMaterial;
class EntityIterator;
class EqnMap;
class Grid;
class StdPDE;
struct ResultInfo;
template <class TYPE> class GradientFieldOp;
template <typename T> class NodeStoreSol;
}  // namespace CoupledField


namespace CoupledField
{


  //! Operator for calculating the electric force 
  //! This operator class calculates the electric force in an element
  //! \f$ F_{p,r} = \epsilon_{0} E \cdot J^{-1} \frac{\delta J}{\delta r}
  //! E \left| J \right| - \frac{\epsilon_{0} E \cdot E}{2} \cdot \frac{\delta
  //! \left| J \right|}{\delta r} \f$
  class ElecForceOp : public BaseForceOp
  {

  public:

    //! This is a static const Double

    //! WARN: This violates the ISO C++ standard. Only integral types
    //!          can be static and const!
    //! \todo eps0 violates the ISO C++ standard. Only integral types
    //!       can be static and const!
    //static const Double  eps0 = 8.854187817e-12;
    
    //! Constructor

    //! \param ptGrid     (input) Pointer to grid
    //! \param EPotential (input) Pointer to vector containing the electric
    //!                           potential for all nodes of domain
    ElecForceOp(Grid * ptGrid, 
                StdPDE * ptPDE,
                shared_ptr<EqnMap> eqnMap,
                NodeStoreSol<Double> & EPotential,
                UInt dim,
                std::map<RegionIdType,BaseMaterial*>& matData,
                shared_ptr<ResultInfo> result,
                bool isaxi,
                bool coordUpdate = false );

    //! Destructor
    virtual ~ElecForceOp();


  protected:
  
    //! returns the scalar material value, used for force computation
    virtual Double GetMatVal(RegionIdType actRegion);

    //! computes the field quantity
    virtual void ComputeField(Vector<Double> & Field, 
                              const EntityIterator& ent, 
                              const Vector<Double> & lCoord);

    //! I'm a class attribute (please add documentation for me)
    GradientFieldOp<Double> * gradFieldOp_;


  };

} // end of namespace

#endif
