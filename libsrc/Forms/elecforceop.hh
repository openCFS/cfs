#ifndef FILE_ELECFORCEOP_2003
#define FILE_ELECFORCEOP_2003

#include "Forms/baseForceOp.hh"
#include "Forms/gradfieldop.hh"


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

    //! Warning: This violates the ISO C++ standard. Only integral types
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
                NodeEQN * ptEQN,
                NodeStoreSol<Double> & EPotential,
                UInt dim,
                std::map<RegionIdType,BaseMaterial*>& matData,
                Boolean isaxi);

    //! Destructor
    virtual ~ElecForceOp();


  protected:
  
    //! returns the scalar material value, used for force computation
    virtual Double GetMatVal(RegionIdType actRegion);

    //! computes the field quantity
    virtual void ComputeField(Vector<Double> & Field, const Elem * ptElement,
                              const Vector<Double> & lCoord);

    //! I'm a class attribute (please add documentation for me)
    GradientFieldOp<Double> * gradFieldOp_;


  };

} // end of namespace

#endif
