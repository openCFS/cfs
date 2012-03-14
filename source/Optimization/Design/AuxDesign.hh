#ifndef AUX_DESIGN_HH_
#define AUX_DESIGN_HH_

#include <stddef.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Utils/StdVector.hh"

namespace CoupledField{
  /** This is an extension to the standard DesignSpace. It adds variables which exceed the
   *  design*elements schema of DesignSpace.
   *  It is required for a slack-variable min alpha for min max problems and is also the
   *  base class for ShapeDesign */
class Condition;
class Objective;

class AuxDesign : public DesignSpace
{
  public:

    AuxDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD);

    virtual ~AuxDesign() { } ;

    /** @see DesignSpace::ReadDesignFromExtern() */
    virtual int ReadDesignFromExtern(const double* space_in);

    /** overwrites DesignSpace::CompareDesign() */
    virtual bool CompareDesign(const double* space_in);

    /** writes design to the vector, prepending with shape variables */
    virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const;

    /** write gradient out to the vector, prepending with shape gradient */
    virtual void WriteDenseGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs,
                               DesignElement::Access access, Condition* constraint = NULL, bool scaling = true) const;

    /** write the aux gradient part */
    void WriteAuxGradientToExtern(StdVector<double>& out, Condition* constraint, bool scale = true) const;

    /** same as in DesignSpace, setting elements to zero, but also aux elements */
    virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);

    virtual void WriteBoundsToExtern(double* x_l, double* x_u) const;

    virtual unsigned int GetNumberOfVariables() const;

    int GetNumberOfAuxParameters() const { return aux_design_.GetSize(); }

    void AddAuxDerivatives(Objective* objective, Condition* constraint, StdVector<double>& d, double weight);



  protected:

    bool alsomatopt_;

    StdVector<BaseDesignElement> aux_design_;

    /** factor for scaling the aux_parameters, this scales the gradient compared to material parameters gradient */
    double scaling_;

  };

}

#endif /*AUX_DESIGN_HH_*/
