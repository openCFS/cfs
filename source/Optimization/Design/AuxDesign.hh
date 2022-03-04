#ifndef AUX_DESIGN_HH_
#define AUX_DESIGN_HH_

#include <stddef.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Utils/StdVector.hh"

namespace CoupledField{
  /** This is an extension to the standard DesignSpace. It adds variables which exceed the
   *  design*elements schema of DesignSpace.
   *  It is required for a slack-variable min alpha for min max problems and is also the
   *  base class for ShapeDesign.
   *  The class handles the case with slack and slack and alpha. The latter is for the bounds 'alpha+slack' and 'alpha-slack' for ban gap maximization */
class Condition;
class Objective;

class AuxDesign : public DesignSpace
{
  public:

    /** @param naux 1 for slack, 2 for slack and alpha, -1 to identify automatically by pn */
    AuxDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method = ErsatzMaterial::NO_METHOD, int naux = -1);

    virtual ~AuxDesign() { } ;

    /** only for slack variable
     * @see DesignSpace::PostInit() */
    virtual void PostInit(int objectives, int constraints);

    /** @see DesignSpace::ReadDesignFromExtern() */
    virtual int ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent = true);

    /** overwrites DesignSpace::CompareDesign() */
    virtual bool CompareDesign(const double* space_in);

    /** writes design to the vector, prepending with shape variables */
    virtual int WriteDesignToExtern(double* space_out, bool scaling = true) const;

    /** write gradient out to the vector, appending with shape gradient
     * Sparse and dense! */
    virtual void WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling = true);

    /** same as in DesignSpace, setting elements to zero, but also aux elements */
    virtual void Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design = DesignElement::DEFAULT);

    virtual void WriteBoundsToExtern(double* x_l, double* x_u) const;

    virtual unsigned int GetNumberOfVariables() const;

    virtual int GetNumberOfAuxParameters() const { return aux_design_.GetSize(); }

    void AddAuxDerivative(Function* f, unsigned int index, double value);

    void AddAuxDerivatives(Objective* objective, Condition* constraint, StdVector<double>& d, double weight);

    /** @see DesignSpace::HasSlackVariable() */
    bool HasSlackVariable() const { return slack_ != NULL; }

    bool HasAlphaVariable() const { assert(!(slack_ == NULL && alpha_ != NULL)); return alpha_ != NULL; }

    BaseDesignElement* GetSlackDesign() { return &(aux_design_[0]); }

    BaseDesignElement* GetAlphaDesign() { return &(aux_design_[1]); }

    /** @see DesignSpace::GetSlackVariable() */
    double GetSlackVariable() const { assert(slack_ != NULL); return aux_design_[0].GetDesign(); }

    double GetAlphaVariable() const { assert(alpha_ != NULL); return aux_design_[1].GetDesign(); }

    /** see DesignSpace::ToInfo() */
    virtual void ToInfo(ErsatzMaterial* em);

    /** @see DesignSpace::GetDesignElement() */
    virtual BaseDesignElement* GetDesignElement(unsigned int idx);

    /** design element with only Aux idx */
    BaseDesignElement* GetAuxDesignElement(unsigned int idx);

  protected:

    /** write the aux gradient part */
    void WriteAuxGradientToExtern(StdVector<double>& out, Function* f, bool scale = true) const;

    /** sparse version of WriteAuxGradientToExtern */
    void WriteSparseAuxGradientToExtern(StdVector<double>& out, Function* f, bool scale = true) const;

    /** compute the offset of aux design. It depends on export_fe_design_, DesignSpace::data and */
    unsigned int AuxDesignOffset() const;

    /** is DesignSpace::data seen from an external optimizer?
     * This covers WriteDesignToExtern(), CompareDesign(), ...
     * If not we either do shape optimization or we do shape mapping where we actually use DesignSpace::data but do not export it */
    bool export_fe_design_;

    /** are the aux parameters located last? Only for non shape opt but with standard simp and shape mapping. false does not mean first
     * but can be intermediate.
      @see AuxDesignOffset() */
    bool tailing_aux_design_;

    StdVector<BaseDesignElement> aux_design_;

    /** factor for scaling the aux_parameters, this scales the gradient compared to material parameters gradient */
    double scaling_;

    /** contains the slack design if slack design */
    PtrParamNode slack_;

    /** contains the alpha design which is only allowed in combination with slack!. */
    PtrParamNode alpha_;

  };

}

#endif /*AUX_DESIGN_HH_*/
