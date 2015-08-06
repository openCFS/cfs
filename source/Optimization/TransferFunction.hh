#ifndef TRANSFERFUNCTION_HH_
#define TRANSFERFUNCTION_HH_

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Optimization.hh"

namespace CoupledField
{
  class SinglePDE;

  /** This defines a transfer function, where the standard SIMP is a variant */
  class TransferFunction
  {
    public:
    /** The transfer functions for a design-variable x:
     * <ul>
     * <li>SIMP_TPYE: tf(x) = x^param, tf'(x)=param*x^(param-1)</li>
     * <li>IDENTIT: tf(x) = x, tf'(x)=1</li>
     * <li>RAMP: ... very slow :(</li>
     * <li>FIXED: tf(x) = param == 0 ? 1.0 : param, tf'(x)=0</li>
     * <li>FULL:  tf(x) = 1, tf'(x)=0</li>
     * <li>HEAVISIDE:  tf(x) = (1-exp(-beta * x))^param</li>
     * <li>TANH:  tf(x) =  1 - 1/(exp(2*beta*(x-param)) + 1) scaled for x in [0:1] and y in [0:1]</li>
     * </ul> */
    typedef enum { NO_TYPE = -1, SIMP_TYPE, IDENTITY, RAMP, FIXED, FULL, HEAVISIDE, TANH } Type;

      /** dummy function for StdVector */
      TransferFunction();
    
      /** @param pn The pointer to a transfer function
       * @param default_type if no design is given in the xml file use the default_type but checks for NO_TYPE */
      TransferFunction(PtrParamNode pn, DesignElement::Type default_type = DesignElement::NO_TYPE);
     
      /** E.G. for the stresses we temporarily construct a own transfer function */
      TransferFunction(Optimization::Application app, TransferFunction::Type tf_type, double param, DesignElement::Type design = DesignElement::NO_TYPE);
    
      /** applies the transformation
       * @param de contains the design value.
       * @param access if SMART and the filter is accordingly defined the filtered design is the base for penalization*/
      double Transform(const DesignElement* de, DesignElement::Access access, bool lower_bimat = false) const;

      /** applies the transformation
       * @param value is the design value. de may be NULL. It is only used in logging and if type is FULL! */
      double Transform(double value, bool lower_bimat = false, const DesignElement* de = NULL) const;

      /** applies the first derivative of the transformation
       * @see Transform() */
      double Derivative(const DesignElement* de, DesignElement::Access access, bool lower_bimat = false) const;

      /** applies the first derivative of the transformation
       * @see Transform() */
      double Derivative(double value, bool lower_bimat = false) const;

      /** Gives the standard, non-mass Application to find the transfer-function. Note that the application for transfer functions
       * does not coincide with pde application due to the laplace stuff */
      static Optimization::Application Default(const SinglePDE* pde);

      /** see the other Default */
      static Optimization::Application Default(DesignElement::Type type, const SinglePDE* = NULL);

      Optimization::Application GetApplication() { return application_; }
      
      DesignElement::Type GetDesign() { return design_; }
      
      Type GetType() const { return type_; }
      
      double GetParam() const { return param_; }

      /** sets the disable stuff */
      void Enable(bool enable);

      bool IsPenalized() const;

      /** Optional for Heaviside and tanh if physical lower bound are desired */
      void SetScaling(double val) { scaling_ = val; }

      void SetOffset(double val) { offset_ = val; }

      /** gives debug information */
      std::string ToString();

      /** Writes the attributes, not the base element */
      void ToInfo(PtrParamNode in) const;
      
      static Enum<Type> type;

      DesignElement::Type design_;


    private:
      /** our design type (DENSITY, POLARIZATION)  */
      //DesignElement::Type design_;
       
      /** our own type of transfer function (SIMP, IDENTIY). Might be IDENTY if disabled. */
      Type type_;
      
      /** our real type of transfer function, only set if disabled, else NO_TYPE. */
      Type orgType_; 
       
      /** on of the ELAST, MECH, PIEZO_COUPLING */ 
      Optimization::Application application_;
       
      /** the exponent for SIMP, not used in IDENTIY */
      double param_;

      /** heaviside and tanh have also beta */
      double beta_;

      /** Heaviside and tanh have this scaling if the design is physical. Set in DesignSpace::DetermineLowerBound.
       * @see offset_ */
      double scaling_;

      /** @see scaling_ */
      double offset_;

      static void SetEnums();
  };
} // end of namespace

#endif /*TRANSFERFUNCTION_HH_*/
