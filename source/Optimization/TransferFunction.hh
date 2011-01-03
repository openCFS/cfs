#ifndef TRANSFERFUNCTION_HH_
#define TRANSFERFUNCTION_HH_

#include <string>
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
     * <li>FULL:  td(x) = 1, tf'(x)=0</li>
     * </ul> */
    typedef enum { NO_TYPE = -1, SIMP_TYPE, IDENTITY, RAMP, FIXED, FULL } Type;

      /** dummy function for StdVector */
      TransferFunction();
    
      /** @param pn The pointer to a transfer function
       * @param default_type if no design is given in the xml file use the default_type but checks for NO_TYPE */
      TransferFunction(PtrParamNode pn, DesignElement::Type default_type = DesignElement::NO_TYPE);
     
      /** E.G. for the stresses we temporarily construct a own transfer function */
      TransferFunction(Optimization::Application app, TransferFunction::Type tf_type, double param, DesignElement::Type design = DesignElement::NO_TYPE);
    
      /** applies the transformation
       * @param de containts the design value
       * @param access if SMART and the filter is accordingly defined the filtered design is the base for penalization*/
      double Transform(const DesignElement* de, DesignElement::Access access, double external_value = -13.456) const;

      /** applies the first derivative of the transformation
       * @see Transform() */
      double Derivative(const DesignElement* de, DesignElement::Access access) const;

      /** Gives the standard, non-mass Application to find the transfer-function. Note that the application for transfer functions
       * does not coincide with pde application due to the laplace stuff */
      static Optimization::Application Default(const SinglePDE* pde);

      /** see the other Default */
      static Optimization::Application Default(DesignElement::Type type);

      Optimization::Application GetApplication() { return application_; }
      
      DesignElement::Type GetDesign() { return design_; }
      
      Type GetType() { return type_; }
      
      double GetParam() { return param_; }

      /** @return true for SIMP with p != 1 and RAMP != 0 */
      bool IsPenalized() const;

      /** sets the disable stuff */
      void Enable(bool enable) 
      { 
        // try to handle to much toggling cases
        if(enable) 
        {
          type_ = orgType_;
          assert(type_ != NO_TYPE);
        } 
        else
        {  
           // only disable if we are enabled
           assert(type_ != NO_TYPE);
          orgType_ = type_;
          type_ = NO_TYPE;
        }
      }
     
      /** gives debug information */
      std::string ToString();

      /** Writes the attributes, not the base element */
      void ToInfo(PtrParamNode in) const;
      
      static Enum<Type> type;

    private:
      /** our design type (DENSITY, POLARIZATION)  */
      DesignElement::Type design_;
       
      /** our own type of transfer function (SIMP, IDENTIY). Might be IDENTY if disabled. */
      Type type_;
      
      /** our real type of transder function, only set if disabled, else NO_TYPE. */
      Type orgType_; 
       
      /** on of the ELAST, MECH, PIEZO_COUPLING */ 
      Optimization::Application application_;
       
      /** the exponent for SIMP, not used in IDENTIY */
      double param_;

      static void SetEnums();
  };
} // end of namespace

#endif /*TRANSFERFUNCTION_HH_*/
