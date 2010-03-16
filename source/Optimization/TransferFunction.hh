#ifndef TRANSFERFUNCTION_HH_
#define TRANSFERFUNCTION_HH_

#include <string>
#include "Optimization/Optimization.hh"

namespace CoupledField
{

  /** This defines a transfer function, where the standard SIMP is a variant */
  class TransferFunction
  {
    public:
      /** dummy function for StdVector */
      TransferFunction();
    
      /** @param pn The pointer to a transfer function
       * @param default_type if no design is given in the xml file use the default_type but checks for NO_TYPE */
      TransferFunction(PtrParamNode pn, DesignElement::Type default_type = DesignElement::NO_TYPE);
     
      /** The transfer functions for a design-variable x:
       * <ul>
       * <li>SIMP_TPYE: tf(x) = x^param, tf'(x)=param*x^(param-1)</li>
       * <li>IDENTIT: tf(x) = x, tf'(x)=1</li>
       * <li>RAMP: ... very slow :(</li>
       * <li>FIXED: tf(x) = param == 0 ? 1.0 : param, tf'(x)=0</li>
       * <li>FULL:  td(x) = 1, tf'(x)=0</li>
       * </ul> */
      typedef enum { NO_TYPE = -1, SIMP_TYPE, IDENTITY, RAMP, FIXED, FULL } Type;
    
      /** applies the transformation */
      double Transform(DesignElement* de); 

      /** applies the first derivative of the transformation */
      double Derivative(DesignElement* de);
     
      Optimization::Application GetApplication() { return application_; }
      
      DesignElement::Type GetDesign() { return design_; }
      
      Type GetType() { return type_; }

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
