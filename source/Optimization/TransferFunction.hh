#ifndef TRANSFERFUNCTION_HH_
#define TRANSFERFUNCTION_HH_

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Optimization/Design/BaseDesignElement.hh"
#include "Optimization/Context.hh"
#include "Optimization/Tune.hh"

namespace CoupledField
{
  class SinglePDE;
  class DesignElement;
  class MathParser;

  /** This defines a transfer function, where the standard SIMP is a variant */
  class TransferFunction
  {
    public:
    /** The transfer functions for a design-variable x:
     * <ul>
     * <li>SIMP_TPYE: tf(x) = x^param, tf'(x)=param*x^(param-1)</li>
     * <li>IDENTITY: tf(x) = x, tf'(x)=1</li>
     * <li>RAMP: ... very slow :(</li>
     * <li>FIXED: tf(x) = param == 0 ? 1.0 : param, tf'(x)=0</li>
     * <li>FULL:  tf(x) = 1, tf'(x)=0</li>
     * <li>HEAVISIDE:  tf(x) = (1-exp(-beta * x))^param</li>
     * <li>TANH:  tf(x) =  1 - 1/(exp(2*beta*(x-param)) + 1) scaled for x in [0:1] and y in [0:1]</li>
     * <li>HASHIN_SHTRIKMAN: tf(x) = x/(3-2*x) as in (7) in Bendsoe, Sigmund; Material interpolation schemes
     *                        Upper Hashin-Shtrikman bound for 2D with Poisson's ratio 0.3</li>
     * <li>EXPRESSION: muParser expression in expression child element with rho as input, param as p and beta as b</li>
     * </ul> */
    typedef enum { NO_TYPE = -1, SIMP_TYPE, IDENTITY, RAMP, FIXED, FULL, HEAVISIDE, TANH, HASHIN_SHTRIKMAN, EXPRESSION } Type;

      /** dummy function for TransferFunction */
      TransferFunction();
    
      /** @param pn The pointer to a transfer function
       * @param default_type if no design is given in the xml file use the default_type but checks for NO_TYPE */
      TransferFunction(PtrParamNode pn, BaseDesignElement::Type default_type = BaseDesignElement::NO_TYPE);
     
      /** E.G. for the stresses we temporarily construct a own transfer function */
      TransferFunction(App::Type app, TransferFunction::Type tf_type, double param, BaseDesignElement::Type design = BaseDesignElement::NO_TYPE);

      /** explicit copy constructor to handle math parser properly */
      TransferFunction(const TransferFunction& tf);

      virtual ~TransferFunction();

      /** applies the transformation
       * @param de contains the design value.
       * @param access if SMART and the filter is accordingly defined the filtered design is the base for penalization*/
      double Transform(const DesignElement* de, BaseDesignElement::Access access, bool lower_bimat = false) const;

      /** applies the transformation
       * @param value is the design value. de may be NULL. It is only used in logging and if type is FULL! */
      double Transform(double value, bool lower_bimat = false, const BaseDesignElement* de = NULL) const;

      /** applies the first derivative of the transformation
       * @see Transform() */
      double Derivative(const DesignElement* de, BaseDesignElement::Access access, bool lower_bimat = false) const;

      /** applies the first derivative of the transformation
       * @see Transform() */
      double Derivative(double value, bool lower_bimat = false) const;

      /** Gives the standard, non-mass App::Type to find the transfer-function. Note that the application for transfer functions
       * does not coincide with pde application due to the laplace stuff */
      static App::Type Default(const Context* ctxt);

      /** see the other Default */
      static App::Type Default(BaseDesignElement::Type type, const Context* = NULL);

      App::Type GetApplication() { return application_; }
      
      BaseDesignElement::Type GetDesign() { return design_; }
      
      Type GetType() const { return type_; }
      
      double GetParam() const { return param_; }

      /** sets the disable stuff */
      void Enable(bool enable);

      bool IsPenalized() const;

      /** Optional for Heaviside and tanh if physical lower bound are desired */
      void SetScaling(double val) { scaling_ = val; }

      void SetOffset(double val) { offset_ = val; }

      /** check if tune is set, and if so, register */
      void RegisterTune(Optimization* opt);

      /** gives debug information */
      std::string ToString();

      /** Writes the attributes, not the base element
       * @param skip_relation no design and no application */
      void ToInfo(PtrParamNode in, bool skip_relation=false) const;
      
      static Enum<Type> type;

      TransferFunction& operator=(const TransferFunction& other);

    private:
      /** for constructor and copy constructor */
      void InitParser(const std::string& func_expr, const std::string& deriv_expr);

      /** common for copy constructor and assignment operator */
      void Copy(const TransferFunction& other);

      static void SetEnums();

      BaseDesignElement::Type design_ = BaseDesignElement::DEFAULT;

      /** our own type of transfer function (SIMP, IDENTIY). Might be IDENTY if disabled. */
      Type type_ = IDENTITY;
      
      /** our real type of transfer function, only set if disabled, else NO_TYPE. */
      Type orgType_ = NO_TYPE;
       
      /** type of application */
      App::Type application_ = App::NO_APP;
       
      /** e.g. the exponent for SIMP, not used in IDENTIY. Can be tuned! */
      double param_ = 0.0;

      /** heaviside and tanh have also beta */
      double beta_ = -1.0;

      /** Heaviside and tanh have this scaling if the design is physical. Set in DesignSpace::DetermineLowerBound.
       * @see offset_ */
      double scaling_ = 1.0;

      /** @see scaling_ */
      double offset_ = 0.0;

      /** Init() if we tune the panalty param and Register() that objection which is actually used (copy constructor hell!) */
      Tune tune_;

      /** for expression */
      MathParser* parser_ = nullptr;

      /** math parser handle */
      unsigned int function_handle_ = 0;
      unsigned int derivative_handle_ = 0;

      /** the expression uses this value via reference.
       * We use this value "temporarily", hence setting it in Transfer() and Derivative() allows these function still to be const
       * We need Copy() to have this reference for each instance! */
      mutable double expression_rho_ = 0.0;
  };



} // end of namespace

#endif /*TRANSFERFUNCTION_HH_*/
