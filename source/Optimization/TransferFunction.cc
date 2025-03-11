#include <assert.h>
#include <cmath>
#include <stddef.h>
#include <ostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "General/Exception.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Tune.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/mathParser/mathParser.hh"

using namespace CoupledField;

DEFINE_LOG(trans, "transferFunction")

Enum<TransferFunction::Type> TransferFunction::type;

TransferFunction::TransferFunction()
{
  LOG_DBG(trans) << "TF::TF - empty";

  // do nothing - is for StdVector only
  // this constructor is also used for providing a identity transfer function during parametric material optimization
  if(type.map.size() == 0)
    SetEnums();
}


TransferFunction::TransferFunction(App::Type app, TransferFunction::Type tf_type, double param, DesignElement::Type design)
{
  LOG_DBG(trans) << "TF::TF(app,...)";
  assert(type.map.size() > 0);

  this->type_ = tf_type;
  this->application_ = app;
  this->design_ = design;
  this->param_ = param;

  assert(type_ != EXPRESSION);
}


TransferFunction::TransferFunction(PtrParamNode pn, DesignElement::Type default_design)
{
  LOG_DBG(trans) << "TF::TF(pn,...)";
  // initialize the static Enum the first time
  if(type.map.size() == 0)
    SetEnums();
  type_ = type.Parse(pn->Get("type")->As<std::string>());
  application_ = Optimization::application.Parse(pn->Get("application")->As<std::string>());
  if(pn->Has("design"))
    design_ = DesignElement::type.Parse(pn->Get("design")->As<std::string>());
  else
  {
    if(default_design == DesignElement::NO_TYPE) 
      throw Exception("Set the 'design' attribute for 'transferFunction' when using multiple design variables.");
    design_ = default_design;
  }

  if(pn->Has("param"))
  {
    if(pn->Get("param")->As<string>() == "tune")
    {
      tune_.Init(pn->Get("tune", ParamNode::PASS), Tune::PENALTY);
      param_ = tune_.start;
    }
    else
      param_ = pn->Get("param")->As<double>();
  }
  else
    param_ = 1.0;

  // one could also tune beta - so implement to have two BETA for Tune
  beta_ = pn->Has("beta") ? pn->Get("beta")->As<double>() : -1.0;
  
  // validate param
  if(!pn->Has("param") && (type_ == SIMP_TYPE || type_ == RAMP))
    throw Exception("transfer function '" + type.ToString(type_) + "' requires a parameter");
  
  if((type_ == HEAVISIDE || type_ == TANH) && (!pn->Has("param") || !pn->Has("beta")))
    throw Exception("transfer function '" + type.ToString(type_) + "' requires a 'param' and 'beta' to be set");

  if(type_ == TANH && (param_ < 0.0 || param_ > 1.0 )) // ignore exotic design variables!
    throw Exception("'param' for transfer function 'tanh' out of bound [0:1]");

  if(type_ == EXPRESSION)
  {
    if(!pn->Has("expression"))
      throw Exception("'child 'expression' required for corresponding transfer function");
    InitParser(pn->Get("expression/function")->As<string>(), pn->Get("expression/derivative")->As<string>());
  }

  // validate design/application
  if(design_ == DesignElement::POLARIZATION && application_ != App::PIEZO_COUPLING)
    throw Exception("transfer functions for 'polarization' can only be '" 
        + Optimization::application.ToString(App::PIEZO_COUPLING) + "'");

  LOG_DBG(trans) << "TF::TF " << ToString();
}

TransferFunction::TransferFunction(const TransferFunction& tf)
{
  LOG_DBG(trans) << "TF::TF(tf) copy " << ToString();
  Copy(tf);
}

TransferFunction::~TransferFunction()
{
  LOG_DBG(trans) << "TF::~TF " << ToString() << " parser_=" << (parser_ ? "!" : "") << "NULL"
                 << " fh=" << function_handle_ << " df=" << derivative_handle_;
  if(parser_)
  {

    // parser_->ReleaseHandle(function_handle_);
    // parser_->ReleaseHandle(derivative_handle_);
    parser_ = NULL;
  }
}

TransferFunction& TransferFunction::operator=(const TransferFunction& other)
{
  LOG_DBG(trans) << "TF::=";
  Copy(other);
  return *this;
}

void TransferFunction::Copy(const TransferFunction& tf)
{
  LOG_DBG(trans) << "TF::TF(tf) copy " << ToString();

  this->application_ = tf.application_;
  this->beta_ = tf.beta_;
  this->design_ = tf.design_;
  this->param_ = tf.param_;
  this->offset_ = tf.offset_;
  this->orgType_ = tf.orgType_;
  this->scaling_ = tf.scaling_;
  this->type_ = tf.type_;
  this->parser_ = nullptr;
  this->tune_ = tf.tune_;

  if(type_ == EXPRESSION) {
    InitParser(tf.parser_->GetExpr(tf.function_handle_), tf.parser_->GetExpr(tf.derivative_handle_));
  }

}


void TransferFunction::InitParser(const string& func_expr, const string& deriv_expr)
{
  assert(function_handle_ == 0);
  assert(derivative_handle_ == 0);

  parser_ = domain->GetMathParser();
  function_handle_  = parser_->GetNewHandle();
  derivative_handle_= parser_->GetNewHandle();
  parser_->SetValue(function_handle_, "p", param_);
  parser_->SetValue(derivative_handle_, "p", param_);
  parser_->SetValue(function_handle_, "b", beta_);
  parser_->SetValue(derivative_handle_, "b", beta_);
  parser_->RegisterExternalVar(function_handle_, "rho", &expression_rho_);
  parser_->RegisterExternalVar(derivative_handle_, "rho", &expression_rho_);
  parser_->SetExpr(function_handle_, func_expr);
  parser_->SetExpr(derivative_handle_, deriv_expr);
  LOG_DBG(trans) << "TF::IP fh=" << function_handle_ << " dh=" << derivative_handle_ << " vals: " << parser_->ToString(function_handle_);
}

void TransferFunction::RegisterTune(Optimization* opt)
{
  LOG_DBG(trans) << "RT set=" << tune_.IsSet();
  if(tune_.IsSet())
    tune_.Register(&param_, opt);
}


App::Type TransferFunction::Default(const Context* ctxt)
{
  switch(ctxt->ToApp()) // will fail for piezo ?!
  {
  case App::MECH:
    return App::MECH;
  case App::BUCKLING:
    return App::BUCKLING;
  case App::ELEC:
    return App::ELEC;
  case App::HEAT:
  case App::ACOUSTIC:
    return App::LAPLACE;
  case App::LBM:
    return App::LBM;
  case App::MAG:
    return App::MAG;
  default:
    throw Exception("invalid");
  }
}

/** see the other Default */
App::Type TransferFunction::Default(DesignElement::Type type, const Context* ctxt)
{
  switch(type)
  {
  case DesignElement::DENSITY:
    if(ctxt)
      return Default(ctxt);

  case DesignElement::EMODUL:
  case DesignElement::EMODULISO:
  case DesignElement::GMODUL:
  case DesignElement::POISSON:
  case DesignElement::POISSONISO:
  case DesignElement::ROTANGLE:
  case DesignElement::ROTANGLEFIRST:
  case DesignElement::ROTANGLESECOND:
  case DesignElement::ROTANGLETHIRD:
  case DesignElement::STIFF1:
  case DesignElement::STIFF2:
  case DesignElement::STIFF3:
  case DesignElement::MECH_11:
  case DesignElement::MECH_12:
  case DesignElement::MECH_13:
  case DesignElement::MECH_14:
  case DesignElement::MECH_15:
  case DesignElement::MECH_16:
  case DesignElement::MECH_22:
  case DesignElement::MECH_23:
  case DesignElement::MECH_24:
  case DesignElement::MECH_25:
  case DesignElement::MECH_26:
  case DesignElement::MECH_33:
  case DesignElement::MECH_34:
  case DesignElement::MECH_35:
  case DesignElement::MECH_36:
  case DesignElement::MECH_44:
  case DesignElement::MECH_45:
  case DesignElement::MECH_46:
  case DesignElement::MECH_55:
  case DesignElement::MECH_56:
  case DesignElement::MECH_66:
  case DesignElement::SHEAR1:
  case DesignElement::MULTIMATERIAL:
  case DesignElement::INTERPOLATION:
    return App::MECH;

  case DesignElement::ACOU_DENSITY:
    return App::LAPLACE;

  case DesignElement::POLARIZATION:
    return App::PIEZO_COUPLING;

  case DesignElement::RHS_DENSITY:
    return App::MAG;

  default:
    throw Exception("invalid request for transfer function");
  }
}

void TransferFunction::ToInfo(PtrParamNode in, bool skip_relation) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  if(!skip_relation)
  {
    in->Get("application")->SetValue(Optimization::application.ToString(application_));
    in->Get("design")->SetValue(DesignElement::type.ToString(design_));
  }
  switch(type_)
  {
  case HEAVISIDE:
  case TANH:
  case EXPRESSION:
    in->Get("beta")->SetValue(tune_.IsSet() ? "tune" : std::to_string(beta_));
    if(tune_.IsSet())
      tune_.ToInfo(in->Get("tune"));
    // intentionally no break
  case SIMP_TYPE:
  case RAMP:
  case FIXED:
    in->Get("param")->SetValue(param_);
    break;
  default:
    break;
  }
  if(type_ == EXPRESSION)
  {
    in->Get("func")->SetValue(parser_->GetExpr(function_handle_));
    in->Get("deriv")->SetValue(parser_->GetExpr(derivative_handle_));
  }

  if(scaling_ != 1.0 || offset_ != 0.0)
  {
    assert(type_ == HEAVISIDE || type_ == TANH || type_ == EXPRESSION);
    in->Get("scaling")->SetValue(scaling_);
    in->Get("offset")->SetValue(offset_);
  }
}

void TransferFunction::Enable(bool enable)
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


bool TransferFunction::IsPenalized() const
{
  switch(type_)
  {
  case SIMP_TYPE:
    return param_ != 1.0;

  case RAMP:
    return param_ != 0.0;

  case FIXED:
  case FULL:
  case IDENTITY:
    return false;

  case HEAVISIDE:
  case TANH:
  case HASHIN_SHTRIKMAN:
  case EXPRESSION: // we assume so
    return true;

  case NO_TYPE:
    EXCEPTION("wrong");
  }

  return false; // stupid C++
}

std::string TransferFunction::ToString()
{
  std::ostringstream os;
  os << "TransferFunction[type= " << type.ToString(type_) 
     << "; application=" << Optimization::application.ToString(application_)
     << "; design=" << DesignElement::type.ToString(design_)
     << "; param=" << param_  
     << (orgType_ == NO_TYPE ? "; enabled" : "; DISABLED!") << "]";
  return os.str();   
}

double TransferFunction::Transform(const DesignElement* de, DesignElement::Access access, bool lower_bimat) const
{
  double value = de->GetValue(BaseDesignElement::DESIGN, access);
  LOG_DBG2(trans) << "T de=" << de->elem->elemNum << " a=" << de->access.ToString(access) << " p=" << de->GetPlainValue(de->DESIGN) << " v=" << value;
  return this->Transform(value, lower_bimat, de);
}



double TransferFunction::Transform(double value, bool lower_bimat, const BaseDesignElement* bde) const
{
  LOG_DBG2(trans) << "T v=" << value << " lb=" << lower_bimat << " param=" << param_ << " bde=" << (bde ? "set" : "null");

  if(lower_bimat)
  {
    value = 1.0 - value;
    if(bde) // with design simp param tune and bimaterial we have rounding issues, where value is 1+1e16 -> becomes negative
      value = std::max(bde->GetLowerBound(), value);
  }
  double result = -1.0;
  switch(type_)
  {
  case NO_TYPE: // if disabled
    result = 1.0;
    break;

  case IDENTITY:
    result = value;
    break;

  case SIMP_TYPE:
    assert(param_ >= 0);
    result = std::pow(value, param_);
    break;

  case RAMP:
    assert(param_ >= 0);
    result = value / (1.0 + param_ * (1.0 - value));
    break;

  case HASHIN_SHTRIKMAN:
    result = value / (3-2*value);
    break;

  case FIXED:
    result = param_;
    break;

  case FULL:
    assert(bde != NULL);
    result = bde->GetUpperBound();
    break;

  case HEAVISIDE:
    assert(!lower_bimat); // first check what we do!
    // some options and the derivatives
    // plot (1-exp(-20*x)), 20*x*exp(-20*x), 4*(1-exp(-10*x))**3 * 10*x*exp(-10*x), (1-exp(-10*x))**4, 1-exp(-20*x**6), 20*x**6*6*x**5*exp(-20*x**6)
    assert(beta_ >= 0.0);
    // we optionally scale the stuff when we have physical design
    result = scaling_ * (std::pow(1.0 - std::exp(-1.0 * beta_ * value), param_)) + offset_;
    break;

  case TANH:
    assert(!lower_bimat); // first check what we do!
    assert(beta_ >= 0.0);
    // tf(x) =  1 - 1/(exp(2*beta*(x-param)) + 1)
    // we optionally scale the stuff when we have physical design
    // tf(x) =  scaling * (1 - 1/(exp(2*beta*(x-param)) + 1)) + offset
    result = scaling_ * (1.0 - 1.0/(std::exp(2.0 * beta_ * ( value - param_)) + 1.0)) + offset_;
    break;

  case EXPRESSION:
    //assert(parser_ != NULL);
    #pragma omp critical
    {
      expression_rho_ = value; // via pointer this lives as rho in the expression
      // parser_->SetValue(function_handle_, "rho", value); // RegisterExternalVar does not work
      double tmp = parser_->Eval(function_handle_);
      result = offset_ + scaling_ * tmp;
      //result = function_parser_.Eval();
      LOG_DBG3(trans) << "T: v=" << value << " t=" << tmp << " -> " << result;
    }
    break;

  default:
    assert(false);
    break;
  }

  LOG_DBG2(trans) << "T v=" << value << " lbm=" << lower_bimat << " org=" << (lower_bimat ? 1.0-value : 88.88)
                  << " type=" << type.ToString(type_) << " param=" << param_ << " (" << &param_ << ") -> " << result;
  assert(!std::isnan(result));

  return result;
}


double TransferFunction::Derivative(const DesignElement* de, DesignElement::Access access, bool lower_bimat) const
{
  double value = de->GetValue(BaseDesignElement::DESIGN, access);
#ifdef CHECK_INDEX
  if(de->GetType() != design_ && (design_ == BaseDesignElement::DEFAULT && de->GetDesignSpace() != NULL && de->GetDesignSpace()->design.GetSize() > 1))
    throw Exception("type mismatch for the transfer function");
#endif
  return this->Derivative(value, lower_bimat);
}

double TransferFunction::Derivative(double value, bool lower_bimat) const
{
  switch(type_)
  {
    case NO_TYPE:
    case IDENTITY:
      if(!lower_bimat)
        return 1.0;
      else
        return -1.0;

    case SIMP_TYPE:
      if (!lower_bimat)
        return param_ * std::pow(value, param_ - 1.0);
      else
        return - param_ * std::pow(1.0-value, param_ - 1.0);

    case RAMP:
    {
      double x = value;
      if(!lower_bimat)
        return (1.0 + param_) / std::pow(1 + param_ * (1.0 - x), 2);
      else
        return -1.0 / (param_*x+1)-(param_*(1-x)) / std::pow(param_*x+1,2);
    }

    case HASHIN_SHTRIKMAN:
    {
      assert(!lower_bimat);
      double den = (3.0 - 2.0 * value);
      return (2*value)/(den*den) + 1/den;
    }

    case HEAVISIDE:
    {
      assert(!lower_bimat);
      // f = (1-hs)^hp,
      assert(beta_ > 0.0);
      double hs = std::exp(-1.0 * beta_ * value);

      return scaling_ * (param_ * std::pow(1.0 - hs, param_ - 1.0) * beta_ * hs);
    }
    case TANH:
    {
      assert(!lower_bimat);
      // tf(x)  =  1 - 1/(exp(2*beta*(x-param)) + 1)
      // tf'(x) =  (exp(2*beta*(x-param)+1)^-2 * 2 * beta * exp(2*beta*(x-param))
      double e = std::exp(2.0 * beta_ * ( value - param_));
      return scaling_ * (1.0/((e+1.0)*(e+1.0)) * 2.0 * beta_ * e);
    }

    case FIXED:
    case FULL:
      return 0.0; // the derivative of a constant is 0

    case EXPRESSION:
    {
      assert(!lower_bimat);
      double result = 0;
      #pragma omp critical
      {
        assert(parser_ != NULL);
        expression_rho_ = value;
        result = scaling_ * parser_->Eval(derivative_handle_);
      }
      return result;
    }
    default:
      assert(false);
      return -1.0;
  }
}
     
void TransferFunction::SetEnums()
{
  type.SetName("TransferFunction::Type");
  type.Add(NO_TYPE, "no_type");  
  type.Add(SIMP_TYPE, "simp");
  type.Add(IDENTITY, "identity");
  type.Add(RAMP, "ramp");
  type.Add(FIXED, "fixed");
  type.Add(FULL, "full");
  type.Add(HEAVISIDE, "heaviside");
  type.Add(TANH, "tanh");
  type.Add(HASHIN_SHTRIKMAN, "hashin_shtrikman");
  type.Add(EXPRESSION, "expression");
}     
