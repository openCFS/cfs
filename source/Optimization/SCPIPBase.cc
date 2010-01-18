#include <cassert>

#ifdef USE_4_CFS
  #include "Optimization/scpip30.hh"
  #include "Optimization/SCPIPBase.hh"
  #include "General/exception.hh"
  #include "DataInOut/Logging/cfslog.hh"

  using namespace CoupledField;

  DECLARE_LOG(scpip_base)
  DEFINE_LOG(scpip_base, "scpip_base")
#else
  #include "scpip30.hh"
  #include "SCPIPBase.hh"
  #include "BasicException.hh"
  #include "BasicEnum.hh"
  #include "BasicStdVector.hh"

  // in CFS++ we use a wrapped version of the unofficial
  // Boost Logging Template library by John Torjo.
  // For C++SCPIP this is simply disabled.
  struct nullstream:
    std::ostream {
    nullstream(): std::ios(0), std::ostream(0) {}
  };

  nullstream nil;

  #define LOG_TRACE(log_name)   std::cout // ((void)0)
  #define LOG_TRACE2(log_name)  std::cout
  #define LOG_DBG(log_name)     nil
  #define LOG_DBG2(log_name)    nil
  #define LOG_DBG3(log_name)    nil
#endif

SCPIPBase::SCPIPBase()
{
  // set enums so Set*Value() works
  SetEnums();
  call_scale_parameters_ = false;
  use_obj_scaling = false;
  use_g_scaling = false;
  obj_scaling = 1.0;

  f_evals = grad_evals = 0;

  // --- scpip settings

  // we do inverse communication
  mode = 2;

  // from the example
  spstrat = 1;

  // set to defined values so we can check if nothing is forgotten */
  InitVariables();

  // dense cholesky solver = 1, sparse chol solver = 2
  linsys = 1;

  // allocate the fixed part tho we can set the default parameters.
  // the problem part cannot be allocated here as we have to call
  // abstract methods and the childs constructor is not executed yet.
  AllocateFixed();

  // set the default parameters here so Set*Value() can be called savely
  SetDefaultParameters();

  statistic_ = new Statistic(this);
}

SCPIPBase::~SCPIPBase()
{
  if(statistic_ != NULL) { delete statistic_; statistic_ = NULL; }
}

void SCPIPBase::Initialize()
{
  // set the problem dependend data
  AllocateProblem();

  // define dense constraint gradient structure
  SetDenseConstraintGradient();
}

void SCPIPBase::SetDefaultParameters()
{
  icntl[1-1] = 1;   // method of moving asymptotes
  icntl[2-1] = 1;   // current svanberg
  icntl[3-1] = 100; // maximum allowed number of iterations
  icntl[4-1] = 1;   // output level -> none
  icntl[5-1] = 10;  // max function calls in line-search
  icntl[6-1] = 0;   // convergence check ->  KT
  icntl[11-1] = 0;  // warmstart
  icntl[12-1] = 0;  // ip matrix not fixed
  icntl[13-1] = 1;  // dynamic allocation

  rcntl[1-1] = 1.e-7; // final accuray of costraint violation
  rcntl[2-1] = 1.e30; // infinity
  rcntl[3-1] = 1.e30; // constraint limit
  rcntl[4-1] = 0.01;  // for relaxed convergence: relative objective value change
  rcntl[5-1] = 0.01;  // for relaxed convergence: absolut objective value change
  rcntl[6-1] = 0.01;  // for relaxed convergence: relative iteration move

  nout = 7;
}

int SCPIPBase::SolveProblem(bool fromWarmstart)
{
  // set start parameters here to allow restarted SolveProblem()
  get_starting_point(n, x.GetPointer());

  // do we start from a warmstart file?
  ierr = fromWarmstart ? -4 : 0;

  // scpip counts its iteration in info[20-1].
  int last_iter = -1;

  // the dimensions are kept in the StdVector size, we don't make them class attributes.
  // here we allocate variables to give their pointers to fortran
  int ielpar = iern.GetSize();
  assert(ielpar == (int) iecn.GetSize());
  assert(ielpar == (int) iederv.GetSize());
  int ieleng = n * mie;
  assert(ieleng == ielpar || ieleng == ielpar-1);

  int eqlpar = eqrn.GetSize();
  assert(eqlpar == (int) eqcn.GetSize());
  assert(eqlpar == (int) eqcoef.GetSize());
  int eqleng = n * meq;
  assert(eqleng == eqlpar || eqleng == eqlpar-1);

  int rdim = r_scp.GetSize();
  int rsubdim = r_sub.GetSize();
  int idim = i_scp.GetSize();
  int isubdim = i_sub.GetSize();

  // is 'nlp_scaling_method' set to 'user-scaling' ? -> all but m are out parameters
  if(call_scale_parameters_)
  {
    double old_scaling = obj_scaling;
    get_scaling_parameters(obj_scaling, use_g_scaling, m, g_scaling.GetPointer());
    // check if there was already 'obj_scaling_factor' and assert it has not changed
    if(use_obj_scaling && obj_scaling != old_scaling)
      throw "There was another objective scaling in 'obj_scaling_factor' than in get_scaling_parameters()";
    use_obj_scaling = true;

    // scale the shifts
    for(int i = 0; i < m; i++)
    {
      LOG_DBG(scpip_base) << " shift[" << i << "].shift -> " << shift[i].shift << " * "
                          << g_scaling[i] << " = " << shift[i].shift * g_scaling[i];
      shift[i].shift *= g_scaling[i];
    }
  }
  
  for(;;) // we break out in success and error
  {
    // this might be dynamic!
    int spiwdim = spiw.GetSize();
    int spdwdim = spdw.GetSize();

    scpip30_(&n, &mie, &meq, &iemax, &eqmax, x.GetPointer(), x_l.GetPointer(), x_u.GetPointer(),
             &f_org, h_org.GetPointer(), g_org.GetPointer(), df.GetPointer(),
             y_ie.GetPointer(), y_eq.GetPointer(), y_l.GetPointer(), y_u.GetPointer(),
             icntl.GetPointer(), rcntl.GetPointer(),info.GetPointer(), rinfo.GetPointer(),
             &nout, r_scp.GetPointer(), &rdim, r_sub.GetPointer(),
             &rsubdim, i_scp.GetPointer(), &idim, i_sub.GetPointer(), &isubdim, active.GetPointer(),
             &mode, &ierr, iern.GetPointer(), iecn.GetPointer(), iederv.GetPointer(), &ielpar, &ieleng,
             eqrn.GetPointer(), eqcn.GetPointer(), eqcoef.GetPointer(), &eqlpar, &eqleng,
             &mactiv, spiw.GetPointer(), &spiwdim, spdw.GetPointer(), &spdwdim, &spstrat, &linsys);

    LOG_TRACE(scpip_base) << "scpip30 returns: ierr=" << ierr << " info[20-1]=" << info[20-1];

    // call callback for iteration. This might cause an user break!
    if(!intermediate_callback(info[20-1], info[20-1] != last_iter))
      ierr = User_Requested_Stop;
    last_iter = info[20-1];

    if(ierr == 0)
    {
      CallFinalizeSolution();
      break;
    }
    if(ierr == -1)
    {
      LOG_DBG2(scpip_base) << "x: " << x.ToString();
      EvaluateFunctionValues();
    }
    if(ierr == -2)
    {
      LOG_DBG2(scpip_base) << "x: " << x.ToString();
      if(!EvaluateGradients())
      {
        ierr = Gradients_Return_False;
        CallFinalizeSolution(); // allow comit of the last iteration to reuse it
        break; // no output for resart!
      }
      LOG_DBG2(scpip_base) << "df: " << df.ToString();
    }
    if(ierr == -3)
    {
      AllocateDynamic();
    }
    if(ierr == -4)
    {
      // we came from a warmstart file run an continue with -5
      ierr = -5;
    }
    if(ierr < -5 || ierr > 0)
    {
      std::cerr << ToString(ierr) << std::endl;
      CallFinalizeSolution();
      break;
    }
  }

  return ierr;
}

void SCPIPBase::InitVariables()
{
  n = m = mie = meq = iemax = eqmax = -1;
}

void SCPIPBase::SetDenseConstraintGradient()
{
   assert(iern.GetSize() != 0 && iecn.GetSize() != 0 && iederv.GetSize() != 0);
   assert(eqrn.GetSize() != 0 && eqcn.GetSize() != 0 && eqcoef.GetSize() != 0);
   assert(m >= 0 && mie >= 0 && meq >= 0);

   for(int ie = 0; ie < mie; ie++)
   {
     for(int e = 0; e < n; e++)
     {
       int index = ie * n + e;
       iern[index] = e + 1; // fortran!
       iecn[index] = ie +1; // fortran!
     }
   }

   for(int eq = 0; eq < meq; eq++)
   {
     for(int e = 0; e < n; e++)
     {
       int index = eq * n + e;
       eqrn[index] = e + 1; // fortran!
       eqcn[index] = eq +1; // fortran!
     }
   }
}

void SCPIPBase::AllocateFixed()
{
  // these are constant arrays. We use StdVector() to have index checks
  icntl.Resize(13, 0);
  rcntl.Resize(6, 0.0);
  info.Resize(23, 0);
  rinfo.Resize(5, 0.0);
}

void SCPIPBase::AllocateProblem()
{
  // get the basic information - mie and meq are not known yet.
  get_nlp_info(n, m, nnz_jac_g);

  // design space
  x.Resize(n, 0.0);

  // lower and upper bound of design space
  x_l.Resize(n, 0.0);
  x_u.Resize(n, 0.0);

  // initialize the constraint stuff
  assert(n != -1);
  assert(m != -1);

  // find mie and meq which is sorted in SCPIP but not in IPOPT interface
  g.Resize(m, 0.0);
  g_unscaled.Resize(m, 0.0);
  y_g.Resize(m, 0.0);

  // this are our constraint scalings, the initial value is replaced
  // via get_scaling_parameters()
  g_scaling.Resize(m);
  for(int i = 0; i < m; i++) g_scaling[i] = 1.0;

  assert(nnz_jac_g == m * n); // dense!!
  jac_g.Resize(nnz_jac_g, 0.0);

  // temporary only
  StdVector<double> g_l(m);
  StdVector<double> g_u(m);

  get_bounds_info(n, x_l.GetPointer(), x_u.GetPointer(), m, g_l.GetPointer(), g_u.GetPointer());

  shift.Resize(m);

  // now determine mie and meq and set the constraint shif and factors;
  mie = meq = 0;
  for(int i = 0; i < m; i++)
  {
    ConstraintShift* cs = &shift[i];

    cs->number = i;
    cs->factor = 1.0;

    cs->lower = g_l[i];
    cs->upper = g_u[i];
    cs->equal = cs->lower == cs->upper;
    // equality constraint ?
    if(cs->equal)
    {
      meq++;
      cs->shift = cs->lower;
      LOG_TRACE2(scpip_base) << "constraint " << i << " becomes equality constraint " << meq << " with shift " << cs->shift;
    }
    else
    {
      // inequality constraint!
      mie++;

      // one value will be +/- infinity -> search the smaller abs-value
      if(abs(cs->lower) < abs(cs->upper))
      {
        cs->shift = abs(cs->lower); // x < c -> x - c < 0
        cs->factor = -1.0;
      }
      else
      {
        cs->shift = abs(cs->upper);
      }
      LOG_TRACE2(scpip_base) << "constraint " << i << " becomes inequality constraint " << meq << " with shift " << cs->shift << " and factor " << cs->factor;
    }
  }
  assert(meq + mie == m);

  // 1 is minimum for SCPIP
  iemax = mie > 1 ? mie : 1;
  eqmax = meq > 1 ? meq : 1;

  // the constraint values
  h_org.Resize(iemax, 0.0);
  g_org.Resize(eqmax, 0.0);

  // inequality constraints considered active
  active.Resize(iemax, 0);
  // starting value
  mactiv = m;

  // objective constraint
  df.Resize(n, 0.0);

  // diverse lagrange multipliers
  y_ie.Resize(iemax, 0.0);
  y_eq.Resize(eqmax, 0.0);
  y_l.Resize(n, 0.0);
  y_u.Resize(n, 0.0);

  // we have only dense gradients
  // there is no NULL pointer, hence the minimal size is 1
  int ielpar = std::max(n * mie, 1);
  iern.Resize(ielpar, 0);
  iecn.Resize(ielpar, 0);
  iederv.Resize(ielpar, 0);

  int eqlpar = std::max(n * meq, 1);
  eqrn.Resize(eqlpar, 0);
  eqcn.Resize(eqlpar, 0);
  eqcoef.Resize(eqlpar, 0);

  // "fixed" Working arrays
  //           30*N+11*IEMAX+8+10*EQMAX
  r_scp.Resize(30*n+11*iemax+8+10*eqmax, 0.0);
  //           22*N+41*IEMAX+27*EQMAX+2*IELPAR+EQLPAR
  r_sub.Resize(22*n+41*iemax+27*eqmax+2*ielpar+eqlpar, 0.0);
  //           5*N+5*IEMAX+2*EQMAX+3
  i_scp.Resize(5*n+5*iemax+2*eqmax+3, 0);
  //           2*N+3*IEMAX+2*EQMAX+IELPAR
  i_sub.Resize(2*n+3*iemax+2*eqmax+ielpar, 0);

  // "dynamic" working arrays. We assume linsys=1 or 2
  assert(linsys == 1 || linsys == 2);
  int spiwdim = 1;
  int spdwdim = -1;
  
  switch(spstrat)
  {
    case 1:
      spdwdim = std::max((mactiv + meq) * (mactiv + meq), 1);
      if(linsys == 2)
      {
        spiwdim = 5 * (n * mie + n * meq) + 2 * n + 1;
        spdwdim = 10 * n * (mie + meq);
      }
      break;
    case 2:
      // scpip does not support spstrat == 2 && linsys == 2
      // assert(linsys == 1);
      spdwdim = n * n; // big!!
      break;
    default:
      throw Exception("spstrat not handled");
  }
  
  // set at least to one, or else we cannot call GetPointer() in call to scpip
  if(spdwdim <= 0) spdwdim = 1;
  
  spdw.Resize(spdwdim, 0.0);
  spiw.Resize(spiwdim, 0);
}

void SCPIPBase::AllocateDynamic()
{
  if((int) spiw.GetSize() != info[6-1])
  {
    std::cout << "SCPIP: request to change spiwdim from " << spiw.GetSize() << " to " << info[6-1] << std::endl;
    spiw.Resize(info[6-1], 0);
  }
  if((int) spdw.GetSize() != info[7-1])
  {
    std::cout << "SCPIP: request to change spdwdim from " << spdw.GetSize() << " to " << info[7-1] << std::endl;
    spdw.Resize(info[7-1], 0.0);
    if(info[7-1] > 1e6)
    {
      // this might be a bad case that cost a lot of computational time, so we issue a warning
      std::cout << "SCPIP: WARNING: this problem might take a very long time to compute" << std::endl;
      std::cout << "SCPIP: WARNING: insert a trivially fulfilled constraint to gain speed" << std::endl;
    }
  }
}

bool SCPIPBase::intermediate_callback(int iter, bool next_iter)
{
  DumpActive();

  if(next_iter) {
     std::cout << std::endl;
     PrintInfo(std::cout);
  }

  return true;
}


void SCPIPBase::EvaluateFunctionValues()
{
  f_evals++;

  eval_f(n, x.GetPointer(), f_org_unscaled);
  f_org = f_org_unscaled * obj_scaling;

  LOG_DBG2(scpip_base) << "eval_f: " << f_org_unscaled << " * " << obj_scaling
                       << " -> " << f_org;

  // now the constraints. These are sorted for equality and inequality in
  // SCPIP and normalized. We don't have this in the IPOPT frontend
  if(m > 0)
    eval_g(n, x.GetPointer(), m, g_unscaled.GetPointer());

  // the constraint scaling is done before the normalization shifting
  // if we do not scale it is set to 1.0
  for(int i = 0; i < m; i++)
  {
    g[i] = g_unscaled[i] * g_scaling[i];
    LOG_DBG3(scpip_base) << "eval_g[" << i << "]: " << g_unscaled[i] << " * "
                         <<  g_scaling[i] << " -> " << g[i];
  }

  // do the shifting!
  int ie = 0;
  int eq = 0;
  for(unsigned int c = 0; c < shift.GetSize(); c++)
  {
    ConstraintShift& cs = shift[c];
    // the cs.shift is constraint scaled in SetSolver()
    double val = g[c] - cs.shift;

    if(cs.equal)
    {
      g_org[eq] = val;
      LOG_DBG2(scpip_base) << "eval_g[" << c << "]: eq=" << eq << " " << g[c] << " -> " << g_org[eq];
      eq++;
    }
    else
    {
      val *= cs.factor;
      h_org[ie] = val;
      LOG_DBG2(scpip_base) << "eval_g[" << c << "]: ie=" << ie << " " << g[c] << " -> " << h_org[ie];
      ie++;
    }
  }
  assert(ie + eq == m);
}

bool SCPIPBase::EvaluateGradients()
{
  assert(n == (int) x.GetSize() && n == (int) df.GetSize());
  bool ok = eval_grad_f(n, x.GetPointer(), df.GetPointer());
  if(!ok) return false; // could not evaluate, maybe because scaling out of range

  if(use_obj_scaling)
    for(int i = 0; i < n; i++) df[i] *= obj_scaling;

  // evaluate the grad g temporarily unsorted
  assert(nnz_jac_g == m * n && (int) jac_g.GetSize() == nnz_jac_g);
  if(m > 0)
  {
    ok = eval_jac_g(n, x.GetPointer(), m, nnz_jac_g, jac_g.GetPointer());
    if(!ok) return false; // might be scaling violation for autoscale stuff
  }

  // the constraint scaling has nothing to do with the normalization
  if(use_g_scaling)
  {
    for(int c = 0; c < m; c++)
    {
       double scale = g_scaling[c];
       for(int e = c * n; e < (c+1) * n; e++)
         jac_g[e] *= scale;
    }
  }

  grad_evals++;

  int ie = 0;
  int eq = 0;
  for(unsigned int c = 0; c < shift.GetSize(); c++)
  {
    double* ipopt = jac_g.GetPointer() + c * n;
    ConstraintShift& cs = shift[c];
    if(cs.equal)
    {
      CopyConstraintGradient(ipopt, eqcoef.GetPointer() + eq * n, cs);
      eq++;
    }
    else
    {
      CopyConstraintGradient(ipopt, iederv.GetPointer() + ie * n, cs);
      ie++;
    }
  }
  assert(ie + eq == m);
  return true;
}

void SCPIPBase::CopyConstraintGradient(const double* ipopt, double* scpip, ConstraintShift& cs)
{
  for(int i = 0; i < n; i++)
  {
    double val = ipopt[i];
    // the shift is not in the gradient, but maybe the inequality normalization
    val *= cs.factor;
    scpip[i] = val;
    LOG_DBG3(scpip_base) << "grad_g[" << cs.number << "]:" << i << " " << ipopt[i] << " -> " << scpip[i];
  }
}

void SCPIPBase::CallFinalizeSolution()
{
  // collect the data from y_ie and y_eq into y_g.
  int ie = 0;
  int eq = 0;
  for(unsigned int i = 0; i < shift.GetSize(); i++)
    y_g[i] = shift[i].equal ? y_eq[eq++] : y_ie[ie++];
  assert(m == ie + eq);

  // call virtual method
  finalize_solution(ierr, n, x.GetPointer(), y_l.GetPointer(), y_u.GetPointer(), m,
                    g_unscaled.IsEmpty() ? NULL : g_unscaled.GetPointer(),
                    y_g.IsEmpty() ? NULL : y_g.GetPointer(), f_org_unscaled);

}




std::string SCPIPBase::ToString(int ierr)
{
  std::ostringstream os;
  switch(ierr)
  {
    case Solve_Succeeded:
              os << "solution obtained";
              break;

    case -1:  os << "reverse communication: function values are requested";
              break;

    case -2:  os << "reverse communication: gradients are requested";
              break;

    case -3:  os << "reverse communication: current values for dynamic allocation";
              break;

    case Maximum_Iterations_Exceeded:
              os << "maximum number of iterations reached (" << (info[20-1]+1) << ")";
              break;

    case 9:   os << "rsubdim (" << r_sub.GetSize() << ") too small. Required: " << info[4-1];
              break;

    case 16:  os << "The jacobian matrices are not stored correctly. For at least one "
                 << "column the components are out of order";

              assert((int) iecn.GetSize() == mie * n);
              for(int i = 0; i < mie * n; i++)
                std::cout << "iern[" << i << "]=" << iern[i]
                          << "\tiecn[" << i << "]=" << iecn[i]
                          << "\tiederv[" << i << "]=" << iederv[i] << std::endl;

              assert((int) eqcn.GetSize() == meq * n);
              for(int i = 0; i < meq * n; i++)
                std::cout << "eqrn[" << i << "]=" << eqrn[i]
                          << "\teqcn[" << i << "]=" << eqcn[i]
                          << "\teqcoef[" << i << "]=" << eqcoef[i] << std::endl;

              break;
              
    case LineSearch_Max_Iter:
              os << "Linesearch could not find a descent direction within maximum allowed iterations!";
              break;

    case 23:  os << "The norm of the gradient of the lagrangian is close to 0 and the "
                 << "maximum of the artificial variables is >= 1. Together it is very likely,"
                 << "that the feasible region is empty!";
              break;

    case Subproblem_Max_Iter:
              os << "The subproblem could not be solved within maximum allowed iterations!";
              break;

    case User_Requested_Stop:
              os << "User Requested Stop: not by SCPIP but own code";
              break;

    case Gradients_Return_False:
              os << "Evaluation of the gradients failed";
              break;

    default:  os << "error: " << ierr;
  }

  return os.str();
}

void SCPIPBase::DumpActive()
{
  for(unsigned int i = 0; i < active.GetSize(); i++)
    std::cout << "active[" << i << "] = " << active[i] << std::endl;

  std::cout << "mactiv = " << mactiv << std::endl;
}

void SCPIPBase::PrintInfo(std::ostream& os)
{
  os << "current iteration number: " << info[20-1] << std::endl;
  os << "number of evaluations of lagrangian function values: " << info[1-1] << std::endl;
  os << "number of evaluations of lagrangian gradients: " << info[2-1] << std::endl;
  os << "number of iterations for the solution of the last subproblem: " << info[21-1] << std::endl;
  os << "actually chosen spstrat: " << info[22-1] << std::endl;
  os << "actually chosen linsys: " << info[23-1] << std::endl;
  os << "residual of the subproblem of the last main iteration: " << rinfo[1-1] << std::endl;
  os << "maximum violation of constraints: " << rinfo[2-1] << std::endl;
  os << "stepsize in last main iteration: " << rinfo[3-1] << std::endl;
  os << "norm of the difference of the last two iteration points: " << rinfo[4-1] << std::endl;
  os << "norm of the gradient of the lagrangian w.r.t. x at last iteration: " << rinfo[5-1] << std::endl;
}

void SCPIPBase::SetStringValue(const std::string& key, const std::string& value)
{
  // first we check for ipop mimicry settings, then we do the semi automated scpip
  // parameter handling
  if(key == "nlp_scaling_method")
  {
    if(value == "user-scaling")
    {
      call_scale_parameters_ = true;
      return;
    }
    else throw value + " is not valid for " + key + ". Know only 'user-scaling'.";
  }

  // is either scpip parameter or a wrong one

  // icntl_ knows the C-position in icntl parameter block
  int idx = icntl_.Parse(key);

  // find the int value behind value
  int ival;
  switch(idx)
  {
    case 1-1: ival = opt_meth_.Parse(value);
              break;

    case 4-1: ival = output_level_.Parse(value);
              break;

    case 6-1: ival = converge_.Parse(value);
              break;

    default: throw Exception(key + " not handled as string value");
  }

  icntl[idx] = ival;
}

void SCPIPBase::SetIntegerValue(const std::string& key, int value)
{
  if(key == "spstrat")
  {
    spstrat = value;
    return;
  }
  if(key == "linsys")
  {
    linsys = value;
    return;
  }

  // icntl has the index of the parameter in C-style
  int idx = icntl_.Parse(key);
  icntl[idx] = value;
}

void SCPIPBase::SetNumericValue(const std::string& key, double value)
{
  // check first for ipopt mimicry
  if(key == "obj_scaling_factor")
  {
    use_obj_scaling = true;
    obj_scaling = value;
    if(value == 0.0) throw "It makes no sense to scale the objective by 0.0";
    return;
  }

  // rcntl has the index of the parameter in C-style
  int idx = rcntl_.Parse(key);
  rcntl[idx] = value;
}


void SCPIPBase::SetEnums()
{
  icntl_.SetName("ICNTL");
  icntl_.Add(1-1, "optimization_method");
  icntl_.Add(3-1, "max_iter");
  icntl_.Add(4-1, "output_level");
  icntl_.Add(5-1, "maximum_linesearch_function_calls");
  icntl_.Add(6-1, "convergence_criteria");
  icntl_.Add(11-1, "warmstart");

  rcntl_.SetName("RCNTL");
  rcntl_.Add(1-1, "max_kuhn_tucker_constraint_violation");
  rcntl_.Add(2-1, "infinity");
  rcntl_.Add(3-1, "inequality_constraint_active_limit");
  rcntl_.Add(4-1, "relaxed_relative_objective_change");
  rcntl_.Add(5-1, "relaxed_absolut_objective_change");
  rcntl_.Add(6-1, "relaxed_relative_iteration_move");

  opt_meth_.SetName("desired optimization method");
  opt_meth_.Add(1, "moving_asymptotes");
  opt_meth_.Add(2, "sequential_convex_programming");

  output_level_.SetName("desired output level");
  output_level_.Add(1, "no_output");
  output_level_.Add(2, "only_final_convergence_analysis");
  output_level_.Add(3, "one_line_of_intermediate_results");
  output_level_.Add(4, "more_detailed_and_intermediate_results");

  converge_.SetName("convergence criteria");
  converge_.Add(0, "kuhn_tucker");
  converge_.Add(1, "relaxed");
}

