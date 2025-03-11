#include <stddef.h>
#include <algorithm>
#include <cassert>

#ifdef USE_4_CFS
  #include "DataInOut/Logging/LogConfigurator.hh"
  #include "General/Exception.hh"
  #include "Optimization/Optimizer/SCPIPBase.hh"
  #include "Optimization/Optimizer/scpip30.hh"
  #include "DataInOut/ProgramOptions.hh"
  using namespace CoupledField;

  DEFINE_LOG(scpip_base, "scpip_base")
#else
  #include "BasicEnum.hh"
  #include "BasicException.hh"
  #include "BasicStdVector.hh"
  #include "SCPIPBase.hh"
  #include "scpip30.hh"

  // in openCFS we use a wrapped version of the unofficial
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
  ieleng = 0;
  eqleng = 0;
  f_org = 0.0;

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

  // define constraint gradient structure
  SetConstraintSparsityPattern();

}

void SCPIPBase::SetDefaultParameters()
{
  // this are the default parameter from Sonja in ScpSolver.cc
  // ioptions[13] = {2, 1, 500, 3, 3, 0, 0, 1, 0, 0, 0, 0, 0}
  // doptions[10] = {1.0e-6, 1.0e30, 1.0e30, 1.0e-6, 1.0e-1, 1.0e-1, 0.e0, 1.0e30, 1.0e30, 1.0e30};


  icntl[1-1] = 1;   // method of moving asymptotes
  icntl[2-1] = 1;   // current svanberg
  icntl[3-1] = 100; // maximum allowed number of iterations
  icntl[4-1] = 1;   // output level -> none
  icntl[5-1] = 10;  // max function calls in line-search
  icntl[6-1] = 0;   // convergence check ->  KT
  icntl[11-1] = 0;  // warmstart
  icntl[12-1] = 0;  // ip matrix not fixed
  icntl[13-1] = 1;  // dynamic allocation - overwritten in FeasSCP

  rcntl[1-1] = 1.e-7; // final accuray of costraint violation
  rcntl[2-1] = 1.e30; // infinity
  rcntl[3-1] = 1.e30; // constraint limit
  rcntl[4-1] = 0.01;  // for relaxed convergence: relative objective value change
  rcntl[5-1] = 0.01;  // for relaxed convergence: absolut objective value change
  rcntl[6-1] = 0.01;  // for relaxed convergence: relative iteration move

  nout = 7;
}

int SCPIPBase::solve_problem(bool fromWarmstart)
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
  assert(ieleng == ielpar || ieleng == ielpar-1); // the array has always min 1 entry to hace a pointer

  int eqlpar = eqrn.GetSize();
  assert(eqlpar == (int) eqcn.GetSize());
  assert(eqlpar == (int) eqcoef.GetSize());
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

    LOG_DBG2(scpip_base) << "sp1: n = " << n << ", mie = " << mie << ", meq = " << meq << ", iemax = " << iemax << ", eqmax = " << eqmax;
    LOG_DBG3(scpip_base) << "sp1: initial guess = " << x.ToString();
    LOG_DBG3(scpip_base) << "sp1: lower bounds = " << x_l.ToString();
    LOG_DBG3(scpip_base) << "sp1: upper bounds = " << x_u.ToString();
    LOG_DBG2(scpip_base) << "sp1: function value = " << f_org;
    LOG_DBG3(scpip_base) << "sp1: inequality constraints = " << h_org.ToString();
    LOG_DBG3(scpip_base) << "sp1: equality constraints = " << g_org.ToString();
    LOG_DBG3(scpip_base) << "sp1: gradient = " << df.ToString();
    LOG_DBG3(scpip_base) << "sp1: lagrange multipliers (ineq.) = " << y_ie.ToString();
    LOG_DBG3(scpip_base) << "sp1: lagrange multipliers (eq.) = " << y_eq.ToString();
    LOG_DBG3(scpip_base) << "sp1: lagrange multipliers (l.b.) = " << y_l.ToString();
    LOG_DBG3(scpip_base) << "sp1: lagrange multipliers (u.b.) = " << y_u.ToString();
    LOG_DBG2(scpip_base) << "sp1: icntl = " << icntl.ToString();
    LOG_DBG2(scpip_base) << "sp1: rcntl = " << rcntl.ToString();
    LOG_DBG2(scpip_base) << "sp1: info = " << info.ToString();
    LOG_DBG2(scpip_base) << "sp1: rinfo = " << rinfo.ToString();
    LOG_DBG2(scpip_base) << "sp1: nout = " << nout;
    LOG_DBG3(scpip_base) << "sp1: r_scp = " << r_scp.ToString();
    LOG_DBG3(scpip_base) << "sp1: r_dim = " << rdim;
    LOG_DBG3(scpip_base) << "sp1: r_sub = " << r_sub.ToString();
    LOG_DBG3(scpip_base) << "sp1: r_subdim = " << rsubdim;
    LOG_DBG3(scpip_base) << "sp1: i_scp = " << i_scp.ToString();
    LOG_DBG3(scpip_base) << "sp1: i_dim = " << idim;
    LOG_DBG3(scpip_base) << "sp1: i_sub = " << i_sub.ToString();
    LOG_DBG3(scpip_base) << "sp1: i_subdim = " << isubdim;
    LOG_DBG3(scpip_base) << "sp1: active = " << active.ToString();
    LOG_DBG2(scpip_base) << "sp1: mode = " << mode << ", ierr = " << ierr;
    LOG_DBG3(scpip_base) << "sp1: iern = " << iern.ToString();
    LOG_DBG3(scpip_base) << "sp1: iecn = " << iecn.ToString();
    LOG_DBG3(scpip_base) << "sp1: iederv = " << iederv.ToString();
    LOG_DBG3(scpip_base) << "sp1: ielpar = " << ielpar << ", ieleng = " << ieleng;
    LOG_DBG3(scpip_base) << "sp1: eqrn = " << eqrn.ToString();
    LOG_DBG3(scpip_base) << "sp1: eqcn = " << eqcn.ToString();
    LOG_DBG3(scpip_base) << "sp1: eqcoef = " << eqcoef.ToString();
    LOG_DBG3(scpip_base) << "sp1: eqlpar = " << eqlpar << ", eqleng = " << eqleng << ", mactiv = " << mactiv;
    LOG_DBG3(scpip_base) << "sp1: spiw = " << spiw.ToString();
    LOG_DBG3(scpip_base) << "sp1: spiwdim = " << spiwdim;
    LOG_DBG3(scpip_base) << "sp1: spdw = " << spdw.ToString();
    LOG_DBG3(scpip_base) << "sp1: spdwdim = " << spdwdim;
    LOG_DBG2(scpip_base) << "sp1: spstrat = " << spstrat << ", linsys = " << linsys;
    
// we must NOT link libscpip30.a when using libscpip40i.a because of double function names.
// but USE_SCPIP needs to be enabled for USE_FEAS_SCP because it is a subclass!

#ifndef USE_FEAS_SCP
    scpip30(&n,                      // 1
             &mie,                    // 2
             &meq,                    // 3
             &iemax,                  // 4
             &eqmax,                  // 5
             x.GetPointer(),          // 6
             x_l.GetPointer(),        // 7
             x_u.GetPointer(),        // 8
             &f_org,                  // 9
             h_org.GetPointer(),      // 10
             g_org.GetPointer(),      // 11
             df.GetPointer(),         // 12
             y_ie.GetPointer(),       // 13
             y_eq.GetPointer(),       // 14
             y_l.GetPointer(),        // 15
             y_u.GetPointer(),        // 16
             icntl.GetPointer(),      // 17
             rcntl.GetPointer(),      // 18
             info.GetPointer(),       // 19
             rinfo.GetPointer(),      // 20
             &nout,                   // 21
             r_scp.GetPointer(),      // 22
             &rdim,                   // 23
             r_sub.GetPointer(),      // 24
             &rsubdim,                // 25
             i_scp.GetPointer(),      // 26
             &idim,                   // 27
             i_sub.GetPointer(),      // 28
             &isubdim,                // 29
             active.GetPointer(),     // 30
             &mode,                   // 31
             &ierr,                   // 32
             iern.GetPointer(),       // 33
             iecn.GetPointer(),       // 34
             iederv.GetPointer(),     // 35
             &ielpar,                 // 36
             &ieleng,                 // 37
             eqrn.GetPointer(),       // 38
             eqcn.GetPointer(),       // 39
             eqcoef.GetPointer(),     // 40
             &eqlpar,                 // 41
             &eqleng,                 // 42
             &mactiv,                 // 43
             spiw.GetPointer(),       // 44
             &spiwdim,                // 45
             spdw.GetPointer(),       // 46
             &spdwdim,                // 47
             &spstrat,                // 48
             &linsys);                // 49
#endif
    
    LOG_DBG3(scpip_base) << "after call to scpip: n = " << n << ", mie = " << mie << ", meq = " << meq
        << ", iemax = " << iemax << ", eqmax = " << eqmax << ", initial guess = " << x.ToString()
        << ", lower bounds = " << x_l.ToString() << ", upper bounds = " << x_u.ToString()
        << ", function value = " << f_org << ", inequality constraints = " << h_org.ToString() << ", equality constraints = " << g_org.ToString()
        << ", gradient = " << df.ToString() << ", lagrange multipliers (ineq.) = " << y_ie.ToString() << ", lagrange multipliers (eq.) = " << y_eq.ToString()
        << ", lagrange multipliers (l.b.) = " << y_l.ToString() << ", lagrange multipliers (u.b.) = " << y_u.ToString()
        << ", icntl = " << icntl.ToString() << ", rcntl = " << rcntl.ToString()
        << ", info = " << info.ToString() << ", rinfo = " << rinfo.ToString()
        << ", nout = " << nout << ", r_scp = " << r_scp.ToString() << ", r_dim = " << rdim
        << ", r_sub = " << r_sub.ToString() << ", r_subdim = " << rsubdim
        << ", i_scp = " << i_scp.ToString() << ", i_dim = " << idim
        << ", i_sub = " << i_sub.ToString() << ", i_subdim = " << isubdim
        << ", active = " << active.ToString() 
        << ", mode = " << mode << ", ierr = " << ierr 
        << ", iern = " << iern.ToString() << ", iecn = " << iecn.ToString() << ", iederv = " << iederv.ToString()
        << ", ielpar = " << ielpar << ", ieleng = " << ieleng
        << ", eqrn = " << eqrn.ToString() << ", eqcn = " << eqcn.ToString() << ", eqcoef = " << eqcoef.ToString()
        << ", eqlpar = " << eqlpar << ", eqleng = " << eqleng
        << ", mactiv = " << mactiv
        << ", spiw = " << spiw.ToString() << ", spiwdim = " << spiwdim
        << ", spdw = " << spdw.ToString() << ", spdwdim = " << spdwdim
        << ", spstrat = " << spstrat << ", linsys = " << linsys;

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

void SCPIPBase::PrepareConstraintPattern(bool initial_call)
{
  assert(m >= 0 && iemax >= 0 && eqmax >= 0);

  if(initial_call)
  {

  }

  ieleng = 0;
  int ie_cnt = 0;
  for(const ConstraintShift& cs : shift)
  {
    if(!cs.equal)
    {
      assert(active[ie_cnt] == 0 || active[ie_cnt] == 1);

      if(active[ie_cnt] == 1)
        ieleng += cs.nnz;

      ie_cnt++;
    }
    LOG_DBG(scpip_base) << "PCP cs " << cs.ToString() << " iec=" << ie_cnt << " ieleng=" << ieleng;
  }
}

void SCPIPBase::SetDenseConstraintGradientPattern()
{
   assert(iern.GetSize() != 0 && iecn.GetSize() != 0 && iederv.GetSize() != 0);
   assert(eqrn.GetSize() != 0 && eqcn.GetSize() != 0 && eqcoef.GetSize() != 0);
   assert(m >= 0 && mie >= 0 && meq >= 0);

   int ie_conseq = 0;
   for(int ie = 0; ie < mie; ie++)
   {
     // the active set wants consecutive constraint numbers which is really strange
     // as such a constraint has a different number for the constraint evaluation
     // and it's gradient!!
     assert(active[ie] == 1 || active[ie] == 0);
     if(active[ie] == 1)
     {
       for(int e = 0; e < n; e++)
       {
         int index = ie_conseq * n + e;
         iern[index] = e + 1; // fortran!
         iecn[index] = ie_conseq +1; // fortran!
       }
       ie_conseq++;
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
  rcntl.Resize(7, 0.0);
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
  for(int i = 0; i < m; i++)
    g_scaling[i] = 1.0;

  // handle dense and sparse constraint jacobians
  jac_g.Resize(nnz_jac_g, 0.0);
  jac_g_size.Resize(m);

  if(nnz_jac_g != m * n) // sparse?
  {
    assert(nnz_jac_g < m * n && nnz_jac_g > 0);
    get_sparsity_pattern_size(m, jac_g_size.GetPointer());
  }
  else // dense
  {
    for(int i = 0; i < m; i++)
      jac_g_size[i] = n;
  }

  // temporary only
  StdVector<double> g_l(m);
  StdVector<double> g_u(m);

  get_bounds_info(n, x_l.GetPointer(), x_u.GetPointer(), m, g_l.GetPointer(), g_u.GetPointer());

  shift.Resize(m);

  // now determine mie and meq and set the constraint shift and factors;
  mie = meq = 0;
  ieleng = 0; // sum the total number of nnz up
  eqleng = 0; // do it on the fly as there might be sparse jacobians
  for(int i = 0; i < m; i++)
  {
    ConstraintShift* cs = &shift[i];

    cs->number = i;
    cs->nnz = jac_g_size[i];

    cs->factor = 1.0;

    cs->lower = g_l[i];
    cs->upper = g_u[i];
    cs->equal = cs->lower == cs->upper;
    // equality constraint ?
    if(cs->equal)
    {
      meq++;
      eqleng += cs->nnz;
      cs->shift = cs->lower;
      LOG_DBG(scpip_base) << "AP: i=" << i << " -> equality : meq=" << meq << " s=" << cs->shift << " nnz=" << cs->nnz << " eqleng=" << eqleng;
    }
    else
    {
      // inequality constraint!
      mie++;
      ieleng += cs->nnz;

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
      LOG_DBG(scpip_base) << "AP: i=" << i << " -> inequality : meq=" << meq << " s=" << cs->shift << " nnz=" << cs->nnz << " ieleng=" << ieleng;
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
  active.Resize(iemax, 1);
  // starting value
  mactiv = m;

  // objective constraint
  df.Resize(n, 0.0);

  // lagrange multipliers
  y_ie.Resize(iemax, 0.0);
  y_eq.Resize(eqmax, 0.0);
  y_l.Resize(n, 0.0);
  y_u.Resize(n, 0.0);

  // take care of sparse constraint gradients
  // there is no NULL pointer, hence the minimal size is 1
  ielpar = std::max(ieleng, 1);
  iern.Resize(ielpar, 0);
  iecn.Resize(ielpar, 0);
  iederv.Resize(ielpar, 0);

  int eqlpar = std::max(eqleng, 1);
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


  int spiwdim = 1;
  int spdwdim = 1;
  
  switch(spstrat)
  {
    case 0:
      linsys = 0;
      break;
    case 1:
      switch(linsys){
      case 1:
        spiwdim = 1;
        spdwdim = std::max(4*(meq)*(meq),1);
        break;
      case 2:
        spiwdim = 136 * (ieleng + eqleng) + 32 * (mie + meq) + 3*n + 22;
        spdwdim = 70 * (ieleng + eqleng) + 2*(mie + meq);  // estimation is 10 * (ieleng + eqleng), this is too small usually
        break;
      case 3:
        spiwdim = 1;
        spdwdim = 1;
        break;
      default:
        throw Exception("SCPIP: linsys not handled");
      }
      break;
    case 2:
      switch(linsys){
      case 1:
        spdwdim = 1;
        spiwdim = n*n;
        break;
      case 2:
        throw Exception("SCPIP: spstart / linsys combination not supported");
      case 3:
        spiwdim = 1;
        spdwdim = 1;
        break;
      default:
        throw Exception("SCPIP: linsys not handled");
      }
      break;
    default:
      throw Exception("spstrat not handled");
  }
  
  spdw.Resize(spdwdim, 0.0);
  spiw.Resize(spiwdim, 0);
}

void SCPIPBase::AllocateDynamic()
{
  if((int) spiw.GetSize() != info[6-1])
  {
    if(progOpts->DoDetailedInfo())
     std::cout << "SCPIP: request to change spiwdim from " << spiw.GetSize() << " to " << info[6-1] << std::endl;
    spiw.Resize(2*info[6-1], 0); // X
    spiw.Resize(info[6-1], 0);
  }
  if((int) spdw.GetSize() != info[7-1])
  {
    if(progOpts->DoDetailedInfo())
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

int SCPIPBase::get_sparsity_pattern_size(int m, int* jac_g_dim)
{
  throw Exception("get_sparsity_pattern_size() not overloaded but get_nlp_info() reported sparse Jacobians.");
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
    LOG_DBG3(scpip_base) << "EFV: eval_g[" << i << "]: " << g_unscaled[i] << " * "
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
      LOG_DBG2(scpip_base) << "EFV: eval_g[" << c << "]: eq=" << eq << " " << g[c] << " -> " << g_org[eq];
      eq++;
    }
    else
    {
      val *= cs.factor;
      h_org[ie] = val;
      LOG_DBG2(scpip_base) << "EFV: eval_g[" << c << "]: ie=" << ie << " " << g[c] << " -> " << h_org[ie];
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

  // evaluate the grad g temporarily unsorted. Be aware of sparsity
  assert(nnz_jac_g <= m * n && (int) jac_g.GetSize() == nnz_jac_g);
  if(m > 0)
  {
    ok = eval_jac_g(n, x.GetPointer(), m, nnz_jac_g, jac_g.GetPointer());
    if(!ok) return false; // might be scaling violation for scale stuff
  }

  // the constraint scaling has nothing to do with the normalization
  if(use_g_scaling)
  {
    int idx = 0;
    for(int c = 0; c < m; c++)
    {
       double scale = g_scaling[c];
       for(int e = 0; e < jac_g_size[c]; e++)
       {
         jac_g[idx] *= scale;
         idx++;
       }
    }
  }

  grad_evals++;

  int ie = 0;
  int eq = 0;
  (void) eq; // "use" variable to avoid unsued warning as we only use it in an assert
  int ie_nnz = 0; // counter
  int eq_nnz = 0;
  for(unsigned int c = 0; c < shift.GetSize(); c++)
  {
    double* own = jac_g.GetPointer() + ie_nnz + eq_nnz;
    ConstraintShift& cs = shift[c];
    if(cs.equal)
    {
      CopyConstraintGradient(own, eqcoef.GetPointer() + eq_nnz, cs);
      eq++;
      eq_nnz += cs.nnz;
    }
    else
    {
      if(active[ie] == 1)
      {
        CopyConstraintGradient(own, iederv.GetPointer() + ie_nnz, cs);
        ie_nnz += cs.nnz;
      }
      ie++;
    }
    LOG_DBG(scpip_base) << "EG: c=" << c << " " << (cs.equal ? "equal" : "inequality")
                        << " nnz=" << cs.nnz << " active=" << (cs.equal ? -1 : active[ie-1]);
  }

  assert(ie + eq == m);
  return true;
}

void SCPIPBase::CopyConstraintGradient(const double* ipopt, double* scpip, ConstraintShift& cs)
{
  for(int i = 0; i < cs.nnz; i++)
  {
    double val = ipopt[i];
    // the shift is not in the gradient, but maybe the inequality normalization
    val *= cs.factor;
    scpip[i] = val;
    LOG_DBG3(scpip_base) << "CCG: grad_g[" << cs.number << "]:" << i << " " << ipopt[i] << " -> " << scpip[i];
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


std::string SCPIPBase::ConstraintShift::ToString() const
{
  std::ostringstream os;
  os << "n=" << number << " s=" << shift << "nnz=" << nnz;
  return os.str();
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

    case 6:   os << "For at least one component, the lower bound is greater or equal to the upper bound. Equal disables the interior point subproblem solver.";
              break;

    case 9:   os << "rsubdim (" << r_sub.GetSize() << ") too small. Required: " << info[4-1];
              break;

    case 15: os << "The Jacobian matrices are not stored columnwise";
              break;

    case 16:  os << "The jacobian matrices are not stored correctly. For at least one "
                 << "column the components are out of order";

              for(int i = 0; i < ieleng; i++)
                std::cout << "iern[" << i << "]=" << iern[i]
                          << "\tiecn[" << i << "]=" << iecn[i]
                          << "\tiederv[" << i << "]=" << iederv[i] << std::endl;

              for(int i = 0; i < eqleng; i++)
                std::cout << "eqrn[" << i << "]=" << eqrn[i]
                          << "\teqcn[" << i << "]=" << eqcn[i]
                          << "\teqcoef[" << i << "]=" << eqcoef[i] << std::endl;

              break;
              
    case LineSearch_Max_Iter:
              os << "Linesearch could not find a descent direction within maximum allowed iterations!";
              break;

    case Infeasible:
              os << "The norm of the gradient of the Lagrangian is close to 0 and the "
                 << "maximum of the artificial variables is greater or equal 1. Together it is very likely,"
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

    case 161: os << "The Jacobian matrices are not stored correctly. More columns than constraints are provided.";
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
  LOG_DBG(scpip_base) << "SSV key=" << key << " idx=" << idx << " value=" << value << " ival=" << ival;
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

  LOG_DBG(scpip_base) << "SIV key=" << key << " idx=" << idx << " value=" << value;
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
  LOG_DBG(scpip_base) << "SNV key=" << key << " idx=" << idx << " value=" << value;
}


void SCPIPBase::SetEnums()
{
  icntl_.SetName("ICNTL");
  icntl_.Add(1-1, "optimization_method");
  icntl_.Add(2-1, "asymptotes_strategy");
  icntl_.Add(3-1, "max_iter");
  icntl_.Add(4-1, "output_level");
  icntl_.Add(5-1, "maximum_linesearch_function_calls");
  icntl_.Add(6-1, "convergence_criteria");
  icntl_.Add(7-1, "linesearch_queue");
  icntl_.Add(8-1, "icntl_8");
  icntl_.Add(11-1, "warmstart");
  icntl_.Add(12-1, "fix_interior_point_matrix");

  rcntl_.SetName("RCNTL");
  rcntl_.Add(1-1, "max_kuhn_tucker_constraint_violation");
  rcntl_.Add(2-1, "infinity");
  rcntl_.Add(3-1, "inequality_constraint_active_limit");
  rcntl_.Add(4-1, "relaxed_relative_objective_change");
  rcntl_.Add(5-1, "relaxed_absolut_objective_change");
  rcntl_.Add(6-1, "relaxed_relative_iteration_move");
  rcntl_.Add(7-1, "max_constraint_violation");

  opt_meth_.SetName("desired optimization method");
  opt_meth_.Add(1, "moving_asymptotes");
  opt_meth_.Add(2, "sequential_convex_programming");

  output_level_.SetName("desired output level");
  output_level_.Add(1, "no_output");
  output_level_.Add(2, "only_final_convergence_analysis");
  output_level_.Add(3, "one_line_of_intermediate_results");
  output_level_.Add(4, "more_detailed_and_intermediate_results");
  output_level_.Add(5, "additional_intermediate_results");

  converge_.SetName("convergence criteria");
  converge_.Add(0, "kuhn_tucker");
  converge_.Add(1, "relaxed");
}

