#include "Optimization/Optimizer/FeasSCP.hh"
#include "Optimization/Optimizer/SCPIPBase.hh"
#include "Optimization/Optimizer/scpip40i.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"

#include <boost/algorithm/string.hpp>
#include <cstring>
#include <time.h>


DECLARE_LOG(feas_scp)
DEFINE_LOG(feas_scp, "feas_scp")

using namespace CoupledField;


// static int default_ioptions[13] = {2, 1, 500, 3, 3, 0, 0, 1, 0, 0, 0, 0, 0}
// static double default_doptions[10] = {1.0e-6, 1.0e30, 1.0e30, 1.0e-6, 1.0e-1, 1.0e-1, 0.e0, 1.0e30, 1.0e30, 1.0e30};

// constructor
FeasSCP::FeasSCP(Optimization* opt, PtrParamNode pn, Optimization::Optimizer type) : SCPIP(opt, pn, type)
{
  // handle the ipopt options simple, we might remove it anyway soon when using the C++ Ipopt!
  for(int i = 0; i < 200; i++)
    ipopt_cargs[i] = ' ';

  // set the ipopt subsolver parameters
  ParamNodeList ops = pn->Has("feasSCP") ? pn->Get("feasSCP")->GetList("ipopt_option") : ParamNodeList();
  ipopt_args.Resize(std::max((int) ops.GetSize(), 1), 0.0); // minimum size to have a valid pointer!

  for(unsigned int i = 0; i < ops.GetSize(); i++)
  {
    ipopt_args[i] = ops[i]->Get("value")->As<double>();
    std::string name = ops[i]->Get("name")->As<std::string>();
    assert(ops.GetSize() <= 9);
    // in the fortran file we extract the names and trim them as it is non-trivial
    // to pass an array of c-strings to fortran :(
    strncpy(ipopt_cargs + i * 20, name.c_str(), name.length()); // yes, we can do pointer arithmetic :)
  }
}

FeasSCP::~FeasSCP()
{

}


void FeasSCP::ToInfo(PtrParamNode pn)
{
  SCPIP::ToInfo(pn);
  for(unsigned int i = 0; i < ipopt_args.GetSize(); i++)
  {
    PtrParamNode pn_ = pn->Get("ipopt", ParamNode::APPEND);
    char full[] = "                    "; // 20 characters + tailing NULL
    strncpy(full, ipopt_cargs + 20 * i, 20);
    std::string name(full);
    boost::algorithm::trim(name);
    pn_->Get("name")->SetValue(name);
    pn_->Get("value")->SetValue(ipopt_args[i]);
  }
}

void FeasSCP::AllocateProblem()
{
  SCPIPBase::AllocateProblem();

  // ie_idx is not set yet!
  assert(m == optimization->constraints.view->GetNumberOfActiveConstraints());
  linear.Reserve(m); // better waste some space than expensive push backs
  mf_idx.Reserve(m);
  for(int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    if(g->GetBound() != Condition::EQUAL)
    {
      linear.Push_back(g->IsLinear() ? 0 : 1);

      if(!g->IsFeasibilityConstraint() && !mf_idx.IsEmpty())
        throw Exception("Feasibility constraints need to be ordered at the end");

      if(g->IsFeasibilityConstraint())
        mf_idx.Push_back(i);
    }
  }
  optimization->constraints.view->Done(); // mandatory after traversing the view

  // now we now mf
  mf = mf_idx.GetSize();
  // correct mie
  assert(mie >= mf);
  mie -= mf;

  if(mf > 0 && order_ != BY_ELEMENT)
    throw Exception("For feasibility constraints feasSCP/order needs to 'by_element'");

  LOG_DBG2(feas_scp) << "PI: linear=" << linear.ToString();
  LOG_DBG2(feas_scp) << "PI: mf_idx=" << mf_idx.ToString();

  // overwrite the default value to fixed memory allocation - suggested by Sonja
  // enablubg this stuff valgrind reports problems!
  // icntl[13-1] = 1;
  // so we need to provide a lot of memory
  // spiw.Resize(2*n+5*iemax+6);
  // spdw.Resize(iemax * iemax);

  // resize the stuff according to ScpSolver.cpp
  ielpar = n * mie + (1 + mf) * 9;
  iern.Resize(ielpar, 0);
  iecn.Resize(ielpar, 0);
  iederv.Resize(ielpar, 0);

  r_scp.Resize(44*n+16*iemax+2*ielpar+20+2*mie +100*mie, 0.0);
  r_sub.Resize(22*n+41*iemax+2*ielpar+28 +100, 0.0);
  i_scp.Resize(7*n+8*mie+3*ielpar+12+100*mie, 0);
  i_sub.Resize(2*n+3*iemax+ielpar+2+100, 0);


  // current number of constraints
  // remove inequality constraints from mie, add to mf. iemax

  // fixme!
  lact = m;
  setact.Resize(m+1, 1);
  // setact[0] = 1;
  // for(int i = 0; i < m; i++);
    // m_setact[i+1] = m_nie - m_nobj + i + 1; FIXME

}

/** Problem solve */
int FeasSCP::solve_problem(bool fromWarmstart)
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
  //assert(ieleng == ielpar || ieleng == ielpar-1); // the array has always min 1 entry to have a pointer

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
      LOG_DBG(feas_scp) << "sp: shift[" << i << "].shift -> " << shift[i].shift << " * " << g_scaling[i] << " = " << shift[i].shift * g_scaling[i];
      shift[i].shift *= g_scaling[i];
    }
  }

  // this confirms SCPIPBase based configuration against ScpSolver configuration
  // assert(ielpar >= n * 8); // m_ielpar = m_nvars*8; //+m_nvars*2*m_nvars + 20;
  // assert(iemax >= 3 * m);  // m_iemax = 3*m_nie;
  LOG_DBG2(feas_scp) << "sp: rdim=" << rdim << " n=" << n << " iemax=" << iemax << " ielpar=" << ielpar << " m=" << m << " mie=" << mie << " mf=" << mf;
  assert(rdim >= 44 * n + 16 * iemax + 2 * ielpar + 20 + 2 * mie + 100 * mie); // m_rdim = 44*m_nvars+16*m_iemax+2*m_ielpar+20+2*m_nie +100*m_nie
  assert(rsubdim >= 22 * n + 41 * iemax + 2 * ielpar + 28 + 100); // m_rsubdim = 22*m_nvars+41*m_iemax+2*m_ielpar+28 +100;
  assert(idim >= 7 * n + 8 * mie + 3 * ielpar + 12 + 100 * mie); // m_idim = 7*m_nvars+8*m_nie+3*m_ielpar+12+100*m_nie
  assert(isubdim >= 2 * n + 3 * iemax + ielpar + 2 + 100); // m_isubdim = 2*m_nvars+3*m_iemax+m_ielpar+2+100

  // double eps = GetLowerBoundsMVars(0);
  // double scale;
  for(;;) // we break out in success and error
  {
    // this might be dynamic!
    int spiwdim = spiw.GetSize();
    int spdwdim = spdw.GetSize();
    int ipopt_nargs = ipopt_args.GetSize();


    assert(iern.GetSize() == iecn.GetSize());
    assert(iern.GetSize() == iederv.GetSize());
    for(unsigned int i = 0; i < iern.GetSize(); i++)
      LOG_DBG3(feas_scp) << " sp: iern[" << i << "]=" << iern[i] << "\tiecn[" << i << "]=" << iecn[i] << "\tiederv[" << i << "]=" << iederv[i];

    LOG_DBG3(feas_scp) << "sp1: n = " << n; // 1
    LOG_DBG3(feas_scp) << "sp1: mie = " << mie; // 2
    LOG_DBG3(feas_scp) << "sp1: meq = " << meq; // 3
    LOG_DBG3(feas_scp) << "sp1: iemax = " << iemax; // 4
    LOG_DBG3(feas_scp) << "sp1: eqmax = " << eqmax; // 5
    LOG_DBG3(feas_scp) << "sp1: initial guess = " << x.ToString(); // 6
    LOG_DBG3(feas_scp) << "sp1: lower bounds = " << x_l.ToString(); // 7
    LOG_DBG3(feas_scp) << "sp1: upper bounds = " << x_u.ToString(); // 8
    LOG_DBG3(feas_scp) << "sp1: function value = " << f_org; // 9
    LOG_DBG3(feas_scp) << "sp1: inequality constraints = " << h_org.ToString(); // 10
    LOG_DBG3(feas_scp) << "sp1: equality constraints = " << g_org.ToString(); // 11
    LOG_DBG3(feas_scp) << "sp1: gradient = " << df.ToString(); // 12
    LOG_DBG3(feas_scp) << "sp1: lagrange multipliers (ineq.) = " << y_ie.ToString(); // 13
    LOG_DBG3(feas_scp) << "sp1: lagrange multipliers (eq.) = " << y_eq.ToString(); // 14
    LOG_DBG3(feas_scp) << "sp1: lagrange multipliers (l.b.) = " << y_l.ToString(); // 15
    LOG_DBG3(feas_scp) << "sp1: lagrange multipliers (u.b.) = " << y_u.ToString(); // 16
    LOG_DBG3(feas_scp) << "sp1: icntl = " << icntl.ToString(); // 17
    LOG_DBG3(feas_scp) << "sp1: rcntl = " << rcntl.ToString(); // 18
    LOG_DBG3(feas_scp) << "sp1: info = " << info.ToString(); // 19
    LOG_DBG3(feas_scp) << "sp1: rinfo = " << rinfo.ToString(); // 20
    LOG_DBG3(feas_scp) << "sp1: nout = " << nout; // 21
    LOG_DBG3(feas_scp) << "sp1: r_scp = " << r_scp.ToString(); // 22
    LOG_DBG3(feas_scp) << "sp1: r_dim = " << rdim; // 23
    LOG_DBG3(feas_scp) << "sp1: r_sub = " << r_sub.ToString(); // 24
    LOG_DBG3(feas_scp) << "sp1: r_subdim = " << rsubdim; // 25
    LOG_DBG3(feas_scp) << "sp1: i_scp = " << i_scp.ToString(); // 26
    LOG_DBG3(feas_scp) << "sp1: i_dim = " << idim; // 27
    LOG_DBG3(feas_scp) << "sp1: i_sub = " << i_sub.ToString(); // 28
    LOG_DBG3(feas_scp) << "sp1: i_subdim = " << isubdim; // 29
    LOG_DBG3(feas_scp) << "sp1: active = " << active.ToString(); // 30
    LOG_DBG3(feas_scp) << "sp1: mode = " << mode; // 31
    LOG_DBG3(feas_scp) << "sp1: ierr = " << ierr; // 32
    LOG_DBG3(feas_scp) << "sp1: iern = " << iern.ToString(); // 33
    LOG_DBG3(feas_scp) << "sp1: iecn = " << iecn.ToString(); // 34
    LOG_DBG3(feas_scp) << "sp1: iederv = " << iederv.ToString(); // 35
    LOG_DBG3(feas_scp) << "sp1: ielpar = " << ielpar; // 36
    LOG_DBG3(feas_scp) << "sp1: ieleng = " << ieleng; // 37
    LOG_DBG3(feas_scp) << "sp1: eqrn = " << eqrn.ToString(); // 38
    LOG_DBG3(feas_scp) << "sp1: eqcn = " << eqcn.ToString(); // 39
    LOG_DBG3(feas_scp) << "sp1: eqcoef = " << eqcoef.ToString(); // 40
    LOG_DBG3(feas_scp) << "sp1: eqlpar = " << eqlpar; // 41
    LOG_DBG3(feas_scp) << "sp1: eqleng = " << eqleng; // 42
    LOG_DBG3(feas_scp) << "sp1: mactiv = " << mactiv; // 43
    LOG_DBG3(feas_scp) << "sp1: spiw = " << spiw.ToString(); // 44
    LOG_DBG3(feas_scp) << "sp1: spiwdim = " << spiwdim; // 45
    LOG_DBG3(feas_scp) << "sp1: spdw = " << spdw.ToString(); // 46
    LOG_DBG3(feas_scp) << "sp1: spdwdim = " << spdwdim; // 47
    LOG_DBG3(feas_scp) << "sp1: spstrat = " << spstrat; // 48
    LOG_DBG3(feas_scp) << "sp1: linsys = " << linsys; // 49
    LOG_DBG3(feas_scp) << "sp1: linear = " << linear.ToString(); // 50
    LOG_DBG3(feas_scp) << "sp1: mf = " << mf; // 51
    LOG_DBG3(feas_scp) << "sp1: lact=" << lact; // 52
    LOG_DBG3(feas_scp) << "sp1: setact = " << setact.ToString(); // 53


    // int nEqMax = 1;
    scpip40i_(&n,                      // 1
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
              &linsys,                 // 49
              linear.GetPointer(),     // 50
              &mf,                     // 51
              &lact,                   // 52
              setact.GetPointer(),     // 53
              &ipopt_nargs,            // 54
              ipopt_args.GetPointer(), // 55
              ipopt_cargs);            // 56

    LOG_TRACE(feas_scp) << "scpip30 returns: ierr=" << ierr << " info[20-1]=" << info[20-1];

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
      LOG_DBG2(feas_scp) << "x: " << x.ToString();
      EvaluateFunctionValues();
    }
    if(ierr == -2)
    {
      LOG_DBG2(feas_scp) << "x: " << x.ToString();
      // TODO active constraints!
      if(!EvaluateGradients())
      {
        ierr = Gradients_Return_False;
        CallFinalizeSolution(); // allow comit of the last iteration to reuse it
        break; // no output for resart!
      }
      LOG_DBG2(feas_scp) << "df: " << df.ToString();
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






