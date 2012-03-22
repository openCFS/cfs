#include "Optimization/Optimizer/FeasSCP.hh"
#include "Optimization/Optimizer/SCPIPBase.hh"
#include "Optimization/Optimizer/scpip40nlph.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"

#include <time.h>


DECLARE_LOG(feas_scp)
DEFINE_LOG(feas_scp, "feas_scp")

using namespace CoupledField;


// static int default_ioptions[13] = {2, 1, 500, 3, 3, 0, 0, 1, 0, 0, 0, 0, 0}
// static double default_doptions[10] = {1.0e-6, 1.0e30, 1.0e30, 1.0e-6, 1.0e-1, 1.0e-1, 0.e0, 1.0e30, 1.0e30, 1.0e30};

// constructor
FeasSCP::FeasSCP(Optimization* opt, PtrParamNode pn, Optimization::Optimizer type) : SCPIP(opt, pn, type)
{

}

FeasSCP::~FeasSCP()
{

}

void FeasSCP::AllocateProblem()
{
  SCPIPBase::AllocateProblem();

  linear.Resize(iemax);
  linear.Init(1);

  assert(false);
  mf = optimization->GetDesign()->elements * 2; // FIXME

  // fixme!
  lact = m;
  setact.Resize(m+1);
  setact[0] = 1;
  for(int i = 0; i < m; i++);
    // m_setact[i+1] = m_nie - m_nobj + i + 1; FIXME

}

/** Problem solve */
int FeasSCP::SolveProblem(bool fromWarmstart)
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
      LOG_DBG(feas_scp) << " shift[" << i << "].shift -> " << shift[i].shift << " * "
                          << g_scaling[i] << " = " << shift[i].shift * g_scaling[i];
      shift[i].shift *= g_scaling[i];
    }
  }

  // this confirms SCPIPBase based configuration against ScpSolver configuration
  assert(ielpar == n * 8); // m_ielpar = m_nvars*8; //+m_nvars*2*m_nvars + 20;
  assert(iemax == 3 * m);  // m_iemax = 3*m_nie;
  assert(rdim == 44 * n + 16 * iemax + 2 * ielpar + 20 + 2 * m + 100 * m); // m_rdim = 44*m_nvars+16*m_iemax+2*m_ielpar+20+2*m_nie +100*m_nie
  assert(rsubdim >= 22 * n + 41 * iemax + 2 * ielpar + 28 + 100); // m_rsubdim = 22*m_nvars+41*m_iemax+2*m_ielpar+28 +100;
  assert(idim >= 7 * n + 8 * m + 3 * ielpar + 12 + 100 * m); // m_idim = 7*m_nvars+8*m_nie+3*m_ielpar+12+100*m_nie
  assert(isubdim >= 2 * n + 3 * iemax + ielpar + 2 + 100); // m_isubdim = 2*m_nvars+3*m_iemax+m_ielpar+2+100

  // double eps = GetLowerBoundsMVars(0);
  // double scale;
  for(;;) // we break out in success and error
  {
    // this might be dynamic!
    int spiwdim = spiw.GetSize();
    int spdwdim = spdw.GetSize();

    LOG_DBG3(feas_scp) << "before call to scpip: n = " << n << ", mie = " << mie << ", meq = " << meq
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
              setact.GetPointer());    // 53

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
    if(ierr == -17)
    {
      /*
      int flag = 0;
      int m_elements = GetNumSDPBlocks();
      int nTmp = m_nie - m_mf;
      int nMSize = m_blks[2];
      concavity_d_3x3_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);

       */
      LOG_DBG2(feas_scp) << "! COMPUTATION OF CONSTRAINTS MIE-MF,... MIE Missing!!";
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
    if((ierr < -5 && ierr != -17) || ierr > 0)
    {
      std::cerr << ToString(ierr) << std::endl;
      CallFinalizeSolution();
      break;
    }
  }

  return ierr;
}





