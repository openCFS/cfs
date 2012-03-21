#include "Optimization/Optimizer/FeasSCP.hh"
#include "Optimization/Optimizer/SCPIPBase.hh"
#include "Optimization/Optimizer/scpip30.hh"

#include <time.h>


DECLARE_LOG(feas_scp)
DEFINE_LOG(feas_scp, "feas_scp")

using namespace CoupledField;



//=============================================================
//      AScpSolver implementation
//=============================================================

int FeasSCP::m_nMajor = 0;
int FeasSCP::m_nFEvals = 0;

static int default_ioptions[13] = {2, 1, 500, 3, 3, 0, 0, 1, 0, 0, 0, 0, 0};
static double default_doptions[10] = {1.0e-6, 1.0e30, 1.0e30, 1.0e-6, 1.0e-1, 1.0e-1, 0.e0, 1.0e30, 1.0e30, 1.0e30};

// constructor
FeasSCP::FeasSCP(Optimization* opt, PtrParamNode pn) :
      BaseOptimizer(opt, pn, Optimization::SNOPT_SOLVER),
    m_status(0),
    m_iemax(0),
    m_nvars(0),
    m_nobj(0),
    m_nie(0),
    m_neq(0),
    m_nsdp(0),
    m_xinit(0),
    m_blks(0),
    m_lbv(0),
    m_ubv(0),
    m_doptions(0),
    m_ioptions(0),
    m_f_org(0),
    m_h_org(0),
    m_g_org(0),
    m_df(0),
    m_y_ie(0),
    m_y_eq(0),
    m_y_l(0),
    m_y_u(0),
    m_ninfo(0),
    m_rinfo(0),
    m_nout(7),
    m_r_scp(0),
    m_rdim(0),
    m_r_sub(0),
    m_rsubdim(0),
    m_i_scp(0),
    m_idim(0),
    m_i_sub(0),
    m_isubdim(0),
    m_active(0),
    m_mactiv(0),
    m_mode(2),
    m_iern(0),
    m_iecn(0),
    m_iederv(0),
    m_ielpar(0),
    m_ieleng(0),
    m_eqrn(0),
    m_eqcn(0),
    m_eqcoef(0),
    m_eqlpar(1),
    m_eqleng(0),
    m_spiw(0),
    m_spiwdim(0),
    m_spdw(0),
    m_spdwdim(0),
    m_spstrat(1),
    m_linsys(1),
    m_linear(0),
    m_rScaleF(1.),
    m_rScaleE(100.),
    m_bAutoScale(true)
//-------------------------------------------------------------
// Create empty SCP problem structure.
//-------------------------------------------------------------
{
  int i=0, nObj = 1;

  m_nMajor = 1;
  m_nFEvals = 1;

  // set default options
  m_doptions = new double[10];
  for(i=0; i<7; i++)
    m_doptions[i] = default_doptions[i];

  m_ioptions = new int[13];
  for(i=0; i<13; i++)
    m_ioptions[i] = default_ioptions[i];

  ////////////////////// FMOK specific code /////////////////////////////////

  // nObj = pFMOProblem->Objective().NumObj();
  nObj = 1;

  // SetNumVars(pFMOProblem->NumVars() + 1);//((nObj > 1) ? 1 : 0)); //number of variables
  m_nvars = optimization->GetDesign()->GetNumberOfVariables(); // we have the slack variable in the problem
  m_nobj = nObj;
  // SetNumIneqConstr(pFMOProblem->NumConstr() + ((nObj > 1) ? nObj : 0));//number of constraints including linear constraints
  m_nie = optimization->constraints.view->GetNumberOfActiveConstraints();
  m_neq = 0;//number of constraints including linear constraints


  // temporary only
  StdVector<double> g_l(m_nie);
  StdVector<double> g_u(m_nie);

  GetBounds(m_nvars, x_l, x_u, m_nie, g_l, g_u);

  assert(x_l.GetSize() == m_nvars);
  assert(x_l[0] == -default_doptions[2]);
  assert(x_u[0] == default_doptions[2]);

  m_mf = optimization->GetDesign()->elements * 2; // only 2D
  // double rLowerEVBound = m_rScaleE*GetLowerBoundsMVars(i);

  m_iemax = m_nie + m_mf;

  //set up start value
  // SetInitialSolution(pFMOProblem->CurrentIterate()); // start value for E-matrices

  // Now, set some specific options (LATER READ FROM FMO Studio)

  // SetIntegerOption(9, 500);
  // SetRealOption(2, 1.0e-1);


#ifndef NOHF5
  if(szConfigFile) {

    // parse the strategy file and set up the problem accordingly
    ICfgDBCore* db = ICfgDBFactory::instance();
    // make it available via the global function interface for other libraries
    ICfgDBCInterface::init(db);
    F95_FMS_CFNAME(cdb_init)();

    int ierr = 0;

    if(db->load(szConfigFile)) {
      throw("Cannot open the config file!\n");
    }

    // Retrieve general options ...
    int nMaxIt = -1;
    string gNameP = string(fms::s_FMOPCfgParamGroupName) + "/" + fms::s_FMOPCfgMaxitName;
    db->getInt(gNameP, nMaxIt);
    SetIntegerOption(2, nMaxIt);
    // scp iterations

    double rPrecision = 1.0e-4;
    gNameP = string(fms::s_FMOPCfgParamGroupName) + "/" + fms::s_FMOPCfgToleranceName;
    db->getDouble(gNameP, rPrecision);
    //initial penalty

    //static const char s_Scp[] = "Scpf1.0";
    //static const char s_scpparams[] = "Algorithmic parameters";
    static const char s_Scp[] = "SCPF1.0";
    static const char s_scpparams[] = "Algorithmic parameters SCPF";
    static const char s_method[] = "method";
    static const char s_asymptotes[] = "asymptotes";
    static const char s_maxls[] = "maxls";
    static const char s_actset[] = "activeset";

    //static const char s_stopping[] = "Stopping criteria";
    static const char s_stopping[] = "Additional stopping criteria SCPF";
    static const char s_relaxed[] = "relaxed";
    static const char s_acc[] = "acc";
    static const char s_relx[] = "relx";
    static const char s_absf[] = "absf";
    static const char s_relf[] = "relf";

    static const char s_NLPIPB[] = "Algorithmic parameters NLPIPB";
    static const char s_accsub[] = "acc";
    static const char s_accqp[] = "accqp";
    static const char s_maxit[] = "maxit";
    static const char s_lssub[] = "lssub";
    static const char s_accactive[] = "accactive";

    static const char s_scaling[] = "Internal scaling";
    static const char s_scaletype[] = "scaletype";
    static const char s_scalef[] = "scalef";
    static const char s_scaleE[] = "scaleE";

    static const char s_DisplOptions[] = "Display options";
    static const char s_Output[] = "output";
    static const char s_OScheme[] = "oscheme";

    // Retrieve solver specific parameters ...
    string sTmp;
    gNameP = string(s_Scp) + "/" + string(s_scpparams) +  "/" + string(s_method);
    db->getString(gNameP, sTmp);


    SetIntegerOption(0, 8);


    if(!sTmp.compare("MMA")) {
      SetIntegerOption(0, 1);
    }
    else if(!sTmp.compare("SCP")) {
      SetIntegerOption(0, 2);
    }
    else if(!sTmp.compare("SCP-filter")) {
      SetIntegerOption(0, 3);
    }


    int nMaxLS = 0;
    gNameP = string(s_Scp) + "/" + string(s_scpparams) +  "/" + string(s_maxls);
    db->getInt(gNameP, nMaxLS);
    SetIntegerOption(4, nMaxLS);

    double rActive = 1.0e-4;
    gNameP = string(s_Scp) + "/" + string(s_scpparams) +  "/" + string(s_actset);
    db->getDouble(gNameP, rActive);
    SetRealOption(2, rActive);


    // Stopping
    bool bRelaxed = 1;
    gNameP = string(s_Scp) + "/" + string(s_stopping) +  "/" + string(s_relaxed);
    db->getBool(gNameP, bRelaxed);
    SetIntegerOption(5, (int) bRelaxed);

    double rRelX = 1.0e-4;
    gNameP = string(s_Scp) + "/" + string(s_stopping) +  "/" + string(s_relx);
    db->getDouble(gNameP, rRelX);
    SetRealOption(5, rRelX);

    double rAbsF = 1.0e-4;
    gNameP = string(s_Scp) + "/" + string(s_stopping) +  "/" + string(s_absf);
    db->getDouble(gNameP, rAbsF);
    SetRealOption(4, rAbsF);

    double rRelF = 1.0e-4;
    gNameP = string(s_Scp) + "/" + string(s_stopping) +  "/" + string(s_relf);
    db->getDouble(gNameP, rRelF);
    SetRealOption(3, rRelF);

    // Scaling
    bool bAutoScale = false;
    gNameP = string(s_Scp) + "/" + string(s_scaling) +  "/" + string(s_scaletype);
    db->getBool(gNameP, bAutoScale);
    m_bAutoScale = bAutoScale;
    if(bAutoScale){
      m_rScaleF = -1.1e0;
      m_rScaleE = 100;
    }else{
      double rScaleF = 1.0e1;
      gNameP = string(s_Scp) + "/" + string(s_scaling) +  "/" + string(s_scalef);
      db->getDouble(gNameP, rScaleF);
      m_rScaleF = rScaleF;

      double rScaleE = 1.0e1;
      gNameP = string(s_Scp) + "/" + string(s_scaling) +  "/" + string(s_scaleE);
      db->getDouble(gNameP, rScaleE);
      m_rScaleE = rScaleE;
    }

    // Output
    string sOutputLevel;
    gNameP = string(s_Scp) + "/" + string(s_DisplOptions) +  "/" + string(s_Output);
    db->getString(gNameP, sOutputLevel);
    if(!sOutputLevel.compare("iterations"))
      SetIntegerOption(3, 3);
    if(!sOutputLevel.compare("sub problems"))
      SetIntegerOption(3, 5);
    else if(!sOutputLevel.compare("verbose"))
      SetIntegerOption(3, 4);

    int nIScheme = -1;
    string gNameIS = string(fms::s_FMOPCfgParamGroupName) + "/" + fms::s_FMOPCfgIterSaveSchemeName;
    db->getInt(gNameIS, nIScheme);
    ICfgDBCore::EIterSaveScheme isType = static_cast<ICfgDBCore::EIterSaveScheme>(nIScheme);

    if(isType == ICfgDBCore::ITER_SAVE_SOLVER_SPECIFIC) {
      int nOScheme = 1;
      gNameP = string(s_Scp) + "/" + string(s_DisplOptions) +  "/" + string(s_OScheme);
      db->getInt(gNameP, nOScheme);
      m_nMajor = nOScheme;
    } else {
      m_nMajor = 1;
    }

    // with this the global functions will also be detached from the cfgDB
    // and the db instance will be deleted

    ICfgDBCInterface::deleteInstance();
  }
#endif



  ////////////////////// end FMOK specific code /////////////////////////////////

  m_ielpar = m_nvars*8; //+m_nvars*2*m_nvars + 20;
  m_iemax = 3*m_nie;

  // Allocate auxiliary structures
  m_h_org = new double[m_iemax];
  m_g_org = new double[1];
  m_df = new double[m_nvars];
  m_y_ie = new double[m_iemax];
  m_y_eq = new double[m_eqlpar];
  m_y_l = new double[m_nvars];
  m_y_u = new double[m_nvars];
  m_ninfo = new int[23];
  m_rinfo = new double[5];
  m_nout = 7;
  m_rdim = 44*m_nvars+16*m_iemax+2*m_ielpar+20+2*m_nie +100*m_nie; 
  m_r_scp = new double[m_rdim];

  m_rsubdim = 22*m_nvars+41*m_iemax+2*m_ielpar+28 +100;   //m_rsubdim = 22*m_nvars+41*m_nie+2*m_ielpar+28; 
  m_r_sub = new double[m_rsubdim];

  m_idim = 7*m_nvars+8*m_nie+3*m_ielpar+12+100*m_nie;
  m_i_scp = new int[m_idim];

  m_isubdim = 2*m_nvars+3*m_iemax+m_ielpar+2+100;   
  m_i_sub = new int[m_isubdim];
  m_active = new int[m_iemax];
  m_mode = 2;
  m_status = 0;
  m_iern = new int[m_ielpar];
  m_iecn = new int[m_ielpar];
  m_iederv = new double[m_ielpar];
  m_ieleng = m_ielpar;
  m_eqrn = new int[m_eqlpar];
  m_eqcn = new int[m_eqlpar];
  m_eqcoef = new double[m_eqlpar];
  m_eqleng = 0;
  m_spiwdim = 1; 
  m_spiw = new int[m_spiwdim];                 
  m_spdwdim = 1;    
  m_spdw = new double[m_spdwdim];  
  m_linear = new int[m_iemax];      

  for (i=0; i<m_iemax; i++){
    m_linear[i] = 1;
  }
  // dummy allocations ...



}

FeasSCP::~FeasSCP() {

  delete[] m_xinit;
  delete[] m_blks;
  delete[] m_doptions;
  delete[] m_ioptions;
  delete[] m_h_org;
  delete[] m_g_org;
  delete[] m_df;
  delete[] m_y_ie;
  delete[] m_y_l;
  delete[] m_y_u;
  delete[] m_y_eq;
  delete[] m_ninfo;
  delete[] m_rinfo;
  delete[] m_r_scp;
  delete[] m_r_sub;
  delete[] m_i_scp;
  delete[] m_i_sub;
  delete[] m_active;
  delete[] m_iern;
  delete[] m_iecn;
  delete[] m_eqrn;
  delete[] m_eqcn;
  delete[] m_eqcoef;
  delete[] m_iederv;
  delete[] m_spiw;                 
  delete[] m_spdw;      
  delete[] m_linear;
}

/** Problem solve */
int FeasSCP::SCPSolve(bool fromWarmstart)
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
              m_linear,                // 50
              &m_mf);                  // 51

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



