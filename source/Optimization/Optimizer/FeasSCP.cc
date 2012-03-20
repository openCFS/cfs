//-------------------------------------------------------------
//       Author:
//               Michael Stingl, 14. Mar, 2008.
//-------------------------------------------------------------

#include "ScpSolver.h"
#include "..\..\fmojob.h"
#include <time.h>

//=============================================================
//      AScpSolver implementation
//=============================================================

void* ASCPSolver::m_pUsrObj;
int ASCPSolver::m_nMajor = 0;
int ASCPSolver::m_nFEvals = 0;

static int default_ioptions[13] = {2, 1, 500, 3, 3, 0, 0, 1, 0, 0, 0, 0, 0};
static double default_doptions[10] = {1.0e-6, 1.0e30, 1.0e30, 1.0e-6, 1.0e-1, 1.0e-1, 0.e0, 1.0e30, 1.0e30, 1.0e30};

// constructor
ASCPSolver::ASCPSolver(void * pUsrObj, char* szConfigFile) : 
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
m_lbmv(0),
m_ubmv(0),
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
m_mf(0),
m_lact(0),
m_setact(0),
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

  m_pUsrObj = pUsrObj;
  AFMOProblem* pFMOProblem = (AFMOProblem*) UsrObj();

  ////////////////////// FMOK specific code /////////////////////////////////

  nObj = pFMOProblem->Objective().NumObj();

  SetNumVars(pFMOProblem->NumVars() + 1);//((nObj > 1) ? 1 : 0)); //number of variables
  SetNumObjs(nObj); //number of variables
  SetNumIneqConstr(pFMOProblem->NumConstr() + ((nObj > 1) ? nObj : 0));//number of constraints including linear constraints
  SetNumEqConstr(0);//number of constraints including linear constraints

  SetNumSDPBlocks(pFMOProblem->NumMVars()); //number of matrix variables 
  SetSizesOfMVariables(pFMOProblem->BlockSizes()()); //size of matrices

  // Set lower and upper bounds for matrix variables
  SetLowerBoundsMVars(pFMOProblem->LowerMBounds()());  // on matrix variables
  SetUpperBoundsMVars(pFMOProblem->UpperMBounds()());

  // Set lower and upper bounds for scalar variables
  SetLowerBoundVars(-default_doptions[2]);  // on scalar variables
  SetUpperBoundVars(default_doptions[2]);

  int nMVars = pFMOProblem->NumMVars();
  int nVOffset = 0;
  m_mf = 0;

  for(i=0; i<nMVars; i++) {
    int nMSize = m_blks[i];
    //double rLowerEVBound = GetLowerBoundsMVars(i);
    double rLowerEVBound = m_rScaleE*GetLowerBoundsMVars(i);
    if(nMSize == 1) {
      SetLowerBoundVar(rLowerEVBound, nVOffset);
      nVOffset += 1;
    } else if(nMSize == 2) {
      SetLowerBoundVar(rLowerEVBound, nVOffset);
      SetLowerBoundVar(rLowerEVBound, nVOffset+2);
	  //m_nie +=  1;
	  m_mf +=  1;
      nVOffset += 3;
    } else if(nMSize == 3) {
      SetLowerBoundVar(rLowerEVBound, nVOffset);
      SetLowerBoundVar(rLowerEVBound, nVOffset+2);
      SetLowerBoundVar(rLowerEVBound, nVOffset+5);
	  //m_nie +=  1;
	  m_mf +=  2;
      nVOffset += 6;
    } else if(nMSize == 6) {
      SetLowerBoundVar(rLowerEVBound, nVOffset);
      SetLowerBoundVar(rLowerEVBound, nVOffset+2);
      SetLowerBoundVar(rLowerEVBound, nVOffset+5);
      SetLowerBoundVar(rLowerEVBound, nVOffset+9);
      SetLowerBoundVar(rLowerEVBound, nVOffset+14);
      SetLowerBoundVar(rLowerEVBound, nVOffset+20);
	  //m_nie += 1;
	  m_mf += 5;
      nVOffset += 21;
    } 
  }
  m_iemax = m_nie + m_mf;

  //set up start value
  SetInitialSolution(pFMOProblem->CurrentIterate()); // start value for E-matrices

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
	SetRealOption(0, rPrecision);
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
    
 //   gNameP = string(s_Scp) + "/" + string(s_scpparams) +  "/" + string(s_asymptotes);
 //   db->getString(gNameP, sTmp);

	//if(!sTmp.compare("Strategy 1")) {
 //     SetIntegerOption(1, 1);
	//}
	//else if(!sTmp.compare("Strategy 2")) {
 //     SetIntegerOption(1, 2);
	//}
	//else if(!sTmp.compare("Strategy 3")) {
 //     SetIntegerOption(1, 3);
	//}
	//else if(!sTmp.compare("Strategy 4")) {
 //     SetIntegerOption(1, 3);
	//}
    
    int nMaxLS = 0;
    gNameP = string(s_Scp) + "/" + string(s_scpparams) +  "/" + string(s_maxls);
    db->getInt(gNameP, nMaxLS);
    SetIntegerOption(4, nMaxLS);
    
    double rActive = 1.0e-4;
    gNameP = string(s_Scp) + "/" + string(s_scpparams) +  "/" + string(s_actset);
    db->getDouble(gNameP, rActive);
	SetRealOption(2, rActive);

 //   double rKKT = 1.0e-4;
 //   gNameP = string(s_Scp) + "/" + string(s_stopping) +  "/" + string(s_acc);
 //   db->getDouble(gNameP, rKKT);
	//SetRealOption(0, rKKT);

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

//NLPIPB
	double racc = 1.0e-4;
	gNameP = string(s_Scp) + "/" + string(s_NLPIPB) +  "/" + string(s_acc);
    db->getDouble(gNameP, racc);
	SetRealOption(7, racc);

	double raccqp = 1.0e-4;
	gNameP = string(s_Scp) + "/" + string(s_NLPIPB) +  "/" + string(s_accqp);
    db->getDouble(gNameP, raccqp);
	SetRealOption(8, raccqp);

	int imaxit = 100;
	gNameP = string(s_Scp) + "/" + string(s_NLPIPB) +  "/" + string(s_maxit);
    db->getInt(gNameP, imaxit);
	SetIntegerOption(8, imaxit);

	int ilssub = 5;
	gNameP = string(s_Scp) + "/" + string(s_NLPIPB) +  "/" + string(s_lssub);
    db->getInt(gNameP, ilssub);
	SetIntegerOption(9, ilssub);

	double ractivesub = 1.0e-6;
	gNameP = string(s_Scp) + "/" + string(s_NLPIPB) +  "/" + string(s_accactive);
    db->getDouble(gNameP, ractivesub);
	SetRealOption(9, ractivesub);


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

  //m_ielpar = m_nvars+m_nvars*2*m_nvars + 20;  // 10 is fixed for number of load cases
  m_ielpar = m_nvars*m_nie+m_mf*9; //+m_nvars*2*m_nvars + 20;
  //m_ielpar = 3*m_nvars+m_nvars*10; 
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
  //m_rdim = (22*m_nvars+41*m_iemax+2*m_ielpar+28)+(((28*m_nvars+11*m_iemax+m_ielpar +18)*2+5*(m_iemax+m_nvars))*10)*5; 
  m_rdim = 44*m_nvars+16*m_iemax+2*m_ielpar+20+2*m_nie +100*m_nie; 
  m_r_scp = new double[m_rdim];

  m_rsubdim = 22*m_nvars+41*m_iemax+2*m_ielpar+28 +100;   //m_rsubdim = 22*m_nvars+41*m_nie+2*m_ielpar+28; 
  m_r_sub = new double[m_rsubdim];

  //m_idim = ((7*m_nvars+8*m_iemax+m_ielpar+10)+2*m_nvars+3*m_iemax+m_ielpar+2+2*(m_iemax+m_nvars))*100;
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
  //m_spiwdim = 2*m_nvars+5*m_iemax+6; 
  m_spiwdim = 1; 
  //m_spiwdim = 2*m_nvars+5*m_nie+6; 
  m_spiw = new int[m_spiwdim];                 
  //m_spdwdim = m_iemax*m_iemax;      
  m_spdwdim = 1;    
  m_spdw = new double[m_spdwdim];  
  m_linear = new int[m_iemax];    
  m_setact = new int[m_nie];

  for (i=0; i<m_iemax; i++){
    m_linear[i] = 1;
  }

  // PARAMS.DAT
  ofstream myfile;
  myfile.open("PARAMS.DAT");

	myfile<<"IMAXITER 200";
    myfile<<"\n";
	myfile<<"DTOL 1.D-7";
    myfile<<"\n";
	myfile<<"DCMAXTOL 1.D-6";
    myfile<<"\n";
	myfile<<"IMERIT 2";
    myfile<<"\n";
	myfile<<"IPRINT 2";
    myfile<<"\n";
	myfile<<"IFILE 1";
    myfile<<"\n";
	myfile<<"ISCALE 2";
    myfile<<"\n";
	myfile<<"IQUASI 0";
    myfile<<"\n";
	myfile<<"ITRON2DERIVS 0";
    myfile<<"\n";
  myfile.close();

}

ASCPSolver::~ASCPSolver() { 

  delete[] m_xinit;
	//cout<<"scpsolver";
  delete[] m_blks;
  delete[] m_lbmv;
  delete[] m_ubmv;
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
  delete[] m_setact;
}

void 
ASCPSolver::SetInitialSolution(double* xinit) {
  /* requires call to SetNumVars first! */
  if (m_nvars == 0)
    return;

  if(m_xinit)
    delete[] m_xinit;
//	sonja - geaendert auf m_nvars+1
  m_xinit = new double[m_nvars];
  for (int i=0; i<m_nvars; i++)
    m_xinit[i] = xinit[i];
}

void 
ASCPSolver::SetLowerBoundsMVars(const double* rLBMVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0)
    return;

  if(m_lbmv)
    delete[] m_lbmv;

  m_lbmv = new double[m_nsdp];

  for (int i=0; i<m_nsdp; i++)
    m_lbmv[i] = rLBMVars[i];   
}

void 
ASCPSolver::SetUpperBoundsMVars(const double* rUBMVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0)
    return;

  if(m_ubmv)
    delete[] m_ubmv;

  m_ubmv = new double[m_nsdp];

  for (int i=0; i<m_nsdp; i++)
    m_ubmv[i] = rUBMVars[i];   
}

void 
ASCPSolver::SetLowerBoundMVars(double rLBMVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0)
    return;

  if(m_lbmv)
    delete m_lbmv;

  m_lbmv = new double[m_nsdp];

  for (int i=0; i<m_nsdp; i++)
    m_lbmv[i] = rLBMVars;   
}

void 
ASCPSolver::SetUpperBoundMVars(double rUBMVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0)
    return;

  if(m_ubmv)
    delete m_ubmv;

  m_ubmv = new double[m_nsdp];

  for (int i=0; i<m_nsdp; i++)
    m_ubmv[i] = rUBMVars;   
}

void 
ASCPSolver::SetLowerBoundMVar(double rLBMVar, int idx) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0 || idx < 0 || idx >= m_nsdp)
    return;

  m_lbmv[idx] = rLBMVar;   
}

void 
ASCPSolver::SetUpperBoundMVar(double rUBMVar, int idx) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0 || idx < 0 || idx >= m_nsdp)
    return;

  m_ubmv[idx] = rUBMVar;   
}

void 
ASCPSolver::SetLowerBoundsVars(const double* rLBVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nvars == 0)
    return;

  if(m_lbv)
    delete[] m_lbv;

  m_lbv = new double[m_nvars];

  for (int i=0; i<m_nvars; i++)
    m_lbv[i] = rLBVars[i];   
}

void 
ASCPSolver::SetUpperBoundsVars(const double* rUBVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nvars == 0)
    return;

  if(m_ubv)
    delete[] m_ubv;

  m_ubv = new double[m_nvars];

  for (int i=0; i<m_nvars; i++)
    m_ubv[i] = rUBVars[i];   
}

void 
ASCPSolver::SetLowerBoundVars(double rLBVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nvars == 0)
    return;

  if(m_lbv)
    delete m_lbv;

  m_lbv = new double[m_nvars];

  for (int i=0; i<m_nvars; i++)
    m_lbv[i] = rLBVars;   
}

void 
ASCPSolver::SetUpperBoundVars(double rUBVars) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nvars == 0)
    return;

  if(m_ubv)
    delete m_ubv;

  m_ubv = new double[m_nvars];

  for (int i=0; i<m_nvars; i++)
    m_ubv[i] = rUBVars;   
}

void 
ASCPSolver::SetLowerBoundVar(double rLBVar, int idx) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nvars == 0 || idx < 0 || idx >= m_nvars)
    return;

  m_lbv[idx] = rLBVar;   
}

void 
ASCPSolver::SetUpperBoundVar(double rUBVar, int idx) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nvars == 0 || idx < 0 || idx >= m_nvars)
    return;

  m_ubv[idx] = rUBVar;   
}

void 
ASCPSolver::SetSizesOfMVariables(const int* nMSizes) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0)
    return;

  if(m_blks)
    delete[] m_blks;

  m_blks = new int[m_nsdp];

  for (int i=0; i<m_nsdp; i++)
    m_blks[i] = nMSizes[i];   
}

void 
ASCPSolver::SetSizeOfMVariables(int nMSizes) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0)
    return;

  if(m_blks)
    delete m_blks;

  m_blks = new int[m_nsdp];

  for (int i=0; i<m_nsdp; i++)
    m_blks[i] = nMSizes;   
}

void 
ASCPSolver::SetSizeOfMVariable(int nMSize, int idx) {
  /* requires call to SetNumSDPBlocks first! */
  if (m_nsdp == 0 || idx < 0 || idx >= m_nsdp)
    return;

  m_blks[idx] = nMSize;   
}

void 
ASCPSolver::SetIntegerOptions(int* nVal) {
  for (int i=0; i<13; i++)
    m_ioptions[i] = nVal[i];
}

void 
ASCPSolver::SetIntegerOption(int idx, int nVal) {
  if(idx < 0 || idx > 11)
    return;

  m_ioptions[idx] = nVal;
}

void 
ASCPSolver::SetRealOptions(double* rVal) {
  for (int i=0; i<10; i++)
    m_doptions[i] = rVal[i];
}

void 
ASCPSolver::SetRealOption(int idx, double rVal) {
  if(idx < 0 || idx > 11)
    return;

  m_doptions[idx] = rVal;
}

/** Problem solve */
bool
ASCPSolver::SCPSolve() {

  double eps = GetLowerBoundsMVars(0);
  double scale;

  if(!m_nvars || !m_nobj) 
    return false;

  m_status = 0;


  int j = m_nvars/6;	
  for (int i=0; i<j; i++){
	  m_xinit[(i+1)*6-6] *= m_rScaleE;
	  m_xinit[(i+1)*6-4] *= m_rScaleE;
	  m_xinit[(i+1)*6-1] *= m_rScaleE;
  }
  if (m_nobj == 1) {	
    m_f_org = ((AFMOProblem*) UsrObj())->Objective().Eval(m_xinit, true, 0);
	if (m_rScaleF < 0){
		scale = 100;
		cout << scale<< endl;
	}else{
		scale = m_rScaleF;
		cout << scale<< endl;
	}


  
//SONJA
	m_nie = m_nie+1;
	m_xinit[m_nvars-1]=1.10;
    m_lbv[m_nvars-1]= 0.0;   
    m_ubv[m_nvars-1]=1.e30;

    m_lact = 1+ m_nobj;
	m_setact[0] = 1;
	for (int i=0; i<m_nobj; i++){
		m_setact[i+1] = m_nie - m_nobj + i + 1;
	}

  }else{
    m_lact = 1+ m_nobj;
	m_setact[0] = 1;
	for (int i=0; i<m_nobj; i++){
		m_setact[i+1] = m_nie - m_nobj + i + 1;
	}
	
  }

# if defined(_WIN32)
  do {
    int nEqMax = 1;

	SCPIP40I(&m_nvars,       
      &m_nie,         
      &m_neq,   
      &m_iemax,       
      &nEqMax,        
      m_xinit,             
      m_lbv,   
      m_ubv,
      &m_f_org, 
      m_h_org,    
      m_g_org,       
      m_df,  
      m_y_ie,             
      m_y_eq,
      m_y_l,
      m_y_u,
      m_ioptions,    
      m_doptions,    
      m_ninfo,
      m_rinfo,   
      &m_nout,
      m_r_scp, 
      &m_rdim, 
      m_r_sub,              
      &m_rsubdim,    
      m_i_scp, 
      &m_idim, 
      m_i_sub,              
      &m_isubdim,    
      m_active,
      &m_mode,
      &m_status,
      m_iern,                
      m_iecn,
      m_iederv,
      &m_ielpar,
      &m_ieleng,
      m_eqrn,                
      m_eqcn,
      m_eqcoef,
      &m_eqlpar,
      &m_eqleng,
      &m_mactiv,
      m_spiw,                    
      &m_spiwdim,
      m_spdw,                 
      &m_spdwdim,                      
      &m_spstrat,
      &m_linsys,
      m_linear,
      &m_mf,
	 &m_lact,
	 m_setact);

	if (m_ioptions[9] == 30){
		if (m_status == -1) {
			scale = 1;
			m_xinit[m_nvars-1] = 1;
		}
	}
	//if (m_ioptions[9] == 50){
	//	if (m_status == -1) {
	//		m_doptions[2] = 10;
	//	}
	//}

    if (m_status == -1) {					// Function evaluation

	  m_nFEvals++;

	  if(IsMajor())
		AFMOJob::TheInstance()->WriteSolution(m_xinit);
      

	  if(!((AFMOProblem*) UsrObj())->XKnown(m_xinit)) { 
		return false;
	  }

      if (m_nobj > 1) {						//Multiple load cases
        m_f_org = m_xinit[m_nvars-1];

        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++){ // Constraints
          m_h_org[i] = ((AFMOProblem*) UsrObj())->Constraint(i).Eval(m_xinit, false);
        }
        int idx = ((AFMOProblem*) UsrObj())->NumConstr();
        for (int i=0; i<((AFMOProblem*) UsrObj())->Objective().NumObj(); i++){ // Objective(s)
          m_h_org[idx+i] = ((AFMOProblem*) UsrObj())->Objective().Eval(m_xinit, false, i);
		  m_h_org[idx+i] = m_h_org[idx+i]*scale;
        }
        int flag = 1;
        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie;
		int nMSize = m_blks[2];				//at the moment all blocks have the same size
		double eps = 0.003333*m_rScaleE;
		if (nMSize == 2){
			CONCAVITY_D_2X2(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 3){
			CONCAVITY_D_3X3(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 6){
			CONCAVITY_D_6X6(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} 
//sonja
		double vol;
		vol = m_rScaleE*m_elements/3-m_elements/3;
		m_h_org[0] = m_h_org[0] - vol;

      } else {								//Single load case         
        m_f_org = ((AFMOProblem*) UsrObj())->Objective().Eval(m_xinit, false, 0);
		m_f_org = m_f_org*scale;
        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++){// Constraints
          m_h_org[i] = ((AFMOProblem*) UsrObj())->Constraint(i).Eval(m_xinit, false);  
        }

		m_h_org[m_nie-1] = m_f_org-m_xinit[m_nvars-1];
		m_f_org = m_xinit[m_nvars-1];

		if (m_h_org[m_nie-1] > 0){			//strict feasibility for the compliance function
			m_xinit[m_nvars-1]=m_xinit[m_nvars-1]+m_h_org[m_nie-1]+0.1; //;
			m_h_org[m_nie-1]=0-0.1; //;
			m_f_org = m_xinit[m_nvars-1];
		}

        int flag = 1;
        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie;
		int nMSize = m_blks[2];
		double eps = 0.003333*m_rScaleE;

		if (nMSize == 2){
			CONCAVITY_D_2X2(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 3){
			CONCAVITY_D_3X3(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 6){
			CONCAVITY_D_6X6(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} 

		double vol;
		vol = m_rScaleE*m_elements/3-m_elements/3;
		m_h_org[0] = m_h_org[0] - vol;

      }


    } else if (m_status == -2) {             // Gradient evaluation
      if (m_nobj > 1) {						 //Multiple load cases
        m_ieleng = 0;     
        m_mactiv = 0;    

        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++) { // Constraints
          if (m_active[i] == 1){									   			
            AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Constraint(i).Grad(m_xinit, false, 0);
            m_mactiv++;
            for(int j=0; j<pGrad->NNZ(); j++) {
              double val = pGrad->Val(j);
              int ind = pGrad->Ind(j)+1;		// C++ starts with 0, FORTRAN with 1
              m_iederv[m_ieleng] = val;
              m_iern[m_ieleng] = ind;
              m_iecn[m_ieleng] = m_mactiv;
              m_ieleng = m_ieleng + 1;
            }
          }
        }
        for (int i=0; i<((AFMOProblem*) UsrObj())->Objective().NumObj(); i++) { // Objective(s)
          if (m_active[i] == 1){	
            m_mactiv++;
            AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Objective().Grad(m_xinit, false, 0, i);
            for(int j=0; j<pGrad->NNZ(); j++) {
              double val = pGrad->Val(j);
              int ind = pGrad->Ind(j)+1;		
              m_iederv[m_ieleng] = val*scale;
              m_iern[m_ieleng] = ind;
              m_iecn[m_ieleng] = m_mactiv;
              m_ieleng = m_ieleng + 1;
            }
          }
        }

        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie;
		int nMSize = m_blks[2];
		double eps = 0.003333*m_rScaleE;
		if (nMSize == 2){
		    CONCAVITY_GRAD_2X2(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 3){
			CONCAVITY_GRAD_3X3(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 6){
			CONCAVITY_GRAD_6X6(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		}
        for (int i=0; i < (m_nvars-1); i++){
          m_df[i] = 0.0;
        }
        m_df[m_nvars-1] = 1.0;

      } else {									//Single load case
        m_ieleng = 0;     
        m_mactiv = 0;   
        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++) { 
          if (m_active[i] == 1){									  
            AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Constraint(i).Grad(m_xinit, false, 0);
            m_mactiv++;
            for(int j=0; j<pGrad->NNZ(); j++) {
              double val = pGrad->Val(j);
              int ind = pGrad->Ind(j)+1;	//hier gehoert das +1 hin, weil ja als row number gespeichert wird!
              m_iederv[m_ieleng] = val;
              m_iern[m_ieleng] = ind;
              m_iecn[m_ieleng] = m_mactiv;
              m_ieleng++;
            }
          }
        }
        AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Objective().Grad(m_xinit, false, 0, 1);

		m_mactiv++;
        for(int j=0; j<pGrad->NNZ(); j++) {
          double val = pGrad->Val(j);
          int ind = pGrad->Ind(j)+1;  
          m_df[ind] = val*scale;

          m_iederv[m_ieleng] = val*scale;
          m_iern[m_ieleng] = ind;
          m_iecn[m_ieleng] = m_mactiv;
          m_ieleng++;
        }
        m_iederv[m_ieleng] = -1.0;
        m_iern[m_ieleng] = m_nvars;
        m_iecn[m_ieleng] = m_mactiv;
        m_ieleng++;

		for (int i=0; i < (m_nvars-1); i++){
          m_df[i] = 0.0;
        }
        m_df[m_nvars-1] = 1.0;



        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie - 1;
		int nMSize = m_blks[2];
		double eps = 0.003333*m_rScaleE;
		if (nMSize == 2){
			CONCAVITY_GRAD_2X2(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 3){
			CONCAVITY_GRAD_3X3(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 6){
			CONCAVITY_GRAD_6X6(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		}        


      }
    } else if (m_status == -17) {			// Stepsize Reduction!
      int flag = 0;
      int m_elements = GetNumSDPBlocks();
      int nTmp = m_nie - m_mf;
	  int nMSize = m_blks[2];
	  if (nMSize == 2){
		 CONCAVITY_D_2X2(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
	  } else if (nMSize == 3){
		 CONCAVITY_D_3X3(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
	  } else if (nMSize == 6){
		 CONCAVITY_D_6X6(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
	  } 

    } else if (m_status > 1) {					 // error case!!!!     
	  cout<<"****   ERROR   **** IERR = "<< m_status;
	  cout<<"\n";
      return 1;  // TO BE REDIFINED!
    } else if (m_status == 1) {					 // error case!!!!     
	  cout<<"****   ERROR   **** IERR = "<< m_status<<endl;
	  cout<<"\n";
      m_status = 0;
    }


  } while (m_status < 0);


# else


  do {
    int nEqMax = 1;

	scpip40i_(&m_nvars,       
      &m_nie,         
      &m_neq,   
      &m_iemax,       
      &nEqMax,        
      m_xinit,             
      m_lbv,   
      m_ubv,
      &m_f_org, 
      m_h_org,    
      m_g_org,       
      m_df,  
      m_y_ie,             
      m_y_eq,
      m_y_l,
      m_y_u,
      m_ioptions,    
      m_doptions,    
      m_ninfo,
      m_rinfo,   
      &m_nout,
      m_r_scp, 
      &m_rdim, 
      m_r_sub,              
      &m_rsubdim,    
      m_i_scp, 
      &m_idim, 
      m_i_sub,              
      &m_isubdim,    
      m_active,
      &m_mode,
      &m_status,
      m_iern,                
      m_iecn,
      m_iederv,
      &m_ielpar,
      &m_ieleng,
      m_eqrn,                
      m_eqcn,
      m_eqcoef,
      &m_eqlpar,
      &m_eqleng,
      &m_mactiv,
      m_spiw,                    
      &m_spiwdim,
      m_spdw,                 
      &m_spdwdim,                      
      &m_spstrat,
      &m_linsys,
      m_linear,
      &m_mf,
	 &m_lact,
	 m_setact);

    if (m_status == -1) {					// Function evaluation

	  m_nFEvals++;

	  if(IsMajor())
		AFMOJob::TheInstance()->WriteSolution(m_xinit);
      

	  if(!((AFMOProblem*) UsrObj())->XKnown(m_xinit)) { 
		return false;
	  }

      if (m_nobj > 1) {						//Multiple load cases
        m_f_org = m_xinit[m_nvars-1];

        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++){ // Constraints
          m_h_org[i] = ((AFMOProblem*) UsrObj())->Constraint(i).Eval(m_xinit, false);
        }
        int idx = ((AFMOProblem*) UsrObj())->NumConstr();
        for (int i=0; i<((AFMOProblem*) UsrObj())->Objective().NumObj(); i++){ // Objective(s)
          m_h_org[idx+i] = ((AFMOProblem*) UsrObj())->Objective().Eval(m_xinit, false, i);
		  m_h_org[idx+i] = m_h_org[idx+i]*scale;
        }
        int flag = 1;
        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie;
		int nMSize = m_blks[2];				//at the moment all blocks have the same size
		double eps = 0.003333*m_rScaleE;
		if (nMSize == 2){
			concavity_d_2x2_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 3){
			concavity_d_3x3_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 6){
			concavity_d_6x6_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} 
//sonja
		double vol;
		vol = m_rScaleE*m_elements/3-m_elements/3;
		m_h_org[0] = m_h_org[0] - vol;

      } else {								//Single load case         
        m_f_org = ((AFMOProblem*) UsrObj())->Objective().Eval(m_xinit, false, 0);
		m_f_org = m_f_org*scale;
        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++){// Constraints
          m_h_org[i] = ((AFMOProblem*) UsrObj())->Constraint(i).Eval(m_xinit, false);  
        }

		m_h_org[m_nie-1] = m_f_org-m_xinit[m_nvars-1];
		m_f_org = m_xinit[m_nvars-1];

		if (m_h_org[m_nie-1] > 0){			//strict feasibility for the compliance function
			m_xinit[m_nvars-1]=m_xinit[m_nvars-1]+m_h_org[m_nie-1]+0.1; //;
			m_h_org[m_nie-1]=0-0.1; //;
			m_f_org = m_xinit[m_nvars-1];
		}

        int flag = 1;
        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie;
		int nMSize = m_blks[2];
		double eps = 0.003333*m_rScaleE;

		if (nMSize == 2){
			concavity_d_2x2_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 3){
			concavity_d_3x3_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} else if (nMSize == 6){
			concavity_d_6x6_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
		} 

		double vol;
		vol = m_rScaleE*m_elements/3-m_elements/3;
		m_h_org[0] = m_h_org[0] - vol;
      }


    } else if (m_status == -2) {             // Gradient evaluation
      if (m_nobj > 1) {						 //Multiple load cases
        m_ieleng = 0;     
        m_mactiv = 0;    

        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++) { // Constraints
          if (m_active[i] == 1){									   			
            AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Constraint(i).Grad(m_xinit, false, 0);
            m_mactiv++;
            for(int j=0; j<pGrad->NNZ(); j++) {
              double val = pGrad->Val(j);
              int ind = pGrad->Ind(j)+1;		// C++ starts with 0, FORTRAN with 1
              m_iederv[m_ieleng] = val;
              m_iern[m_ieleng] = ind;
              m_iecn[m_ieleng] = m_mactiv;
              m_ieleng = m_ieleng + 1;
            }
          }
        }
        for (int i=0; i<((AFMOProblem*) UsrObj())->Objective().NumObj(); i++) { // Objective(s)
          if (m_active[i] == 1){	
            m_mactiv++;
            AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Objective().Grad(m_xinit, false, 0, i);
            for(int j=0; j<pGrad->NNZ(); j++) {
              double val = pGrad->Val(j);
              int ind = pGrad->Ind(j)+1;		
              m_iederv[m_ieleng] = val*scale;
              m_iern[m_ieleng] = ind;
              m_iecn[m_ieleng] = m_mactiv;
              m_ieleng = m_ieleng + 1;
            }
          }
        }

        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie;
		int nMSize = m_blks[2];
		double eps = 0.003333*m_rScaleE;
		if (nMSize == 2){
		    concavity_grad_2x2_(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 3){
			concavity_grad_3x3_(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 6){
			concavity_grad_6x6_(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		}
        for (int i=0; i < (m_nvars-1); i++){
          m_df[i] = 0.0;
        }
        m_df[m_nvars-1] = 1.0;

      } else {									//Single load case
        m_ieleng = 0;     
        m_mactiv = 0;   
        for (int i=0; i<((AFMOProblem*) UsrObj())->NumConstr(); i++) { 
          if (m_active[i] == 1){									  
            AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Constraint(i).Grad(m_xinit, false, 0);
            m_mactiv++;
            for(int j=0; j<pGrad->NNZ(); j++) {
              double val = pGrad->Val(j);
              int ind = pGrad->Ind(j)+1;	//hier gehoert das +1 hin, weil ja als row number gespeichert wird!
              m_iederv[m_ieleng] = val;
              m_iern[m_ieleng] = ind;
              m_iecn[m_ieleng] = m_mactiv;
              m_ieleng++;
            }
          }
        }
        AFMOKSparseGradient *pGrad = ((AFMOProblem*) UsrObj())->Objective().Grad(m_xinit, false, 0, 1);

		m_mactiv++;
        for(int j=0; j<pGrad->NNZ(); j++) {
          double val = pGrad->Val(j);
          int ind = pGrad->Ind(j)+1;  
          m_df[ind] = val*scale;

          m_iederv[m_ieleng] = val*scale;
          m_iern[m_ieleng] = ind;
          m_iecn[m_ieleng] = m_mactiv;
          m_ieleng++;
        }
        m_iederv[m_ieleng] = -1.0;
        m_iern[m_ieleng] = m_nvars;
        m_iecn[m_ieleng] = m_mactiv;
        m_ieleng++;

		for (int i=0; i < (m_nvars-1); i++){
          m_df[i] = 0.0;
        }
        m_df[m_nvars-1] = 1.0;



        int m_elements = GetNumSDPBlocks();
        int nTmp = m_nie - 1;
		int nMSize = m_blks[2];
		double eps = 0.003333*m_rScaleE;
		if (nMSize == 2){
			concavity_grad_2x2_(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 3){
			concavity_grad_3x3_(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		} else if (nMSize == 6){
			concavity_grad_6x6_(m_iederv,m_iern,m_iecn,&m_ieleng,m_active,&nTmp,&m_mactiv,m_xinit,&m_elements,&eps);
		}        


      }
    } else if (m_status == -17) {			// Stepsize Reduction!
      int flag = 0;
      int m_elements = GetNumSDPBlocks();
      int nTmp = m_nie - m_mf;
	  int nMSize = m_blks[2];
	  if (nMSize == 2){
		 concavity_d_2x2_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
	  } else if (nMSize == 3){
		 concavity_d_3x3_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
	  } else if (nMSize == 6){
		 concavity_d_6x6_(m_xinit,&m_elements,m_h_org,&eps,&nTmp,&flag);
	  } 

    } else if (m_status > 1) {					 // error case!!!!     
	  cout<<"****   ERROR   **** IERR = "<< m_status;
	  cout<<"\n";
      return 1;  // TO BE REDIFINED!
    } else if (m_status == 1) {					 // error case!!!!     
	  cout<<"****   ERROR   **** IERR = "<< m_status<<endl;
	  cout<<"\n";
      m_status = 0;
    }


  } while (m_status < 0);


# endif 
  //bool optimal = true;
  //cout<<"return";
  return (m_status == 0);

}


double* 
ASCPSolver::GetSolution() {
  if(!m_nvars || !m_xinit)
    return 0;

  double *res = new double[m_nvars];

  for (int i=0; i<m_nvars; i++)
    res[i] = m_xinit[i]/100; //;/100
// sonja
  //
  ofstream myfile;
  myfile.open("example_scp.txt");
  int j = (m_nvars-1)/6;
  for (int i=0; i<j; i++){
	int k = (i-1)*6;
	double dum = m_xinit[k]+m_xinit[k+2]+m_xinit[k+5];
	myfile<<i;
	myfile<<"     ";
	myfile<<dum/100;
    myfile<<"\n";
  }
  myfile.close();
  return res;
}


