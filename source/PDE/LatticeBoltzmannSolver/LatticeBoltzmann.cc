#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
# ifdef _OPENMP
  #include <omp.h>
#endif

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ResultHandler.hh"
#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/StaticDriver.hh"
#include "Driver/Assemble.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "Driver/FormsContexts.hh"
#include "Optimization/Optimization.hh"
#include "PDE/LatticeBoltzmannSolver/LatticeBoltzmann.hh"
#include "Utils/Timer.hh"
#include "PDE/BasePDE.hh"


namespace CoupledField
{

using std::fstream;
using std::ios;

DECLARE_LOG(lbm)
DEFINE_LOG(lbm, "lbm")

// instantiation of the static elements
Enum<LatticeBoltzmann::Direction>         LatticeBoltzmann::directions;

LatticeBoltzmann::LatticeBoltzmann(int dim, int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, StdVector< StdVector<double> > uin, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency, bool srt, double omega_e, double omega_eps, double omega_q, double alpha_max)
{
  assert(dim == 2 || dim == 3);
  // n_q_: number of discrete directions in this model, e.g. 19 for D3Q19
  if (dim == 2) {
    assert(sizeZ == 1);
    n_q_ = 9;
  }
  else
    n_q_ = 19;

  dim_ = dim;
  sizeX_ = sizeX;
  sizeY_ = sizeY;
  sizeZ_ = sizeZ;
  ux_ = ux;
  uy_ = uy;
  uz_ = uz;
  omega_nu_ = omega; // this relaxation rate is directly related to the fluid's viscosity (for both SRT and MRT models)
  maxIter_ = maxIterations;
  maxTol_ = maxTolerance;
  writeFrequency_ = writeFrequency;
  numWriteResults_ = 0;
  srt_ = srt;
  omega_e_ = omega_e;
  omega_eps_ = omega_eps;
  omega_q_ = omega_q;
  alpha_max_ = alpha_max;

  nNodes_ = sizeX_ * sizeY_ * sizeZ_;

  u_in.Resize(nNodes_);
  u_in = uin;

  moments_.Resize(nNodes_);
  eqMoments_.Resize(nNodes_);

  plot_ = plot;

  //matrix of the probability distributions
  LOG_DBG(lbm) << "Allocating arrays for " << nNodes_ * n_q_ << " PDFs (" << (sizeof(double) * nNodes_ * n_q_ * 2.0 / 1024.0 / 1024.0) << " MiB)";

  pdfs_.Resize(2);
  pdfs_[0].Resize(nNodes_ * n_q_);
  pdfs_[1].Resize(nNodes_ * n_q_);

  // microVelDirections stores information about the 19 microscopic velocities/directions of D3Q19 model
  microVelDirections.Resize(n_q_);

  weights.Resize(n_q_);

  SetEnums();
  InitWeights();
  SetInvDirections();
  TestInvDirections();

  scales.Resize(nNodes_);
  rel.Resize(nNodes_);
  bb.Resize(2 * sizeX_ + 2 * sizeY_ + 2 * sizeZ_);

  cur_  = 0;
  next_ = 1;

  lbmCalls_ = 0;

  SetMicroVelocities();

  if (!srt_)
    InitTransformMatrix();

  //initlialize function pointers in dependence on problem's dimension
  if (dim == 2) {
    prop_coll_step = &LatticeBoltzmann::Prop_coll_step2D;
    prop_coll_velinlet = &LatticeBoltzmann::Prop_coll_velinlet2D;
    prop_coll_bounce_back = &LatticeBoltzmann::Prop_coll_bounce_back2D;
    prop_coll_densoutlet = &LatticeBoltzmann::Prop_coll_densoutlet2D;
  }
  else
  {
    prop_coll_step = &LatticeBoltzmann::Prop_coll_step3D;
    prop_coll_velinlet = &LatticeBoltzmann::Prop_coll_velinlet3D;
    prop_coll_bounce_back = &LatticeBoltzmann::Prop_coll_bounce_back3D;
    prop_coll_densoutlet = &LatticeBoltzmann::Prop_coll_densoutlet3D;
  }
}

LatticeBoltzmann::~LatticeBoltzmann()
{
}

void LatticeBoltzmann::CalcVelocities(const Vector<double>& pdfs, double& ux, double& uy, double& uz)
{

  double density = CalcDensity(pdfs);
  double tmp_ux = 0;
  double tmp_uy = 0;
  double tmp_uz = 0;

  for (int  dir = 0; dir < n_q_; dir++) {
    //store current pdf values in array for better accessing
    tmp_ux += microVelDirections[dir].off_x * pdfs[dir];
    tmp_uy += microVelDirections[dir].off_y * pdfs[dir];
    tmp_uz += microVelDirections[dir].off_z * pdfs[dir];
  }

  ux = tmp_ux / density;
  uy = tmp_uy / density;
  uz = tmp_uz / density;
}

StdVector<double>* LatticeBoltzmann::Iterate(const StdVector<double>& elements, PtrParamNode in)
{
  int it = 0;
  // counts number of written steps when not in optimization
  int count = 0;

  // this flag cannot be set in constructor, since information about optimization is not available when this constructor is called
  if (domain->GetOptimization() != NULL && domain->GetOptimization()->GetOptimizerType() != Optimization::EVALUATE_INITIAL_DESIGN) {
    writeIntermediateResults_ = false;
  }
  else
    writeIntermediateResults_ = true;

  double res = -1.;
  double R = 1.0;
  bool steady_state = false;

  std::ofstream plot;
  if(plot_)
    plot.open(std::string(progOpts->GetSimName() + ".lbm.dat").c_str());

  InitializePdfs();
  SetupDataStructures(elements);

  assert((int) elements.GetSize() == nNodes_);
  for (int i = 0; i < nNodes_; ++i) {
    scales[i] = 1.0 - elements[i];

    LOG_DBG3(lbm) << "Element " << i << " has density " << elements[i] << " and porosity " << scales[i];
  }



  Timer timer;
  timer.Start();

  LOG_DBG(lbm) << "bb = " << ToString(bb);
  LOG_DBG(lbm) << "inlet = " << ToString(inlet);
  LOG_DBG(lbm) << "outlet = " << ToString(outlet);
  LOG_DBG(lbm) << "rel = " << ToString(rel);

  in->Get("converged")->SetValue("running");

  while(it < maxIter_ && !steady_state && R <= 1000)
  {
    LOG_DBG3(lbm) << "---------------------------Iteration " << it << "---------------------------------------------------";
    // -- Combined propagation and collision step -------------------------
    (this->*prop_coll_step)(cur_, next_);
    // -- Bounce back step ------------------------------------------------
    (this->*prop_coll_bounce_back)(next_, bb);
    // -- Inlet condition -------------------------------------------------
    (this->*prop_coll_velinlet)(next_, inlet);
    // -- Outlet condition ------------------------------------------------
    (this->*prop_coll_densoutlet)(next_, outlet);
    if((it == 0 || it % 100 == 0))
    {
      //Calculation of the residual
      R = 0.;

      for (int elem = 0; elem < nNodes_; elem++) {
        //            index = k * m_sizeX * m_sizeY + j * m_sizeX + i;
        for(int  dir = 0; dir < n_q_; dir++)
        {
          res = PDF(next_, elem, dir) - PDF(cur_, elem, dir);
          R += res * res;
        }
      }

      R = sqrt(R);

      if(R <= maxTol_)
        steady_state = true;

      if(plot_)
        plot << it << "\t" << R << "\n";

      in->Get("iterations")->SetValue(it);
      in->Get("residuum")->SetValue(R);
      domain->GetInfoRoot()->ToFile(); // is not written when called too often
    }

    cur_  = (cur_  + 1) % 2;
    next_ = (next_ + 1) % 2;

    it++;

    if (writeIntermediateResults_) {
      if (it % writeFrequency_ == 0) {
        domain->GetDriver()->StoreResults(count,(double) it);
        count++;
      }
    }
//    LOG_DBG3(lbm) << "\n Iteration " << it;
//    for (int elem = 0; elem < nNodes_; elem++) {
//      LOG_DBG3(lbm) << "element " << elem;
//      for(int  dir = 0; dir < n_q_; dir++)
//        LOG_DBG3(lbm) << "dir " << dir << " pdf= " << PDF(next_,elem,dir) << " ";
//    }
  }

  timer.Stop();

  PtrParamNode node = in->Get(ParamNode::PROCESS)->Get("call", ParamNode::APPEND); // write out how many lbm iterations until convergence
  node->Get("number")->SetValue(lbmCalls_);
  node->Get("iterations")->SetValue(it);
  node->Get("residuum")->SetValue(R);
  node->Get("converged")->SetValue(steady_state);

  if(plot_) {
    plot << it << "\t" << R << "\n";
    plot.flush();
  }

  double wt = timer.GetWallTime();
  double performance;
  if (dim_ == 3)
    performance = (sizeX_ - 1) * (sizeY_ - 1) * (sizeZ_ - 1) * it / wt / 1e6;
  else
    performance = (sizeX_ - 1) * (sizeY_ - 1) * it / wt / 1e6;

  node->Get("wall")->SetValue(wt);
  node->Get("cpu")->SetValue(timer.GetCPUTime());
  node->Get("MFLUP_s")->SetValue(performance);
  in->Get("memory")->SetValue(sizeof(double) * nNodes_ * n_q_ * 2.0 / 1024.0 / 1024.0);

  if(R >= 1000)
    EXCEPTION("In LBM iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("internal LBM simulation could not converge: iterations: " << it << " residuum: " << R);

  numIterations_ = it;
  numWriteResults_ = count;

  lbmCalls_++; // first solver call is call number 0 (to match iteration numbering of optimizer)

  if (!srt_)
  {
    Vector<double> pdfs;
    pdfs.Resize(9);
    double dissipation = 0.0;
    for (int elem = 0; elem < nNodes_; elem++) {
      for (int dir = 0; dir < n_q_; dir++)
        pdfs[dir] = PDF(cur_,elem,dir);

      Vector<double> moms;
      moms.Resize(n_q_);
      transformation.Mult(pdfs,moms);
      Vector<double> eqMoms;
      CalcEquilMoments(moms,eqMoms);

      // store moments and equilibrium moments of steady-state solution
      for (int dir = 0; dir < n_q_; dir++) {
        moments_[dir] = moms[dir];
        eqMoments_[dir] = eqMoms[dir];
      }
      Vector<double> f1,f2;
      CalcDarcyForce(moms,elem,f1,f2);
      dissipation += CalcDissipation(moms,eqMoms,f1[3],f1[5]); // FIXME Add body forces
    }
    node->Get("dissipation")->SetValue(dissipation);
  }


  return &(pdfs_[cur_]);
}

void LatticeBoltzmann::InitializePdfs()
{
#pragma omp parallel for default(none)
  for (int elem = 0; elem < nNodes_; elem++) {
    for (int  dir = 0; dir < n_q_; dir++) {
      PDF(0, elem, dir) = weights[dir];
      PDF(1, elem, dir) = weights[dir];
    }
  }
}

void LatticeBoltzmann::SetMicroVelocities()
{
  microVelDirections[Q_0] = PDFDirectionVector(0,0,0);
  microVelDirections[Q_E] = PDFDirectionVector(1,0,0);
  microVelDirections[Q_W] = PDFDirectionVector(-1,0,0);
  microVelDirections[Q_N] = PDFDirectionVector(0,1,0);
  microVelDirections[Q_S] = PDFDirectionVector(0,-1,0);
  microVelDirections[Q_NE] = PDFDirectionVector(1,1,0);
  microVelDirections[Q_NW] = PDFDirectionVector(-1,1,0);
  microVelDirections[Q_SE] = PDFDirectionVector(1,-1,0);
  microVelDirections[Q_SW] = PDFDirectionVector(-1,-1,0);
  if (dim_ == 3) {
    microVelDirections[Q_T] = PDFDirectionVector(0,0,1);
    microVelDirections[Q_B] = PDFDirectionVector(0,0,-1);
    microVelDirections[Q_TN] = PDFDirectionVector(0,1,1);
    microVelDirections[Q_TS] = PDFDirectionVector(0,-1,1);
    microVelDirections[Q_TE] = PDFDirectionVector(1,0,1);
    microVelDirections[Q_TW] = PDFDirectionVector(-1,0,1);
    microVelDirections[Q_BN] = PDFDirectionVector(0,1,-1);
    microVelDirections[Q_BS] = PDFDirectionVector(0,-1,-1);
    microVelDirections[Q_BE] = PDFDirectionVector(1,0,-1);
    microVelDirections[Q_BW] = PDFDirectionVector(-1,0,-1);
  }
}

void LatticeBoltzmann::InitTransformMatrix()
{
  LOG_TRACE(lbm) << "Initializing transformation matrix M and relaxation rates matrix S.";

  // transformation matrix M
  transformation.Resize(n_q_);
  transformation.Init(); // initialize with zeros
  // density mode
  for (int  i = 0; i < n_q_; i++)
  {
    transformation[0][i] = 1.0;
  }
  // energy mode
  transformation[1][0] = -4.0;
  transformation[1][1] = -1.0;
  transformation[1][2] = -1.0;
  transformation[1][3] = -1.0;
  transformation[1][4] = -1.0;
  transformation[1][5] = 2;
  transformation[1][6] = 2;
  transformation[1][7] = 2;
  transformation[1][8] = 2;
  // energy square
  transformation[2][0] = 4.0;
  transformation[2][1] = -2.0;
  transformation[2][2] = -2.0;
  transformation[2][3] = -2.0;
  transformation[2][4] = -2.0;
  transformation[2][5] = 1.0;
  transformation[2][6] = 1.0;
  transformation[2][7] = 1.0;
  transformation[2][8] = 1.0;
  // x momentum
  transformation[3][1] = 1.0;
  transformation[3][3] = -1.0;
  transformation[3][5] = 1.0;
  transformation[3][6] = -1.0;
  transformation[3][7] = -1.0;
  transformation[3][8] = 1.0;
  // x energy flux
  transformation[4][1] = -2.0;
  transformation[4][3] = 2.0;
  transformation[4][5] = 1.0;
  transformation[4][6] = -1.0;
  transformation[4][7] = -1.0;
  transformation[4][8] = 1.0;
  // y momentum
  transformation[5][2] = 1.0;
  transformation[5][4] = -1.0;
  transformation[5][5] = 1.0;
  transformation[5][6] = 1.0;
  transformation[5][7] = -1.0;
  transformation[5][8] = -1.0;
  // y energy flux
  transformation[6][2] = -2.0;
  transformation[6][4] = 2.0;
  transformation[6][5] = 1.0;
  transformation[6][6] = 1.0;
  transformation[6][7] = -1.0;
  transformation[6][8] = -1.0;
  // stress tensor diagonal
  transformation[7][1] = 1.0;
  transformation[7][2] = -1.0;
  transformation[7][3] = 1.0;
  transformation[7][4] = -1.0;
  //stress tensor off-diagonal
  transformation[8][5] = 1.0;
  transformation[8][6] = -1.0;
  transformation[8][7] = 1.0;
  transformation[8][8] = -1.0;

  LOG_DBG3(lbm) << "IT: M=" << transformation.ToString(3,true);

//  std::cout << "---------------Transformation matrix--------------------" << std::endl;
//  std::cout << transformation.ToString(0,true) << std::endl;

  // sanity check
  for (int i = 1; i < n_q_; i++)
  {
    double sum = 0;
    for (int j = 0; j < n_q_; j++)
    {
      sum += transformation[i][j];
    }
    assert(sum == 0);
  }
  invTransformation.Resize(n_q_);
  invTransformation.Init();

  for (int i = 0; i < n_q_; i++)
    invTransformation[i][0] = 2.0;

  invTransformation[0][1] = -2.0;
  invTransformation[0][2] = 2.0;
  invTransformation[1][1] = -0.5;
  invTransformation[1][2] = -1.0;
  invTransformation[1][3] = 3.0;
  invTransformation[1][4] = -3.0;
  invTransformation[1][7] = 4.5;
  invTransformation[2][1] = -0.5;
  invTransformation[2][2] = -1.0;
  invTransformation[2][5] = 3.0;
  invTransformation[2][6] = -3.0;
  invTransformation[2][7] = -4.5;
  invTransformation[3][1] = -0.5;
  invTransformation[3][2] = -1.0;
  invTransformation[3][3] = -3.0;
  invTransformation[3][4] = 3.0;
  invTransformation[3][7] = 4.5;
  invTransformation[4][1] = -0.5;
  invTransformation[4][2] = -1.0;
  invTransformation[4][5] = -3.0;
  invTransformation[4][6] = 3;
  invTransformation[4][7] = -4.5;
  invTransformation[5][1] = 1.0;
  invTransformation[5][2] = 0.5;
  invTransformation[5][3] = 3.0;
  invTransformation[5][4] = 1.5;
  invTransformation[5][5] = 3.0;
  invTransformation[5][6] = 1.5;
  invTransformation[5][8] = 4.5;
  invTransformation[6][1] = 1.0;
  invTransformation[6][2] = 0.5;
  invTransformation[6][3] = -3.0;
  invTransformation[6][4] = -1.5;
  invTransformation[6][5] = 3.0;
  invTransformation[6][6] = 1.5;
  invTransformation[6][8] = -4.5;
  invTransformation[7][1] = 1.0;
  invTransformation[7][2] = 0.5;
  invTransformation[7][3] = -3.0;
  invTransformation[7][4] = -1.5;
  invTransformation[7][5] = -3.0;
  invTransformation[7][6] = -1.5;
  invTransformation[7][8] = 4.5;
  invTransformation[8][1] = 1.0;
  invTransformation[8][2] = 0.5;
  invTransformation[8][3] = 3.0;
  invTransformation[8][4] = 1.5;
  invTransformation[8][5] = -3.0;
  invTransformation[8][6] = -1.5;
  invTransformation[8][8] = -4.5;

  invTransformation /= 18.0;

//  Matrix<Double> identity;
//  identity.Resize(n_q_);
//  invTransformation.Mult(transformation,identity);
//  std::cout << identity.ToString(0,true)<< std::endl;

//  std::cout << "---------------Inverse transformation matrix--------------------" << std::endl;
//  std::cout << invTransformation.ToString(0,true) << std::endl;

  relax_rates.Resize(n_q_);
  relax_rates.Init();
  // relaxation rates matrix S is diagonal
  // entries for density and momentum components are 0 due to collision invariance
  // the first 4 values here, can be chosen arbritrarily between 0 and 2
  relax_rates[1] = omega_e_;
  relax_rates[2] = omega_eps_;
  relax_rates[4] = omega_q_;
  relax_rates[6] = omega_q_;
  // these two values are related to fluid's kinematic viscosity
  relax_rates[7] = omega_nu_;
  relax_rates[8] = omega_nu_;

  LOG_TRACE(lbm) << "MRT relaxation rates are [0," << omega_e_ << "," << omega_eps_ << ",0," << omega_q_ << "," << omega_nu_ << "," << omega_nu_ << "]";
  invM_S.Resize(n_q_);
  invM_S.Init();
  // multiply inverse of M with relaxation rates: M^{-1} * S
  for (int i = 0; i < n_q_; i++)
  {
    for (int j = 0; j < n_q_; j++)
    {
      invM_S[i][j] = relax_rates[j] * invTransformation[i][j];
    }
  }

//    std::cout << "---------------M^-1 S--------------------" << std::endl;
//    std::cout << invM_S.ToString(0,true) << std::endl;

}

void LatticeBoltzmann::SetEnums()
{
  directions.SetName("Q19 directions");
  directions.Add(Q_0,"C");
  directions.Add(Q_E,"E");
  directions.Add(Q_N,"N");
  directions.Add(Q_W,"W");
  directions.Add(Q_S,"S");
  directions.Add(Q_T,"T");
  directions.Add(Q_B,"B");
  directions.Add(Q_NE,"NE");
  directions.Add(Q_NW,"NW");
  directions.Add(Q_SW,"SW");
  directions.Add(Q_SE,"SE");
  directions.Add(Q_TN,"TN");
  directions.Add(Q_BS,"BS");
  directions.Add(Q_TS,"TS");
  directions.Add(Q_BN,"BN");
  directions.Add(Q_TE,"TE");
  directions.Add(Q_BW,"BW");
  directions.Add(Q_TW,"TW");
  directions.Add(Q_BE,"BE");
}

void LatticeBoltzmann::InitWeights()
{
  weights.Resize(n_q_);

  if (dim_ == 3) {
    weights[Q_0] = 1. / 3.;
    weights[Q_E] = 1. / 18.;
    weights[Q_N] = 1. / 18.;
    weights[Q_W] = 1. / 18.;
    weights[Q_S] = 1. / 18.;
    weights[Q_T] = 1. / 18.;
    weights[Q_B] = 1. / 18.;
    weights[Q_NE] = 1. / 36.;
    weights[Q_NW] = 1. / 36.;
    weights[Q_SW] = 1. / 36.;
    weights[Q_SE] = 1. / 36.;
    weights[Q_TN] = 1. / 36.;
    weights[Q_BS] = 1. / 36.;
    weights[Q_TS] = 1. / 36.;
    weights[Q_BN] = 1. / 36.;
    weights[Q_TE] = 1. / 36.;
    weights[Q_BW] = 1. / 36.;
    weights[Q_TW] = 1. / 36.;
    weights[Q_BE] = 1. / 36.;
  }
  else {
    weights[Q_0] = 4.0 / 9.0;
    weights[Q_E] = 1.0 /  9.0;
    weights[Q_N] = 1.0 /  9.0;
    weights[Q_W] = 1.0 /  9.0;
    weights[Q_S] = 1.0 /  9.0;
    weights[Q_NE] = 1.0 / 36.0;
    weights[Q_NW] = 1.0 / 36.0;
    weights[Q_SW] = 1.0 / 36.0;
    weights[Q_SE] = 1.0 / 36.0;
  }
}

void LatticeBoltzmann::SetInvDirections()
{
  invPDFDirections.Resize(n_q_);

  invPDFDirections[Q_0] = Q_0;
  invPDFDirections[Q_E] = Q_W;
  invPDFDirections[Q_N] = Q_S;
  invPDFDirections[Q_W] = Q_E;
  invPDFDirections[Q_S] = Q_N;
  invPDFDirections[Q_NE] = Q_SW;
  invPDFDirections[Q_SW] = Q_NE;
  invPDFDirections[Q_NW] = Q_SE;
  invPDFDirections[Q_SE] = Q_NW;

  if (dim_ == 3) {
    invPDFDirections[Q_T] = Q_B;
    invPDFDirections[Q_B] = Q_T;
    invPDFDirections[Q_BS] = Q_TN;
    invPDFDirections[Q_TS] = Q_BN;
    invPDFDirections[Q_TN] = Q_BS;
    invPDFDirections[Q_BN] = Q_TS;
    invPDFDirections[Q_TE] = Q_BW;
    invPDFDirections[Q_BW] = Q_TE;
    invPDFDirections[Q_TW] = Q_BE;
    invPDFDirections[Q_BE] = Q_TW;
  }
}

void LatticeBoltzmann::TestInvDirections()
{
  assert(GetInvDirection(Q_0) == Q_0);
  assert(GetInvDirection(Q_E) == Q_W);
  assert(GetInvDirection(Q_N) == Q_S);
  assert(GetInvDirection(Q_W) == Q_E);
  assert(GetInvDirection(Q_S) == Q_N);
  assert(GetInvDirection(Q_NE) == Q_SW);
  assert(GetInvDirection(Q_NW) == Q_SE);
  assert(GetInvDirection(Q_SW) == Q_NE);
  assert(GetInvDirection(Q_SE) == Q_NW);

  if (dim_ == 3) {
    assert(GetInvDirection(Q_T) == Q_B);
    assert(GetInvDirection(Q_B) == Q_T);
    assert(GetInvDirection(Q_TN) == Q_BS);
    assert(GetInvDirection(Q_TS) == Q_BN);
    assert(GetInvDirection(Q_BS) == Q_TN);
    assert(GetInvDirection(Q_BN) == Q_TS);
    assert(GetInvDirection(Q_TE) == Q_BW);
    assert(GetInvDirection(Q_BW) == Q_TE);
    assert(GetInvDirection(Q_TW) == Q_BE);
    assert(GetInvDirection(Q_BE) == Q_TW);
  }
}

void LatticeBoltzmann::SetupDataStructures(const StdVector<double>& elements)
{
  rel.Clear();
  bb.Clear();
  inlet.Clear();
  outlet.Clear();

  StdVector<int> tmp(3);

  int n = 0;
  double porosity;

  for(int k = 0; k < sizeZ_; k++)
  {
    for(int j = 0; j < sizeY_; j++)
    {
      for(int i = 0; i < sizeX_; i++)
      {
        porosity = elements[n];

        if (LbmNodeTypeIsFluid(porosity)) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          rel.Push_back(tmp);
//          std::cout << "(" << tmp[0] << "," <<  tmp[1] << "," << tmp[2] << ") is fluid node" << std::endl;
        } else if (LbmNodeTypeIsBB(porosity)) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          bb.Push_back(tmp);
//          std::cout << "(" << tmp[0] << "," <<  tmp[1] << "," << tmp[2] << ") is bb node" << std::endl;
        } else if (LbmNodeTypeIsInlet(porosity)) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          inlet.Push_back(tmp);
//          std::cout << "(" << tmp[0] << "," <<  tmp[1] << "," << tmp[2] << ") is inlet node" << std::endl;
        } else if (LbmNodeTypeIsOutlet(porosity)) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          outlet.Push_back(tmp);
//          std::cout << "(" << tmp[0] << "," <<  tmp[1] << "," << tmp[2] << ") is outlet node" << std::endl;
        } else if (LbmNodeIsObstacle(porosity)) {
          obst.Push_back(GetIndex(i,j,k));
//          std::cout << "(" << tmp[0] << "," <<  tmp[1] << "," << tmp[2] << ") is obstacle node" << std::endl;
        }
        ++n;
      }
    }
  }
  assert(nNodes_ == n);
}


std::string LatticeBoltzmann::ToString(const StdVector<StdVector<int> >& data)
{
  std::stringstream ss;

  for(unsigned int  e = 0; e < data.GetSize(); e++)
  {
    ss << e << ": (";
    for(unsigned int  d = 0; d < data[e].GetSize(); d++)
    {
      ss << data[e][d];
      if(d < data[e].GetSize() - 1)
        ss << ", ";
    }
    ss << ") ";
  }
  return ss.str();
}

void LatticeBoltzmann::CreateOutput(const char * file, int cur)
{
  // for debug purposes
  fstream f;
  f.precision(16);
  f.open(file, ios::out);

  for(int i = 0; i < nNodes_; i++) {
    for(int  j = 0; j < n_q_; j++) {
      f << PDF(cur, i, j) << " ";
    }

    f << std::endl;
  }
  f.close();
}


void LatticeBoltzmann::Prop_step()
{
  int tmp_x, tmp_y, tmp_z;
  // perform a propagation step
  for (int z = 0; z < sizeZ_; ++z) {
    for (int y = 0; y < sizeY_ ; ++y) {
      for (int x = 0; x < sizeX_ ; ++x) {

        for (int  dir = 0; dir < n_q_; dir++) {
          int invDir = GetInvDirection((Direction)dir);
          if (PointsToBoundary(x,y,z,invDir)) { // if the neighbor element that I want to access is outside the domain, keep current value
            tmp_x = x; // here we only set the coordinates
            tmp_y = y;
            tmp_z = z;
          }
          // else: standard propagation (get value from neighbor pdf)
          else {
            tmp_x = microVelDirections[invDir].off_x + x;
            tmp_y = microVelDirections[invDir].off_y + y;
            tmp_z = microVelDirections[invDir].off_z + z;
          }

          PDF(next_,x,y,z,dir) = PDF(cur_, tmp_x, tmp_y,  tmp_z, dir);
        }
      }
    }
  }

  pdfs_[cur_] = pdfs_[next_];
}

// Calculates macroscopic density for given element
double LatticeBoltzmann::CalcDensity(const Vector<double>& pdfs)
{
  double sum = 0;
  for (int  dir = 0; dir < n_q_; dir++) {
    sum += pdfs[dir];
  }

  // debugging
//  if (sum > 1.5 || sum < 1e-8) {
//    std::cout << "i: " << i << " j: " << j << " k: " << j << " sum: " << sum << std::endl;
//    for (int  dir = 0; dir < n_q_; dir++) {
//      std::cout << "dir " << dir << " pdf " <<  pdf(cur, i, j, k, dir) << std::endl;
//    }
//    EXCEPTION("LBM simulation has problems, macroscopic densities are too big!");
//  }
//  assert(sum < 1.5);
//  assert(sum > 1e-8);
  return sum;
}

void LatticeBoltzmann::CalcDarcyForce(const Vector<double>& moments, int elemId, Vector<double>& f1, Vector<double>& f2)
{
  f1.Resize(n_q_);
  f1.Init(0.0);
  f2.Resize(n_q_);
  f2.Init(0.0);

//  double q = 1; // this value was is always set to 1 in accordance to the paper of Liu et al. (2014); but can also be chosen differently
//  double nu_water =1.004e-6; // unit m^2/s
//  double characLength = 1e-3; // 1 mm
//  double nu_lb = 0.1;  // nu = 1/3*(1/omega_nu - 0.5)
//  double time_ref = 1/3.0 * (1/1.9-0.5) * 1e-2/(sizeX_*sizeX_) / nu_water; // reference time unit for conversion from physical units to LBM ones
//  double time_ref = nu_lb * characLength * characLength / nu_water;
//  std::cout << "time_ref =" << time_ref << " = " <<  nu_lb << "*" << characLength << "*" << characLength << "/" << nu_water << std::endl;
  //double alpha_lb = alpha_max_ * time_ref;
//  double alpha_lb = 1.7;

//  std::cout << "alpha_max in LB units: " << alpha_lb << " time scale: " << time_ref << std::endl;
  // alpha = alphaMax * q(1-gamma) / (q + gamma), where gamma is the porosity (0 is solid and 1 is fluid)
//  double alpha = alpha_lb * q * (1 - scales[elemId]) / (q +  scales[elemId]); // TODO Check if mapping between design variable and porosity is correct
  double alpha = CalcResistanceCoeff(elemId);

  double rho = moments[0];
  double ux = moments[3];
  double uy = moments[5];

  double fx = -alpha * rho * ux;
  double fy = -alpha * rho * uy;

  LOG_DBG3(lbm) << " CDF: Element " << elemId << " has LBM density " << rho << " and porosity " << scales[elemId];
  LOG_DBG3(lbm) << "CDF: fx=" << fx << "= -" << alpha << "*" << rho << "*" << ux;
  LOG_DBG3(lbm) << "CDF: alpha= " << alpha << " fx=" << fx << " fy=" << fy;
//  if (scales[elemId] == 0.0) {
//  for(int dir = 0; dir < n_q_; dir++)
//    LOG_DBG3(lbm) << pdfs[dir];
//  }

  f1[3] = fx; // the assembly of f1 and f2 are taken from MRT paper
  f1[4] = -f1[3];
  f1[5] = fy;
  f1[6] = -f1[5];

  f2[1] = 6 * (fx * ux + fy * uy);
  f2[2] = -f2[1];
  f2[7] = 2 * (fx * ux - fy * uy);
  f2[8] = fy * ux + fx * uy;

  //debugging
//  if (scales[elemId] == 0.0) {
//    for(int dir = 0; dir < n_q_; dir++)
//      std::cout << pdfs[dir] << " ";
//  }

  LOG_DBG3(lbm) << "CDF: F1=" << f1.ToString(0,',');
  LOG_DBG3(lbm) << "CDF: F2=" << f2.ToString(0,',') << "\n";
}

double LatticeBoltzmann::CalcResistanceCoeff(int elemId)
{
  return alpha_max_ * (1 - scales[elemId]) / (1 +  scales[elemId]);
}

// m_eq = M * f_eq
// Here, m_eq is calculated via given formulas avoiding matrix multiplication
// density and momentum (velocity) can be extracted directly from moments vector
void LatticeBoltzmann::CalcEquilMoments(const Vector<double>&  moments, Vector<double>& m_eq) const
{
  m_eq.Resize(n_q_);
  double density = moments[0];
  double ux = moments[3] / density;
  double uy = moments[5] / density;
  double uz = 0; // TODO: Change this for 3D simulation

  if (dim_ == 2)
    assert(uz == 0);

  double us = ux * ux + uy * uy + uz * uz;

  m_eq[0] = 1;
  m_eq[1] = -2 + 3 * us;
  m_eq[2] = 1 - 3 * us;
  m_eq[3] = ux;
  m_eq[4] = -ux;
  m_eq[5] = uy;
  m_eq[6] = -uy;
  m_eq[7] = ux * ux - uy * uy;
  m_eq[8] = ux * uy;

  m_eq.ScalarMult(density);

  //std::cout << "equil moments: " << m_eq.ToString(0,',') << std::endl;
}

double LatticeBoltzmann::CalcDissipation(const Vector<double>& moments, const Vector<double>& eqMoments, double fx, double fy)
{
  double nu = (1/omega_nu_ - 0.5) / 3.0; // fluid's kinematic viscosity in LBM terms
  double term1 = (eqMoments[1] - moments[1]) * (eqMoments[1] - moments[1]);
  double term2 = (eqMoments[7] - moments[7]) * (eqMoments[7] - moments[7]);
  double term3 = (eqMoments[8] - moments[8]) * (eqMoments[8] - moments[8]);
  double ux =  moments[5] / moments[0];
  double uy =  moments[7] / moments[0];
  return nu / moments[0] *  (0.25 * (omega_e_ * omega_e_ * term1 + 9 * omega_nu_ * omega_nu_ * term2) + 9 * omega_nu_ * omega_nu_ * term3 ) - fx * ux - fy * uy;
}

void LatticeBoltzmann::CalcAdjointCollMatrix(int elemId, const Vector<double>& moments, Matrix<double>& out)
{
  // S_A = I - S(I - d_mEq/d_m) + d_F1/d_m + (I - S/2) d_F2/d_m

  double rho = moments[0];
  double rho2Inv = 1 / (rho * rho);
  double jx = moments[3];
  double jy = moments[5];
  double ux = jx / rho;
  double uy = jy / rho;
  double u2 = ux * ux + uy * uy;

  // non-zero entries of d_mEq/d_m
  Vector<double> d_mEq_d_rho(n_q_);
  d_mEq_d_rho.Init(0.0);
  d_mEq_d_rho[0] = 1.0;
  d_mEq_d_rho[1] = - 2.0 - 3 * u2;
  d_mEq_d_rho[2] = 1 + 3 * u2;
  d_mEq_d_rho[7] = - ux * ux + uy * uy;
  d_mEq_d_rho[8] = - ux * uy;
  Vector<double> d_mEq_djx(n_q_);
  d_mEq_djx.Init(0.0);
  d_mEq_djx[1] = 6 * ux;
  d_mEq_djx[2] = -6 * ux;
  d_mEq_djx[3] = 1.0;
  d_mEq_djx[4] = -1.0;
  d_mEq_djx[7] = 2 * ux;
  d_mEq_djx[8] = uy;
  Vector<double> d_mEq_djy(n_q_);
  d_mEq_djy.Init(0.0);
  d_mEq_djy[1] = 6 * uy;
  d_mEq_djy[2] = - 6 * uy;
  d_mEq_djy[5] = 1.0;
  d_mEq_djy[6] = -1.0;
  d_mEq_djy[7] = -2*uy;
  d_mEq_djy[8] = ux;

  Vector<double> d_F1_d_jx(n_q_),d_F1_d_jy(n_q_);
  d_F1_d_jx.Init(0.0);
  d_F1_d_jy.Init(0.0);
  double alpha = CalcResistanceCoeff(elemId);
  d_F1_d_jx[3] = -alpha;
  d_F1_d_jx[4] = alpha;
  d_F1_d_jy[5] = -alpha;
  d_F1_d_jy[6] = alpha;

  Vector<double> d_F2_d_jx(n_q_),d_F2_d_jy(n_q_);
  d_F2_d_jx.Init(0.0);
  d_F2_d_jy.Init(0.0);
  d_F2_d_jx[1] = -12 * jx;
  d_F2_d_jx[2] = 12 * jx;
  d_F2_d_jx[7] = -4 * jx;
  d_F2_d_jx[8] = -2 * jy;
  d_F2_d_jx.ScalarMult(alpha/rho);
  d_F2_d_jy[1] = -12 * jy;
  d_F2_d_jy[2] = 12 * jy;
  d_F2_d_jy[7] = 4 * jy;
  d_F2_d_jy[8] = -2 * jx;
  d_F2_d_jy.ScalarMult(alpha/rho);

}
/************************************************** 2D operators *****************************************************/

void LatticeBoltzmann::Prop_coll_step2D(int cur, int next)
{
  int x, y, z = 0;

  double tmp, tmp_ux, tmp_uy, tmp_us, scale, sum;

  int index;

  int tmp_x,tmp_y;

#pragma omp parallel default(none)\
    private(index), \
    private(tmp_ux, tmp_uy, tmp_us, scale, sum, tmp, x, y, tmp_x, tmp_y), shared(next, cur, z)
  {
    Vector<double> pdfs;
    pdfs.Resize(n_q_);
#pragma omp for collapse(2)
    Vector<double> m_eq(n_q_);
    Vector<double> moments(n_q_);
    for (y = 0; y < sizeY_ ; ++y) {
      for (x = 0; x < sizeX_ ; ++x) {

        index= GetIndex(x,y,z);

        // sum: macroscopic density is sum over all discrete distributions of an element
        sum = 0;
        tmp_ux = 0;
        tmp_uy = 0;
        tmp_us = 0;

        // propagation
        for (int  dir = 0; dir < n_q_; dir++) {
          int invDir = GetInvDirection((Direction)dir);
          if (PointsToBoundary(x,y,z,invDir)) { // if the neighbor element that I want to access is outside the domain, keep current value
            tmp_x = x; // here we only set the coordinates
            tmp_y = y;
          }
          // else: standard propagation (get value from neighbor pdf)
          else {
            tmp_x = microVelDirections[invDir].off_x + x;
            tmp_y = microVelDirections[invDir].off_y + y;
          }

          //store current pdf values in array for better accessing
          pdfs[dir] = PDF(cur, tmp_x, tmp_y,  z, dir); // accessed pdf value can be the old one or the neighbor's value
          sum += pdfs[dir];
          tmp_ux += microVelDirections[dir].off_x*pdfs[dir];
          tmp_uy += microVelDirections[dir].off_y*pdfs[dir];
        }

        // macroscopic scaling by design variable
        scale = scales[index];

        Vector<double> collResult;
        if (srt_)
        {
          tmp_ux = scale * tmp_ux / sum;
          tmp_uy = scale * tmp_uy / sum;
          tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
          tmp_ux = 3.0 * tmp_ux;
          tmp_uy = 3.0 * tmp_uy;

          // propagation and collision in one step
          for (int  dir = 0; dir < n_q_; dir++)
          {
            tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy ;
            // no collision on the boundaries
            if (x == 0 || y == 0 || x == sizeX_ - 1 || y == sizeY_ - 1)
              PDF(next, x, y, z, dir) = pdfs[dir];
            else
              PDF(next, x, y, z, dir) = pdfs[dir] + omega_nu_ * (sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us) - pdfs[dir]);
          }
        }
        else // MRT case
        {
          // no collision on the boundaries
          if (x == 0 || y == 0 || x == sizeX_ - 1 || y == sizeY_ - 1)
          {
            for (int  dir = 0; dir < n_q_; dir++)
              PDF(next, x, y, z, dir) = pdfs[dir];
            continue;
          }

          transformation.Mult(pdfs,moments); // transform moments  m = M*f
          CalcEquilMoments(moments, m_eq); // compute equilibirum moments from moments
//          subtrahend.Resize(n_q_);
//          invM_S.Mult(moments - m_eq,subtrahend); // back transformation with relaxation

          Vector<double> momentsAfterCollision; // result of collision step in moment space including porosity model
          Vector<double> term1(n_q_);
          Vector<double> f1, f2;
          CalcDarcyForce(moments,index,f1,f2);
          Vector<double> noneq_moments(n_q_);
          for (int dir = 0; dir < n_q_; dir++)
          {
            noneq_moments[dir] = moments[dir] - m_eq[dir];
          }

          // S * (m - m_eq)
          for (int dir = 0; dir < n_q_; dir++) {
            term1[dir] = relax_rates[dir] * noneq_moments[dir];
          }

          Vector<double> term2(n_q_);
          for (int dir = 0; dir < n_q_; dir++) {
            term2[dir] = (1.0 - 0.5 * relax_rates[dir]) * f2[dir];
          }

          // m* = m - S * (m - m_eq) + F1 + (I - S/2) * F2
          momentsAfterCollision = moments - term1 + f1 + term2;

          collResult.Resize(n_q_);

          invTransformation.Mult(momentsAfterCollision,collResult);

//          for (int dir = 0; dir < n_q_; dir++) {
//            LOG_DBG3(lbm) << "f=" << collResult[dir] << " m* = " << momentsAfterCollision[dir] << "=" << moments[dir] << "-" << relax_rates[dir] << "*(" << moments[dir] << "-" << m_eq[dir] <<")"  << "+" << f1[dir] << "+" << term2[dir];
//          }

          // propagation and collision in one step
          for (int  dir = 0; dir < n_q_; dir++)
            PDF(next, x, y, z, dir) = collResult[dir];
        }

      }
    }
  }
  return;
}


void LatticeBoltzmann::Prop_coll_velinlet2D(int cur, StdVector<StdVector<int> >& inlet)
{
  assert(uz_ == 0);

  int x, y, z = 0;
  double tmp_ux, tmp_uy, tmp_us, sum;
  double tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];

    assert(z == 0);

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur,x,y,z,dir);
    }

    sum = CalcDensity(pdfs);

    tmp_ux = ux_; // velocity at inlet is prescribed
    tmp_uy = uy_;

    if (u_in.GetSize() != 0) {
      tmp_ux = u_in[i][0];
      tmp_uy = u_in[i][1];
    }

    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy;
      PDF(cur, x, y, z, dir)  = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }
  return;
}
//
// Performs a bounce back step.
//
void LatticeBoltzmann::Prop_coll_bounce_back2D(int cur, StdVector<StdVector<int> >& bb)
{
  int x, y, z = 0;
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur, x, y, z, dir);
    }
    for (int  dir = 0; dir < n_q_; dir++) {
      PDF(cur, x, y, z, GetInvDirection((Direction)dir)) = pdfs[dir];
    }
  }

  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::Prop_coll_densoutlet2D(int cur, StdVector<StdVector<int> >& outlet)
{
  double tmp_ux, tmp_uy, tmp_us, sum, tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  int x, y, z = 0;

  for(unsigned int  i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];

    assert (z == 0);

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur,x,y,z,dir);
    }

    CalcVelocities(pdfs, tmp_ux, tmp_uy, tmp);

    sum = 1.0; // the enforced density

    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy;
      PDF(cur, x, y, z, dir) =  sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }
  return;
}


/************************************************** 3D operators *****************************************************/
void LatticeBoltzmann::Prop_coll_step3D(int cur, int next)
{

  int x, y, z;
  double tmp, tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum;

  int index;

  int tmp_x, tmp_y, tmp_z;

#pragma omp parallel default(none)\
    private(index), \
    private(tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum, tmp, x, y, z, tmp_x, tmp_y, tmp_z), \
    shared(next, cur)
{
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  #pragma omp for collapse(3)
  for (z = 0; z < sizeZ_; ++z) {
    for (y = 0; y < sizeY_; ++y) {
      for (x = 0; x < sizeX_; ++x) {
        index= GetIndex(x,y,z);

        // sum: macroscopic density is sum over all discrete distributions of an element
        sum = 0;
        tmp_ux = 0;
        tmp_uy = 0;
        tmp_uz = 0;
        tmp_us = 0;

        for (int  dir = 0; dir < n_q_; dir++) {
          int invDir = GetInvDirection((Direction)dir);
          if (PointsToBoundary(x,y,z,invDir)) { // boundary case
            tmp_x = x;
            tmp_y = y;
            tmp_z = z;
          }
          else { // standard propagation rule
            tmp_x = microVelDirections[invDir].off_x + x;
            tmp_y = microVelDirections[invDir].off_y + y;
            tmp_z = microVelDirections[invDir].off_z + z;
          }

          //store current pdf values in array for better accessing
          pdfs[dir] = PDF(cur, tmp_x, tmp_y, tmp_z, dir);
          sum += pdfs[dir];
          tmp_ux += microVelDirections[dir].off_x*pdfs[dir];
          tmp_uy += microVelDirections[dir].off_y*pdfs[dir];
          tmp_uz += microVelDirections[dir].off_z*pdfs[dir];
        }

        // macroscopic scaling by design variable
        scale = scales[index];

        tmp_ux = scale * tmp_ux / sum;
        tmp_uy = scale * tmp_uy / sum;
        tmp_uz = scale * tmp_uz / sum;
        tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);

        tmp_ux = 3.0 * tmp_ux;
        tmp_uy = 3.0 * tmp_uy;
        tmp_uz = 3.0 * tmp_uz;

        for (int  dir = 0; dir < n_q_; dir++) {
          tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
          if (x == 0 || y == 0 || x == sizeX_ - 1 || y == sizeY_ - 1 || z == 0 || z == sizeZ_ - 1)
            PDF(next, x, y, z, dir) = pdfs[dir];
          else
            PDF(next, x, y, z, dir) = pdfs[dir] + omega_nu_ * ((sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdfs[dir]);
        }

      }
    }
  }
}
  return;
}

void LatticeBoltzmann::Prop_coll_velinlet3D(int cur, StdVector<StdVector<int> >& inlet)
{
  int x, y, z;
  double tmp_ux, tmp_uy, tmp_uz, tmp_us, sum;
  double tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];
    z = inlet[i][2];

    for (int dir = 0; dir < n_q_; dir ++)
      pdfs[dir] = PDF(cur,x,y,z,dir);

    sum = CalcDensity(pdfs);

    tmp_ux = ux_;
    tmp_uy = uy_;
    tmp_uz = uz_;

    if (!u_in.IsEmpty()) { //use parabolic profile
      tmp_ux = u_in[i][0];
      tmp_uy = u_in[i][1];
      tmp_uz = u_in[i][2];
    }
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    tmp_uz = 3.0 * tmp_uz;

    LOG_DBG3(lbm) << "pcv: i=" << i << " tux=" << tmp_ux << " tuy=" << tmp_uy << " tuz=" << tmp_uz << std::endl;


    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
      PDF(cur, x, y, z, dir)  = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


//
// Performs a bounce back step.
//
void LatticeBoltzmann::Prop_coll_bounce_back3D(int cur, StdVector<StdVector<int> >& bb)
{

  int x, y, z;
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];
    z = bb[i][2];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur, x, y, z, dir);
    }
    for (int  dir = 0; dir < n_q_; dir++) {
      PDF(cur, x, y, z, GetInvDirection((Direction)dir)) = pdfs[dir];
    }
  }
  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::Prop_coll_densoutlet3D(int cur, StdVector<StdVector<int> >& outlet)
{
  double tmp_ux, tmp_uy, tmp_uz, tmp_us, sum, tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  int x, y, z;

  for(unsigned int  i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];
    z = outlet[i][2];

    for (int dir = 0; dir < n_q_; dir ++)
      pdfs[dir] = PDF(cur,x,y,z,dir);

    CalcVelocities(pdfs, tmp_ux, tmp_uy, tmp_uz);

    sum = 1.0; // the enforced density
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    tmp_uz = 3.0 * tmp_uz;

    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
      PDF(cur, x, y, z, dir) = sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}

} // end namespace
