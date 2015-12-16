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

//  moments_.Resize(nNodes_);
//  eqMoments_.Resize(nNodes_);
  adjCollision.Resize(nNodes_);
  d_diss_d_m.Resize(nNodes_);
  d_pdrop_d_m.Resize(nNodes_);

  plot_ = plot;

  //matrix of the probability distributions
  LOG_DBG(lbm) << "Allocating arrays for " << nNodes_ * n_q_ << " PDFs (" << (sizeof(double) * nNodes_ * n_q_ * 2.0 / 1024.0 / 1024.0) << " MiB)";

  pdfs_.Resize(2);
  pdfs_[0].Resize(nNodes_ * n_q_);
  pdfs_[1].Resize(nNodes_ * n_q_);

  adjMoments_.Resize(2);
  adjMoments_[0].Resize(nNodes_ * n_q_);
  adjMoments_[1].Resize(nNodes_ * n_q_);
  adjMoments_[0].Init(0.0);
  adjMoments_[1].Init(0.0);

  adjPdfs_.Resize(2);
  adjPdfs_[0].Resize(nNodes_ * n_q_);
  adjPdfs_[1].Resize(nNodes_ * n_q_);
  adjPdfs_[0].Init(0.0);
  adjPdfs_[1].Init(0.0);

  // microVelDirections stores information about the 19 microscopic velocities/directions of D3Q19 model
  microVelDirections.Resize(n_q_);
  weights.Resize(n_q_);

  SetEnums();
  InitWeights();
  SetInvDirections();
  TestInvDirections();

  bounceback.Resize(n_q_);
  bounceback.Init();
  for (int dir = 0; dir < n_q_; dir++) {
    bounceback[dir][GetInvDirection((Direction)dir)] = 1.0;
  }

  scales.Resize(nNodes_);
  rel.Resize(nNodes_);
  bb.Resize(2 * sizeX_ + 2 * sizeY_ + 2 * sizeZ_);

  cur_  = 0;
  next_ = 1;
  adjCur_ = 0;
  adjNext_ = 1;

  lbmCalls_ = 0;

  SetMicroVelocities();

//  if (!srt_)
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
  double jx = 0;
  double jy = 0;
  double jz = 0;

  if (n_q_ == 9) {
    jx = (pdfs[Q_E] - pdfs[Q_W]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_SE] - pdfs[Q_NW]);
    jy = (pdfs[Q_N] - pdfs[Q_S]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_NW] - pdfs[Q_SE]);
  }
  else {
    jx = (pdfs[Q_E] - pdfs[Q_W]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_SE] - pdfs[Q_NW]) + (pdfs[Q_TE] - pdfs[Q_BW]) + (pdfs[Q_BE] - pdfs[Q_TW]);
    jy = (pdfs[Q_N] - pdfs[Q_S]) + (pdfs[Q_NW] - pdfs[Q_SE]) + (pdfs[Q_NE] - pdfs[Q_SW]) + (pdfs[Q_TN] - pdfs[Q_BS]) + (pdfs[Q_BN] - pdfs[Q_TS]);
    jz = (pdfs[Q_T] - pdfs[Q_B]) + (pdfs[Q_TW] - pdfs[Q_BE]) + (pdfs[Q_TE] - pdfs[Q_BW]) + (pdfs[Q_TN] - pdfs[Q_BS]) + (pdfs[Q_TS] - pdfs[Q_BN]);
  }

//  for (int  dir = 0; dir < n_q_; dir++) {
//    //store current pdf values in array for better accessing
//    jx += microVelDirections[dir].off_x * pdfs[dir];
//    jy += microVelDirections[dir].off_y * pdfs[dir];
//    jz += microVelDirections[dir].off_z * pdfs[dir];
//  }

  ux = jx / density;
  uy = jy / density;
  uz = jz / density;
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

//  LOG_DBG(lbm) << "bb = " << ToString(bb);
//  LOG_DBG(lbm) << "inlet = " << ToString(inlet);
//  LOG_DBG(lbm) << "outlet = " << ToString(outlet);
//  LOG_DBG(lbm) << "rel = " << ToString(rel);

  in->Get("converged")->SetValue("running");

  while(it < maxIter_ && !steady_state && R <= 1000)
  {
//    if (srt_)
//      CreateOutput("srt.pdfs",cur_);
//    else
//      CreateOutput("mrt.pdfs",cur_);
    LOG_DBG3(lbm) << "---------------------------Iteration " << it << "---------------------------------------------------";
    // -- Combined propagation and collision step -------------------------
    (this->*prop_coll_step)(cur_, next_);
    // -- Bounce back step ------------------------------------------------
    (this->*prop_coll_bounce_back)(next_);
    // -- Inlet condition -------------------------------------------------
    (this->*prop_coll_velinlet)(next_);
    // -- Outlet condition ------------------------------------------------
    (this->*prop_coll_densoutlet)(next_);

    if((it == 0 || it % 100 == 0)) // check convergence
    {
      R = CalcResidual(cur_,next_,false);

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

  PtrParamNode node = in->Get(ParamNode::PROCESS)->Get("forward", ParamNode::APPEND); // write out how many lbm iterations until convergence
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

  return &(pdfs_[cur_]);
}

void LatticeBoltzmann::InitializePdfs()
{
#pragma omp parallel for default(none)
  for (int elem = 0; elem < nNodes_; elem++) {
    for (int  dir = 0; dir < n_q_; dir++) {
      PDF(0, elem, dir) = weights[dir];
      PDF(1, elem, dir) = weights[dir];
      APDF(0, elem, dir) = weights[dir];
      APDF(1, elem, dir) = weights[dir];
    }
  }

  if (!srt_)
  {
    for (int elem = 0; elem < nNodes_; elem++) {
      Vector<double> pdfs(n_q_);
      Vector<double> moments(n_q_);
      for (int dir = 0; dir < n_q_; dir++)
        pdfs[dir] = PDF(0, elem, dir);
      transformation.Mult(pdfs, moments);
      for (int dir = 0; dir < n_q_; dir++) {
        AMoments(0, elem, dir) = moments[dir];
        AMoments(1, elem, dir) = moments[dir];
      }
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

  adjTransformation.Resize(n_q_);
  invTransformation.Transpose(adjTransformation); // adjoint transformation matrix is tranposed inverse of primal transformation matrix

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

  double porosity;

  for(int elem = 0; elem < nNodes_; elem++)
  {
    porosity = elements[elem];

    if (LbmNodeTypeIsFluid(porosity))
      rel.Push_back(elem);
    else if (LbmNodeTypeIsBB(porosity))
      bb.Push_back(elem);
    else if (LbmNodeTypeIsInlet(porosity))
      inlet.Push_back(elem);
    else if (LbmNodeTypeIsOutlet(porosity))
      outlet.Push_back(elem);
    else if (LbmNodeIsObstacle(porosity))
      obst.Push_back(elem);
  }

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

// Calculates macroscopic density for given pdfs of an element
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
  f1.Init();
  f2.Resize(n_q_);
  f2.Init();

  double alpha = CalcResistanceCoeff(elemId);

  double rho = moments[0];
  double ux = moments[3] / rho;
  double uy = moments[5] / rho;

  double fx = -alpha * rho * ux;
  double fy = -alpha * rho * uy;

  LOG_DBG3(lbm) << " CDF: Element " << elemId << " has LBM density " << rho << " and porosity " << scales[elemId];
  LOG_DBG3(lbm) << "CDF: alpha= " << alpha << " fx=" << fx << " fy=" << fy;
//  if (scales[elemId] == 0.0) {
//  for(int dir = 0; dir < n_q_; dir++)
//    LOG_DBG3(lbm) << pdfs[dir];
//  }

  f1[3] = fx; // the assembly of f1 and f2 are taken from MRT paper
  f1[4] = -fx;
  f1[5] = fy;
  f1[6] = -fy;

  f2[1] = 6.0 * (fx * ux + fy * uy);
  f2[2] = -6.0 * (fx * ux + fy * uy);
  f2[7] = 2.0 * (fx * ux - fy * uy);
  f2[8] = fy * ux + fx * uy;

  //debugging
//  if (scales[elemId] == 0.0) {
//    for(int dir = 0; dir < n_q_; dir++)
//      std::cout << pdfs[dir] << " ";
//  }

  LOG_DBG3(lbm) << "CDF: F1=" << f1.ToString(0,',');
  LOG_DBG3(lbm) << "CDF: F2=" << f2.ToString(0,',') << "\n";
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
  double rho = moments[0];
  double jx = moments[3];
  double jy = moments[5];

  //debugging
  double e_eq = -2.0 * rho + 3.0 / rho * (jx*jx + jy*jy);
  assert(std::fabs(e_eq - eqMoments[1]) < EPS);
  double pxx_eq = (jx*jx-jy*jy) / rho;
  if (std::fabs(pxx_eq - eqMoments[7]) > EPS)
  assert(std::fabs(pxx_eq - eqMoments[7]) < EPS);
  double pxy_eq = jx*jy/rho;
  assert(std::fabs(pxy_eq - eqMoments[8]) < EPS);


//  double term1 = (eqMoments[1] - moments[1]) * (eqMoments[1] - moments[1]); // (e^eq - e)^2
//  double term2 = (eqMoments[7] - moments[7]) * (eqMoments[7] - moments[7]); // (p_xx^eq - p_xx)^2
//  double term3 = (eqMoments[8] - moments[8]) * (eqMoments[8] - moments[8]); // (p_xy^eq - p_xy)^2
  double ux =  jx / rho;
  double uy =  jy / rho;
//
//  return nu * rho *  (0.25 * (omega_e_ * omega_e_ * term1 + 9 * omega_nu_ * omega_nu_ * term2) + 9.0 * omega_nu_ * omega_nu_ * term3 ) - fx * ux - fy * uy;
  double e = moments[1];
  double pxx = moments[7];
  double pxy = moments[8];

  double term1 = omega_e_ / (2.0 * rho) * (e_eq - e);
  double term2 = 3.0 * omega_nu_ / (2.0 * rho) * (pxx_eq - pxx);

  double dxux = 0.5 * (term1 + term2);
  double dyuy = 0.5 * (term1 - term2);
  double dxuydyux = 3.0 * omega_nu_ / rho * (pxy_eq - pxy);

//  return nu * rho * (2.0 * dxux*dxux + 2.0 * dyuy * dyuy + dxuydyux * dxuydyux) - fx * ux - fy * uy;
  double diss = nu * rho * (2.0 * dxux*dxux + 2.0 * dyuy * dyuy + dxuydyux * dxuydyux) - fx * ux - fy * uy;
  return jx * diss / diss;
}

void LatticeBoltzmann::CalcAdjointCollMatrix(int elemId, const Vector<double>& moments)
{
  // S_A = I - S(I - d_mEq/d_m) + d_F1/d_m + (I - S/2) d_F2/d_m

  double rho = moments[0];
//  double jx = moments[3] * scales[elemId];
//  double jy = moments[5] * scales[elemId];
  double jx = moments[3];
  double jy = moments[5];
  double ux = jx / rho;
  double uy = jy / rho;
  double u2 = ux * ux + uy * uy;

  Matrix<double> identity(n_q_,n_q_), relaxation(n_q_,n_q_);
  identity.InitValue(0.0);
  relaxation.InitValue(0.0);

  for (int dir = 0; dir < n_q_; dir++)
  {
    identity[dir][dir] = 1.0;
    relaxation[dir][dir] = relax_rates[dir];
  }

  // non-zero entries of d_mEq/d_m
  Matrix<double> d_mEq_d_m(n_q_,n_q_);
  d_mEq_d_m.InitValue(0.0);
//  d_mEq_d_m[0][0] = 1.0; // 0th row of d_mEq/d_m describes d_mEq/d_rho
//  d_mEq_d_m[1][0] = -2.0 - 3.0 * u2;
//  d_mEq_d_m[2][0] = 1.0 + 3.0 * u2;
//  d_mEq_d_m[7][0] = -(ux * ux - uy * uy);
//  d_mEq_d_m[8][0] = -ux * uy;
//
//  d_mEq_d_m[1][3] = 6.0 * ux; // 3th row of d_mEq/d_m describes d_mEq/d_jx
//  if(elemId == 4) {
//    std::cout << "d_mEq_d_m[1][3]=" << d_mEq_d_m[1][3] << "= 6.0 * " << ux << std::endl;
//  }
//  d_mEq_d_m[2][3] = -6.0 * ux;
//  d_mEq_d_m[3][3] = 1.0;
//  d_mEq_d_m[4][3] = -1.0;
//  d_mEq_d_m[7][3] = 2.0 * ux;
//  d_mEq_d_m[8][3] = uy;
//
//  d_mEq_d_m[1][5] = 6.0 * uy; // 5th row of d_mEq/d_m describes d_mEq/d_jy
//  d_mEq_d_m[2][5] = -6.0 * uy;
//  d_mEq_d_m[5][5] = 1.0;
//  d_mEq_d_m[6][5] = -1.0;
//  d_mEq_d_m[7][5] = -2.0 * uy;
//  d_mEq_d_m[8][5] = ux;

  double por = scales[elemId];
  double por2 = scales[elemId] * scales[elemId];
  d_mEq_d_m[0][0] = 1.0;
  d_mEq_d_m[1][0] = - 2.0 - 3.0 * por2 * u2;
  d_mEq_d_m[2][0] = 1.0 + 3.0 * por2 * u2;
  d_mEq_d_m[7][0] = - por2 * (ux * ux - uy * uy);
  d_mEq_d_m[8][0] = - por2 * ux * uy;

  d_mEq_d_m[1][3] = 6.0 * por2 * ux;
  d_mEq_d_m[2][3] = -6.0 * por2 * ux;
  d_mEq_d_m[3][3] = por;
  d_mEq_d_m[4][3] = -por;
  d_mEq_d_m[7][3] = 2.0 * por2 * ux;
  d_mEq_d_m[8][3] = por2 * uy;

  d_mEq_d_m[1][5] = 6.0 * por2 * uy;
  d_mEq_d_m[2][5] = -6.0 * por2 * uy;
  d_mEq_d_m[5][5] = por;
  d_mEq_d_m[6][5] = -por;
  d_mEq_d_m[7][5] = -2.0 * por2 * uy;
  d_mEq_d_m[8][5] = por2 * ux;


//  Matrix<double> d_F1_d_m(n_q_,n_q_);
//  d_F1_d_m.InitValue(0.0);
//  double alpha = CalcResistanceCoeff(elemId);
//  d_F1_d_m[3][3] = -alpha; // d_F1_d_jx
//  d_F1_d_m[4][3] = alpha;
//  d_F1_d_m[5][5] = -alpha; // d_F1_d_jy
//  d_F1_d_m[6][5] = alpha;
//
//  Matrix<double> d_F2_d_m(n_q_,n_q_);
//  d_F2_d_m.InitValue(0.0);
//  d_F2_d_m[1][0] = 6.0 * u2; // d_F2_d_rho; rho is LBM density
//  d_F2_d_m[2][0] = -6.0 * u2;
//  d_F2_d_m[7][0] = 2.0 * (ux * ux - uy * uy);
//  d_F2_d_m[8][0] = 2.0 * ux * uy;
//
//  d_F2_d_m[1][3] = -12.0 * ux; // d_F2/d_jx
//  d_F2_d_m[2][3] = 12.0 * ux;
//  d_F2_d_m[7][3] = -4.0 * ux;
//  d_F2_d_m[8][3] = -2.0 * uy;
//
//  d_F2_d_m[1][5] = -12.0 * uy; // d_F2/d_jy
//  d_F2_d_m[2][5] = 12.0 * uy;
//  d_F2_d_m[7][5] = 4.0 * uy;
//  d_F2_d_m[8][5] = -2.0 * ux;
//  d_F2_d_m *= alpha;

//  if (elemId == 4) {
//    std::cout << "\n Element 4\n d_mEq_d_m\n:";
//    for (int i = 0; i < n_q_; i++) {
//      for (int j = 0; j < n_q_; j++)
//        if (d_mEq_d_m[i][j] != 0)
//          std::cout << "(" << i << "," << j << "): " << d_mEq_d_m[i][j] << std::endl;
//    }
//    std::cout << "\n d_F1_d_m:\n";
//    for (int i = 0; i < n_q_; i++) {
//      for (int j = 0; j < n_q_; j++)
//        if (d_F1_d_m[i][j] != 0)
//          std::cout << "(" << i << "," << j << "): " << d_F1_d_m[i][j] << std::endl;
//    }
//    std::cout << "\n d_F2_d_m:\n";
//    for (int i = 0; i < n_q_; i++) {
//      for (int j = 0; j < n_q_; j++)
//        if (d_F2_d_m[i][j] != 0)
//          std::cout << "(" << i << "," << j << "): " << d_F2_d_m[i][j] << std::endl;
//    }
//  }

//  Matrix<double> relax_tmp(relaxation);
//  relax_tmp *= 0.5; // matrix is diagonal

  if (inlet.Contains(elemId))
  {
    Matrix<double>& mat = adjCollision[elemId];
    mat.Resize(n_q_);
    mat.Init();
    mat[0][0] = 1.0; // 0th column: d_mEq/d_rho
    mat[1][0] = -2.0 + 3.0 * u2;
    mat[2][0] = 1.0 - 3.0 * u2;
    mat[3][0] = ux;
    mat[4][0] = -ux;
    mat[5][0] = uy;
    mat[6][0] = -uy;
    mat[7][0] = ux * ux - uy * uy;
    mat[8][0] = ux * uy;
  }
  else if (outlet.Contains(elemId))
  {
    Matrix<double>& mat = adjCollision[elemId];
    mat.Resize(n_q_);
    mat.Init();
    mat[1][3] = 6.0 * ux;  // 3rd column: d_mEq_d_jx
    mat[2][3] = -6.0 * ux;
    mat[3][3] = 1.0;
    mat[4][3] = -1.0;
    mat[7][3] = 2.0 * ux;
    mat[8][3] = uy;

    mat[1][5] = 6.0 * uy; // 5th column: d_mEq_d_jy
    mat[2][5] = -6.0 * uy;
    mat[5][5] = 1.0;
    mat[6][5] = -1.0;
    mat[7][5] = -2.0 * uy;
    mat[8][5] = ux;
  }
  else if (bb.Contains(elemId)) // adjoint of node bounce back is node bounce back
  {
    Matrix<double>& mat = adjCollision[elemId]; // bounce back operator in moment space
    mat.Resize(n_q_);
    mat.Init();
//    for (int dir = 0; dir < n_q_; dir++)
//      mat[0][dir] = 1.0;
//    mat[1][0] = -4;
//    mat[1][1] = -1;
//    mat[1][2] = -1;
//    mat[1][3] = -1;
//    mat[1][4] = -1;
//    mat[1][5] = 2;
//    mat[1][6] = 2;
//    mat[1][7] = 2;
//    mat[1][8] = 2;
//    mat[2][0] = 4;
//    mat[2][1] = -2;
//    mat[2][2] = -2;
//    mat[2][3] = -2;
//    mat[2][4] = -2;
//    mat[2][5] = 1;
//    mat[2][6] = 1;
//    mat[2][7] = 1;
//    mat[2][8] = 1;
//    mat[3][1] = -1;
//    mat[3][3] = 1;
//    mat[3][5] = -1;
//    mat[3][6] = 1;
//    mat[3][7] = 1;
//    mat[3][8] = -1;
//    mat[4][1] = 2;
//    mat[4][3] = -2;
//    mat[4][5] = -1;
//    mat[4][6] = 1;
//    mat[4][7] = 1;
//    mat[4][8] = -1;
//    mat[5][2] = -1;
//    mat[5][4] = 1;
//    mat[5][5] = -1;
//    mat[5][6] = -1;
//    mat[5][7] = 1;
//    mat[5][8] = 1;
//    mat[6][2] = 2;
//    mat[6][4] = -2;
//    mat[6][5] = -1;
//    mat[6][6] = -1;
//    mat[6][7] = 1;
//    mat[6][8] = 1;
//    mat[7][1] = 1;
//    mat[7][2] = -1;
//    mat[7][3] = 1;
//    mat[7][4] = -1;
//    mat[8][5] = 1;
//    mat[8][6] = -1;
//    mat[8][7] = 1;
//    mat[8][8] = -1;
//    mat[0][0] = 1.0;
//    mat[1][1] = 1.0;
//    mat[2][2] = 1.0;
//    mat[3][3] = -1.0;
//    mat[4][4] = -1.0;
//    mat[5][5] = -1.0;
//    mat[6][6] = -1.0;
//    mat[7][7] = 1.0;
//    mat[8][8] = 1.0;
    mat[0][0] = 1.0;
    mat[1][1] = 1.0;
    mat[2][2] = 1.0;
    mat[3][3] = 1.0;
    mat[4][4] = 1.0;
    mat[5][5] = 1.0;
    mat[6][6] = 1.0;
    mat[7][7] = 1.0;
    mat[8][8] = 1.0;
  } else
  {
    Matrix<double> dcorr_dm(n_q_,n_q_);
    dcorr_dm.Init();
    double s = scales[elemId];
    double rho2 = rho * rho;

    dcorr_dm[0][0] = -6.938893904e-18*jx*jy*omega_nu_*(s*s - 1.387778781e-17*s + 1.0)/rho2;
    dcorr_dm[1][0] = -dcorr_dm[0][0];
    dcorr_dm[2][0] = 1.387778781e-17*jx*jy*omega_nu_*(s*s - 1.387778781e-17*s + 1.0)/rho2;
    dcorr_dm[5][0] = dcorr_dm[0][0];
    dcorr_dm[6][0] = dcorr_dm[2][0];
    dcorr_dm[7][0] = -dcorr_dm[0][0];

    dcorr_dm[0][3] = -7.703719778e-34*omega_nu_*(1.651319863e16*rho - 9.007199255e15*jy + 1.801439851e16*jy*s + 1.501199876e15*rho*s - 9.007199255e15*jy*s*s)/rho;
    dcorr_dm[1][3] = 7.703719778e-34*omega_nu_*(4.803839603e16*rho - 9.007199255e15*jy + 1.801439851e16*jy*s + 6.004799503e15*rho*s - 9.007199255e15*jy*s*s)/rho;
    dcorr_dm[2][3] = -1.540743956e-33*omega_nu_*(9.007199255e15*jy + 1.501199876e16*rho - 1.801439851e16*jy*s + 3.002399752e15*rho*s + 9.007199255e15*jy*s*s)/rho;
    dcorr_dm[3][3] = -omega_nu_*s;
    dcorr_dm[5][3] = 7.703719778e-34*omega_nu_*(9.007199255e15*jy - 4.203359652e16*rho - 1.801439851e16*jy*s + 3.602879702e16*rho*s + 9.007199255e15*jy*s*s)/rho;
    dcorr_dm[6][3] = -1.540743956e-33*omega_nu_*(9.007199255e15*jy - 4.203359652e16*rho - 1.801439851e16*jy*s + 3.602879702e16*rho*s + 9.007199255e15*jy*s*s)/rho;
    dcorr_dm[7][3] = -7.703719778e-34*omega_nu_*(9.007199255e15*jy*s*s - 1.801439851e16*jy*s + 9.007199255e15*jy - 6.004799503e15*rho)/rho;

    dcorr_dm[0][5] = 7.703719778e-34*omega_nu_*(9.007199255e15*jx + 1.050839913e16*rho - 1.801439851e16*jx*s + 1.501199876e15*rho*s + 9.007199255e15*jx*s*s)/rho;
    dcorr_dm[1][5] = -3.081487911e-33*omega_nu_*(2.251799814e15*jx + 1.050839913e16*rho - 4.503599627e15*jx*s + 1.501199876e15*rho*s + 2.251799814e15*jx*s*s)/rho;
    dcorr_dm[2][5] = 3.081487911e-33*omega_nu_*(1.050839913e16*rho - 4.503599627e15*jx + 9.007199255e15*jx*s + 1.501199876e15*rho*s - 4.503599627e15*jx*s*s)/rho;
    dcorr_dm[5][5] = 6.938893904e-18*omega_nu_*(jx - 2.0*jx*s - 1.441151881e17*rho*s + jx*s*s)/rho;
    dcorr_dm[6][5] = -1.387778781e-17*jx*omega_nu_*(s - 1.0)*(s - 1.0)/rho;
    dcorr_dm[7][5] = -6.938893904e-18*jx*omega_nu_*(s - 1.0)*(s - 1.0)/rho;

    Matrix<double> tmp1 = identity;
    tmp1.Add(-1.0,d_mEq_d_m);
//    Matrix<double> tmp2 = identity;
//    tmp2.Add(-1.0,relax_tmp);
    // S_A = I - S(I - d_mEq/d_m) + d_F1/d_m + (I - S/2) d_F2/d_m
//    adjCollision[elemId] = identity - relaxation * tmp1 + d_F1_d_m + tmp2 * d_F2_d_m;
//    adjCollision[elemId] = identity - relaxation * tmp1;
    adjCollision[elemId] = identity - relaxation * tmp1 + dcorr_dm; // I - S(I-dmEq/dm) + dcorr/dm
//    std::cout << "identity:\n" << identity.ToString(0,true) << "\n";
//    std::cout << "relaxation: \n" << relaxation.ToString(0,true) << "\n";
//    std::cout << "dmEq_dm:\n" << d_mEq_d_m.ToString(0,true) << "\n";
//    std::cout << "dF1_dm:\n" << d_F1_d_m.ToString(0,true) << "\n";
//    std::cout << "0.5*s:\n" << relax_tmp.ToString(0,true) << "\n";
//    std::cout << "dF2_dm:\n" << d_F2_d_m.ToString(0,true) << "\n";
//    std::cout << "I - dmEq_dm:\n" << tmp1.ToString(0,true) << "\n";
//    std::cout << "I - S/2:\n" << tmp2.ToString(0,true) << "\n";
//    std::cout << "(I - S/2) d_F2/d_m: \n" << (tmp2*d_F2_d_m).ToString(0,true) << "\n";
//    std::cout << "S/2:\n" << relax_tmp.ToString(0,true) << "\n";
//    std::cout << "I-S(I - d_mEq/d_m):\n" << (identity-relaxation*tmp1).ToString(0,true) << "\n";
//    Matrix<double> out1 = relaxation*tmp1;
//    Matrix<double> out2 = identity - out1;
//    Matrix<double> out3 = out2 + d_F1_d_m;
//    std::cout << "Version 2: I-S(I - d_mEq/d_m):\n" << out2.ToString(0,true) << "\n";
//    std::cout << "dF1_dm:\n" << d_F1_d_m.ToString(0,true) << "\n";
//    std::cout << "I-S(I - d_mEq/d_m) + d_F1/d_m:\n" << (identity-relaxation*tmp1+d_F1_d_m).ToString(0,true) << "\n";
//    std::cout << " Version 2 I-S(I - d_mEq/d_m) + d_F1/d_m:\n" << out3.ToString(0,true) << "\n";
//    std::cout << "Adjoint collision matrix of element " << elemId << "\n" << adjCollision[elemId].ToString(0,true) << std::endl;
  }
}

void LatticeBoltzmann::d_diss_d_moments(int elemId, const Vector<double>& moments)
{
  Vector<double> result(n_q_);
  result.Init();

  if (!rel.Contains(elemId)) { // at the boundary the derivative of dissipation does not exist
    d_diss_d_m[elemId] = result;
    return;
  }

  double rho = moments[0];
  double e = moments[1];
  double jx = moments[3];
  double jy = moments[5];
  double pxx = moments[7];
  double pxy = moments[8];

  double nu = (1/omega_nu_ - 0.5) / 3.0; // fluid's kinematic viscosity in LBM units
  double se2 = omega_e_ * omega_e_;
  double snu2 = omega_nu_ * omega_nu_;
  double j2 = jx * jx + jy * jy;
  double rho2 = rho * rho;
  double rho3 = rho * rho * rho;
  double rho4 = rho * rho * rho * rho;
  double alpha = CalcResistanceCoeff(elemId);
  double fx = -alpha * jx;
  double fy = -alpha * jy;

  // 0th entry: d_diss/d_rho
  result[0] = 4.0 * nu * se2 - 27.0 * nu * (4.0 * se2 + snu2) * j2 * j2 / (4.0 * rho4)
      + 3.0 * nu * (4.0 * se2 * e * j2 + 3.0 * snu2 * (pxx * (jx * jx - jy * jy) + 4.0 * jx * jy * pxy)) / rho3
      - nu * (4.0 * se2 * (e * e - 12.0 * j2) + 9.0 * snu2 * (pxx * pxx + 4.0 * pxy * pxy)) / (4.0 * rho2) + (fx * jx + fy * jy) / rho2;
//  result[0] = -9.0 * nu * snu2 * pxy * pxy / rho2 - 9.0 * nu * snu2 * pxx * pxx / (4*rho2)+ 36.0 * nu * snu2 * jx * jy * pxy / rho3
//      - 9 * nu * snu2 * pxx * jy2 / rho3 + 9.0 * nu * snu2 * pxx * jx * jx / rho3
//      - 27.0 * nu * snu2 * jy * jy * jy * jy / (4*rho4) - 27 * nu * snu2 * jx *jx * jy * jy / (2.0*rho4) - 27.0 * nu * snu2 * jx * jx * jx * jx / (4*rho4)
//      + nu * se2 + 3.0 * nu * se2 / rho2 * jy2 + 3 * nu * se2 / rho2 * jx * jx2 - nu * se2 * e * e / (4.0*rho2) + 3 * nu * se2 * e * jy * jy / rho3
//      + 3 * nu * se2 * jx * jx * e / rho3 - 27 * nu * se2 * jy * jy * jy * jy / (4*rho4) - 27 * nu * se2 * jx * jx * jy * jy / (2*rho4)
//      - 27 * nu * se2 * jx * jx * jx * jx / (4*rho4) - alpha * jy * jy / rho2 - alpha * jx * jx / rho2;

  // 1th entry: d_diss/d_e
  result[1] = 2.0 * nu * se2 * (-3.0 * j2 + rho * (e + 2.0 * rho)) / rho2;
//  result[1] = nu * se2 / (2.0 * rho) * e - (3.0 * nu * se2 * jy2 + 3.0 * nu * se2 * jx2) / (2 * rho2) + nu * se2;

  // 3th entry: d_diss/d_jx
  result[3] = 9.0 * nu * jx * (4.0 * se2 + snu2) * j2 / rho3 - 3.0 * nu * (4.0 * se2 * e * jx + 3.0 * snu2 * (jx * pxx + 2.0 * jy * pxy)) / rho2
      - (fx + 24.0 * nu * se2 * jx - jx * alpha) / rho;
//  result[3] = (9.0 * nu * snu2 * jx * jy2 + 9.0 * nu * snu2 * jx2 * jx + 9.0 * nu * se2 * jx * jy2 + 9.0 * nu * se2 * jx2 * jx) / rho3
//      + (-18.0 * nu * snu2 * jy * pxy - 9.0 * nu * snu2 * jx * pxx - 3.0 * nu * se2 * e * jx) / rho2
//      + (-6.0 * nu * se2 * jx + 2.0 * alpha * jx) / rho;

  // 5th entry: d_diss/d_jy
  result[5] = 9.0 * nu * jy * (4.0 * se2 + snu2) * j2 / rho3 - 3.0 * nu * (4.0 * se2 * e * jy + 3.0 * snu2 *(-jy * pxx + 2.0 * jx * pxy)) / rho2
      - (fy + 24.0 * nu * se2 * jy - alpha * jy) / rho;
//    result[5] = (-18.0 * nu * snu2 * jx * pxy + 9.0 * nu * snu2 * jy * pxx - 3.0 * nu * se2 * e * jy) / rho2
//        + (9.0 * nu * snu2 * jy2 * jy + 9.0 * nu * snu2 * jx2 * jy + 9.0 * nu * se2 * jy2 * jy + 9.0 * nu * se2 * jx2 * jy) / rho3
//        + (-6.0 * nu * se2 * jy + 2.0 * alpha * jy) / rho;

  // 7th entry: d_diss/d_pxx
  result[7] = 9.0 * nu * snu2 * (-jx * jx + jy * jy + pxx * rho) / (2.0 * rho2);

  // 8th entry: d_diss/d_pxy
  result[8] = 18.0 * nu * snu2 * (-jx * jy + pxy * rho) / rho2;

  result.Init();
  result[3] = 1.0;
  d_diss_d_m[elemId] = result;
}

void LatticeBoltzmann::d_pdrop_d_moments(int elemId, const Vector<double>& moments)
{
  Vector<double> result(n_q_);
  result.Init();

  if (inlet.Contains(elemId)) {
    result[0] = 1.0 / 3.0 + 0.5 * (ux_ * ux_ + uy_ * uy_);
    result /= (double) inlet.GetSize();
  }
  else if (outlet.Contains(elemId)) {
    double rho = 1.0; // enforced density at outlet nodes
    result[3] = -rho * moments[3];
    result[5] = -rho * moments[5];
    result /= (double) outlet.GetSize();
  }
  d_pdrop_d_m[elemId] = result;
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
    Vector<double> pdfs(n_q_);
    Vector<double> eqMoments(n_q_);
    Vector<double> moments(n_q_);
#pragma omp for collapse(2)
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

          Vector<double> res(n_q_), noneq(n_q_);
          // propagation and collision in one step
          for (int  dir = 0; dir < n_q_; dir++)
          {
            tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy ;
            // no collision on the boundaries
            if (x == 0 || y == 0 || x == sizeX_ - 1 || y == sizeY_ - 1)
              PDF(next, x, y, z, dir) = pdfs[dir];
            else {
              PDF(next, x, y, z, dir) = pdfs[dir] + omega_nu_ * (sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us) - pdfs[dir]);
              res[dir] = PDF(next, x, y, z, dir);
              noneq[dir] = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us) - pdfs[dir];
            }
          }
//          if (index == 4) {
//            std::cout << "SRT noneq PDFS:\n" << (-noneq).ToString(2,' ') << std::endl;
//            noneq *= -omega_nu_;
//            std::cout << "SRT pdfs of elem 4:\n" << pdfs.ToString(2,' ') << std::endl;
//            std::cout << "SRT subtrahend:\n" << noneq.ToString(2,' ') << std::endl;
//            std::cout << "SRT result for element 4:\n" << res.ToString(2,' ') << std::endl;
//          }
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

// using porosity model of Pingen [u*=(1-s)u]
	  
          double jx = moments[3];
          moments[3] *= scales[index];
          double jy = moments[5];
          moments[5] *= scales[index];
          CalcEquilMoments(moments, eqMoments); // compute equilibirum moments from moments

          Vector<double> momentsAfterCollision; // result of collision step in moment space including porosity model
//          Vector<double> term1(n_q_);
//          Vector<double> f1, f2;
//          CalcDarcyForce(moments,index,f1,f2);
          moments[3] = jx;
          moments[5] = jy;
          Vector<double> noneq_moments(n_q_);
          noneq_moments = moments - eqMoments;

          // S * (m - m_eq)
//          for (int dir = 0; dir < n_q_; dir++) {
//            term1[dir] = relax_rates[dir] * noneq_moments[dir];
//          }
          Matrix<double> relaxation(n_q_,n_q_);
          relaxation.Init();
          for (int dir = 0; dir < n_q_; dir++)
            relaxation[dir][dir] = relax_rates[dir];
//          term1 = relaxation * noneq_moments; // S * (m - mEq)

//          Vector<double> term2(n_q_);
//          for (int dir = 0; dir < n_q_; dir++) {
//            term2[dir] = (1.0 - 0.5 * relax_rates[dir]) * f2[dir];
//          }

          // m* = m - S * (m - m_eq) + F1 + (I - S/2) * F2
//          momentsAfterCollision = moments - term1 + f1 + term2;
          Vector<double> correction(n_q_); // correction necessary as we insert the porosity model directly into momentum space
          correction.Init();
          double s = 1-scales[index];
          double rho = moments[0];
//          correction[1] = -jx;
//          correction[3] = jx;
//          correction[5] = -jx-jy;
//          correction[6] = jx-jy;
//          correction[7] = jx+jy;
//          correction[8] = -jx+jy;
//          correction *= omega_nu_ * s / 6.0;
//          correction[0] = omega_nu_ * (jx-jy) * (-8.095376221e-18  - 1.156482317e-18*s);
//          correction[2] = omega_nu_ * (1.387778781e-17*jx*(s - 1) - jy*s/6.0 + 6.938893904e-18*jx*jy/rho *(1 + s*s - 1.387778781e-17*s) );
//          correction[4] =  omega_nu_* (1.387778781e-17*jx*(1-s) +jy*s/6.0);
          correction[0] = -omega_nu_*((jx - jy)*(1.156482317e-18*s + 8.095376221e-18) - 6.938893904e-18*jx*jy*(s*s - 1.387778781e-17*s + 1.0)/rho);
          correction[1] = omega_nu_*((1.156482317e-18*s + 8.095376221e-18)*4.0*(jx - jy) - 6.938893904e-18*jx*jy*(s*s - 1.387778781e-17*s + 1.0)/rho);
          correction[2] = -2.0*omega_nu_*((1.156482317e-18*s + 8.095376221e-18)*2.0*(jx - jy) + (6.938893904e-18*jx*jy*(s*s - 1.387778781e-17*s + 1.0))/rho);
          correction[3] = -jx*omega_nu_*s;
          correction[4] = 0.0;
          correction[5] = jx*omega_nu_*2.775557562e-17*(s -1) - 1.0*jy*omega_nu_*s + (6.938893904e-18*jx*jy*omega_nu_*(s*s - 1.387778781e-17*s + 1.0))/rho;
          correction[6] = 1.925929945e-34*jx*jy*omega_nu_*s/rho - 1.387778781e-17*jx*jy*omega_nu_*(s*s + 1.0)/rho - 5.551115124e-17*jx*omega_nu_*(s - 1.0);
          correction[7] = -jx*jy*omega_nu_*(6.938893904e-18*s*s - 9.629649724e-35*s + 6.938893904e-18)/rho;
          correction[8] = 0.0;
          Vector<double> subtrahend(n_q_), tmp(n_q_);
          tmp = relaxation * noneq_moments;
          tmp.Add(-1,correction);
          invTransformation.Mult(tmp, subtrahend);
//          subtrahend -= correction;
//          momentsAfterCollision = moments - (relaxation * noneq_moments); // m* = m - S * (m - m_eq)
//          collResult.Resize(n_q_);
          collResult = pdfs - subtrahend;
//          invTransformation.Mult(momentsAfterCollision,collResult);
          for (int  dir = 0; dir < n_q_; dir++)
            PDF(next, x, y, z, dir) = collResult[dir];

//          if (index == 4) {
//            Vector<double> out(n_q_);
//            invTransformation.Mult(noneq_moments,out);
//            std::cout << "MRT noneq PDFS:\n" << out.ToString(2,' ') << std::endl;
//            std::cout << "MRT pdfs of elem 4:\n" << pdfs.ToString(2,' ') << std::endl;
//            std::cout << "MRT subtrahend:\n" << subtrahend.ToString(2,' ') << std::endl;
//            std::cout << "MRT result for element 4:\n" << collResult.ToString(2,' ') << std::endl;
//          }
        }

      }
    }
  }
  return;
}



void LatticeBoltzmann::Prop_coll_velinlet2D(int cur)
{
  assert(uz_ == 0);

  double tmp_ux, tmp_uy, tmp_us, sum;
  double tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < inlet.GetSize(); i++) {
    int index = inlet[i];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur,index,dir);
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
      PDF(cur, index, dir)  = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }
  return;
}
//
// Performs a bounce back step.
//
void LatticeBoltzmann::Prop_coll_bounce_back2D(int cur)
{
  Vector<double> pdfs(n_q_), res(n_q_);

  for(unsigned int  i = 0; i < bb.GetSize(); i++) {
    int index = bb[i];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur, index, dir);
    }
    bounceback.Mult(pdfs,res);
    for (int  dir = 0; dir < n_q_; dir++) {
      PDF(cur, index, dir) = res[dir];
//      PDF(cur, index, GetInvDirection((Direction)dir)) = pdfs[dir];
    }
  }

  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::Prop_coll_densoutlet2D(int cur)
{
  double tmp_ux, tmp_uy, tmp_us, sum, tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < outlet.GetSize(); i++) {
    int index = outlet[i];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur,index,dir);
    }

    CalcVelocities(pdfs, tmp_ux, tmp_uy, tmp);

    sum = 1.0; // the enforced density

    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy;
      PDF(cur, index, dir) =  sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }
  return;
}

/************************************************** adjoint 2D operators *****************************************************/
void LatticeBoltzmann::AdjointCollision(int cur)
{
  tmpPdfs_.Resize(nNodes_ * n_q_);
  int index;
  int z = 0;
  Vector<double> moments(n_q_);
  for (int x = 0; x < sizeX_ ; x++)
    for (int y = 0; y < sizeY_; y++)
    {
      index = GetIndex(x,y,z);

      for (int dir = 0; dir < n_q_; dir++)
        moments[dir] = AMoments(cur,index,dir);

      Vector<double> momentsAfterCollision(n_q_); // result of collision step in moment space including porosity model
      Matrix<double> collMatrix = adjCollision[index];
      Matrix<double> collMatrixT(n_q_,n_q_);
      collMatrix.Transpose(collMatrixT);
//      momentsAfterCollision = d_diss_d_m[index] + collMatrixT * moments;
      momentsAfterCollision = d_pdrop_d_m[index] + collMatrixT * moments;

      Vector<double> collResult(n_q_);
      Matrix<double> transpose(n_q_, n_q_);
      transformation.Transpose(transpose); // adjoint backstransformation matrix is tranpose of primal transformation matrix
      transpose.Mult(momentsAfterCollision, collResult);

      for (int dir = 0; dir < n_q_; dir++)
        tmpPdfs_[GetPdfIndex(index,dir)] = collResult[dir];
    }
}

void LatticeBoltzmann::AdjointBounceBack(int cur)
{
  Vector<double> pdfs(n_q_), res(n_q_);

  for(unsigned int  i = 0; i < bb.GetSize(); i++) {
    int index = bb[i];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = tmpPdfs_[GetPdfIndex(index, dir)];
    }
    bounceback.Mult(pdfs,res);
    for (int  dir = 0; dir < n_q_; dir++) {
      tmpPdfs_[GetPdfIndex(index, dir)] = res[dir];
    }
  }
}

void LatticeBoltzmann::AdjointPropagation(int cur, int next)
{
  Vector<double> pdfs(n_q_);

  int z = 0;
  int tmp_x, tmp_y, tmp_z = 0;
  Vector<double> tmp(n_q_);
  for (int x = 0; x < sizeX_; x++)
    for (int y = 0; y < sizeY_; y++)
    {
      // propagation
      for (int  dir = 0; dir < n_q_; dir++) {
        if (PointsToBoundary(x,y,z,dir)) { // if the neighbor element that I want to access is outside the domain, keep current value
          tmp_x = x; // here we only set the coordinates
          tmp_y = y;
        }
        // else: standard propagation (get value from neighbor pdf)
        else {
//          tmp_x = microVelDirections[dir].off_x + x;
//          tmp_y = microVelDirections[dir].off_y + y;
          tmp_x = microVelDirections[dir].off_x + x;
          tmp_y = microVelDirections[dir].off_y + y;
        }

        int index = GetIndex(tmp_x,tmp_y,tmp_z);

        tmp[dir] = tmpPdfs_[GetPdfIndex(index,dir)];
      }

      if (srt_)
      {
        for (int dir = 0; dir < n_q_; dir++)
          APDF(next,x,y,z,dir) = tmp[dir];
      } else {
        Vector<double> propResult(n_q_);
        adjTransformation.Mult(tmp,propResult); // transforming propagation result (adjoint distributions) back to moment space
        for (int dir = 0; dir < n_q_; dir++)
          AMoments(next,x,y,z,dir) = propResult[dir];
      }
    }
}

StdVector<double>* LatticeBoltzmann::IterateAdjoint(PtrParamNode info)
{
  Vector<double> pdfs(n_q_);
  double dissipation = GetDissipation();
//  std::cout << "Dissipation: " << dissipation << std::endl;
  for (int elem = 0; elem < nNodes_; elem++) {
    for (int dir = 0; dir < n_q_; dir++)
      pdfs[dir] = PDF(cur_,elem,dir);

    Vector<double> moms(n_q_);
    transformation.Mult(pdfs,moms);
//    Vector<double> eqMoms;
//    CalcEquilMoments(moms,eqMoms);

//    if (rel.Contains(elem))
//    {// we assume that dissipation does not happen on the boundary
//      Vector<double> f1,f2;
//      CalcDarcyForce(moms,elem,f1,f2);
//      dissipation += CalcDissipation(moms,eqMoms,f1[3],f1[5]);
//    }

    CalcAdjointCollMatrix(elem,moms);
    d_diss_d_moments(elem,moms);
    d_pdrop_d_moments(elem,moms);
  }

  int count = numWriteResults_;

  Timer timer;
  timer.Start();

  int it = 0;
  double R = 0.0;
  bool steady_state = false;

  std::ofstream plot;
  if(plot_)
    plot.open(std::string(progOpts->GetSimName() + ".adjLbm.dat").c_str());

  LOG_DBG3(lbm) << "\n steady state pdfs: " << pdfs_.ToString(false) << std::endl;

  while(it < maxIter_ && !steady_state && R <= 1000)
  {
    LOG_DBG3(lbm) << "---------------------------Adjoint Iteration " << it << "---------------------------------------------------";
    AdjointCollision(adjCur_);
    AdjointBounceBack(adjCur_);
    AdjointPropagation(adjCur_,adjNext_);

    if((it == 0 || it % 100 == 0))
    {
      R = CalcResidual(adjCur_,adjNext_,true);

      if(R <= maxTol_)
        steady_state = true;

      if(plot_)
        plot << it << "\t" << R << "\n";

      domain->GetInfoRoot()->ToFile(); // is not written when called too often
    }

    adjCur_  = (adjCur_  + 1) % 2;
    adjNext_ = (adjNext_ + 1) % 2;

    it++;

    if (writeIntermediateResults_) {
      if (it % writeFrequency_ == 0) {
        domain->GetDriver()->StoreResults(count,(double) numIterations_ + it);
        count++;
      }
    }
  }

  if(R >= 1000)
    EXCEPTION("In adjoint LBM iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("internal LBM simulation could not converge: iterations: " << it << " residuum: " << R);

  numWriteResults_ = count;

  timer.Stop();

  PtrParamNode node = info->Get(ParamNode::PROCESS)->Get("adjoint", ParamNode::APPEND); // write out how many lbm iterations until convergence
  node->Get("number")->SetValue(lbmCalls_-1);
  node->Get("iterations")->SetValue(it);
  node->Get("residuum")->SetValue(R);
  node->Get("converged")->SetValue(steady_state);
  node->Get("dissipation")->SetValue(dissipation);

  return &(adjMoments_[adjCur_]);
}

StdVector<double>* LatticeBoltzmann::IterateAdjointSRT(const StdVector<Matrix<double> >& collisionMatrices, const StdVector<Vector<double> >& d_pdrop_d_f)
{
  tmpPdfs_.Resize(nNodes_ * n_q_);
  int z = 0;

  int count = numWriteResults_;

  int it = 0;
  double R = 0.0;
  bool steady_state = false;

  InitializePdfs();

  while(it < maxIter_ && !steady_state && R <= 1000)
  {
    /***************** Adjoint SRT collision ***/
    Vector<double> pdfs(n_q_);
    for (int x = 0; x < sizeX_ ; x++)
      for (int y = 0; y < sizeY_; y++)
      {
        int index = GetIndex(x,y,z);

        for (int dir = 0; dir < n_q_; dir++)
          pdfs[dir] = APDF(adjCur_,index,dir);

        // adjoint collision: f* = d_pdrop_d_f + (d_coll_d_f)^T * f
        Matrix<double> d_coll_d_f = collisionMatrices[index]; // collision matrices are already transposed
        Vector<double> d_pd_d_f = d_pdrop_d_f[index];
        Vector<double> collResult(n_q_), tmp(n_q_);
        d_coll_d_f.Mult(pdfs,tmp);
        for (int dir = 0; dir < n_q_; dir++)
          collResult[dir] = d_pd_d_f[dir] + tmp[dir];

        for (int dir = 0; dir < n_q_; dir++)
          tmpPdfs_[GetPdfIndex(index,dir)] = collResult[dir];

        if (!rel.Contains(index)) {
          for (int dir1 = 0; dir1 < n_q_; dir1++) {
//            for (int dir2 = 0; dir2 < n_q_; dir2++)
//              if (inlet.Contains(index) || outlet.Contains(index))
//                assert(d_coll_d_f[dir1][dir2] == 0.0);
            if (bb.Contains(index) || rel.Contains(index))
              assert(d_pd_d_f[dir1] == 0.0);
          }
        }
      }

//    AdjointBounceBack(adjCur_);
    AdjointPropagation(adjCur_, adjNext_);

    if((it == 0 || it % 100 == 0))
    {
      R = CalcResidual(adjCur_,adjNext_,true);
      if(R <= maxTol_)
        steady_state = true;
    }

    adjCur_  = (adjCur_  + 1) % 2;
    adjNext_ = (adjNext_ + 1) % 2;

    it++;
    if (writeIntermediateResults_) {
      if (it % writeFrequency_ == 0) {
        domain->GetDriver()->StoreResults(count,(double) numIterations_ + it);
        count++;
      }
    }
  }

  adjCur_  = (adjCur_  + 1) % 2;
  adjNext_ = (adjNext_ + 1) % 2;

//  std::cout << "Adjoint SRT simulation reached steady state after " << it << " iterations" << std::endl;
  if(R >= 1000)
    EXCEPTION("In adjoint SRT iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("Adjoint SRT simulation could not converge: iterations: " << it << " residuum: " << R);

  numWriteResults_ = count;

  return &(adjPdfs_[adjCur_]);
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
            PDF(next, index, dir) = pdfs[dir];
          else
            PDF(next, index, dir) = pdfs[dir] + omega_nu_ * ((sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdfs[dir]);
        }

      }
    }
  }
}
  return;
}

void LatticeBoltzmann::Prop_coll_velinlet3D(int cur)
{
  double tmp_ux, tmp_uy, tmp_uz, tmp_us, sum;
  double tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < inlet.GetSize(); i++) {
    int index = inlet[i];
    for (int dir = 0; dir < n_q_; dir ++)
      pdfs[dir] = PDF(cur,index,dir);

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
      PDF(cur, index, dir)  = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


//
// Performs a bounce back step.
//
void LatticeBoltzmann::Prop_coll_bounce_back3D(int cur)
{
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < bb.GetSize(); i++) {
    int index = bb[i];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = PDF(cur,index, dir);
    }
    for (int  dir = 0; dir < n_q_; dir++) {
      PDF(cur, index, GetInvDirection((Direction)dir)) = pdfs[dir];
    }
  }
  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::Prop_coll_densoutlet3D(int cur)
{
  double tmp_ux, tmp_uy, tmp_uz, tmp_us, sum, tmp;

  Vector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < outlet.GetSize(); i++) {
    int index = outlet[i];

    for (int dir = 0; dir < n_q_; dir ++)
      pdfs[dir] = PDF(cur,index,dir);

    CalcVelocities(pdfs, tmp_ux, tmp_uy, tmp_uz);

    sum = 1.0; // the enforced density
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    tmp_uz = 3.0 * tmp_uz;

    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
      PDF(cur, index, dir) = sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}

} // end namespace
