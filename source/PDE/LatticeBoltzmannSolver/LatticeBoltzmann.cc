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

DECLARE_LOG(lattice)
DEFINE_LOG(lattice, "lattice")

// instantiation of the static elements
Enum<LatticeBoltzmann::Direction>         LatticeBoltzmann::directions;

LatticeBoltzmann::LatticeBoltzmann(int dim, int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, StdVector< StdVector<double> > uin, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency)
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
  omega_ = omega;
  maxIter_ = maxIterations;
  maxTol_ = maxTolerance;
  writeFrequency_ = writeFrequency;
  numWriteResults_ = 0;

  nNodes_ = sizeX_ * sizeY_ * sizeZ_;

  u_in.Resize(nNodes_);
  u_in = uin;

  plot_ = plot;

  //matrix of the probability distributions
  LOG_DBG(lattice) << "Allocating arrays for " << nNodes_ << " PDFs (" << (sizeof(double) * nNodes_ * n_q_ * 2.0 / 1024.0 / 1024.0) << " MiB)";

  pdfs.Resize(2);
  pdfs[0].Resize(nNodes_ * n_q_);
  pdfs[1].Resize(nNodes_ * n_q_);

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

//  InitTransformMatrix();

  //initlialize function pointers in dependence on problem's dimension
  if (dim == 2) {
    prop_coll_step = &LatticeBoltzmann::prop_coll_step2D;
    prop_coll_velinlet = &LatticeBoltzmann::prop_coll_velinlet2D;
    prop_coll_bounce_back = &LatticeBoltzmann::prop_coll_bounce_back2D;
    prop_coll_densoutlet = &LatticeBoltzmann::prop_coll_densoutlet2D;
  }
  else
  {
    prop_coll_step = &LatticeBoltzmann::prop_coll_step3D;
    prop_coll_velinlet = &LatticeBoltzmann::prop_coll_velinlet3D;
    prop_coll_bounce_back = &LatticeBoltzmann::prop_coll_bounce_back3D;
    prop_coll_densoutlet = &LatticeBoltzmann::prop_coll_densoutlet3D;
  }
}

LatticeBoltzmann::~LatticeBoltzmann()
{
}

void LatticeBoltzmann::CalcVelocitites(int cur, int i, int j, int k, double& ux, double& uy, double& uz)
{

  double density = CalcDensity(cur, i, j, k);
  double tmp_ux = 0;
  double tmp_uy = 0;
  double tmp_uz = 0;
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for (int  dir = 0; dir < n_q_; dir++) {
    //store current pdf values in array for better accessing
    pdfs[dir] = pdf(cur, i, j, k, dir);
    tmp_ux += microVelDirections[dir].off_x*pdfs[dir];
    tmp_uy += microVelDirections[dir].off_y*pdfs[dir];
    tmp_uz += microVelDirections[dir].off_z*pdfs[dir];
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
  }


  Timer timer;
  timer.Start();

  LOG_DBG(lattice) << "bb = " << ToString(bb);
  LOG_DBG(lattice) << "inlet = " << ToString(inlet);
  LOG_DBG(lattice) << "outlet = " << ToString(outlet);
  LOG_DBG(lattice) << "rel = " << ToString(rel);

  in->Get("converged")->SetValue("running");

  while(it < maxIter_ && !steady_state && R <= 1000)
  {
    // -- Combined propagation and collision step -------------------------
    (this->*prop_coll_step)(cur_, next_);

    // -- Bounce back step ------------------------------------------------
    (this->*prop_coll_bounce_back)(next_, bb);
    // -- Inlet condition -------------------------------------------------
    (this->*prop_coll_velinlet)(next_, inlet, ux_, uy_, uz_);
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
          res = pdf(next_, elem, dir) - pdf(cur_, elem, dir);
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

  return &(pdfs[cur_]);
}

void LatticeBoltzmann::InitializePdfs()
{
#pragma omp parallel for default(none)
  for (int elem = 0; elem < nNodes_; elem++) {
    for (int  dir = 0; dir < n_q_; dir++) {
      pdf(0, elem, dir) = weights[dir];
      pdf(1, elem, dir) = weights[dir];
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
  transformation.Resize(n_q_ * n_q_);
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

  relax_rates.Resize(n_q_ * n_q_);
  relax_rates.Init();
  // relaxation rates matrix S is diagonal
  // entries for density and momentum components are 0 due to collision invariance
  // the first 4 values here, can be chosen arbritrarily between 0 and 2
  relax_rates[1][1] = 1.9;
  relax_rates[2][2] = 1.9;
  relax_rates[4][4] = 1.9;
  relax_rates[6][6] = 1.9;
  // these two values are related to fluid's kinematic viscosity
  relax_rates[7][7] = omega_;
  relax_rates[8][8] = omega_;
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
    invPDFDirections[Q_TN] = Q_BS;
    invPDFDirections[Q_BS] = Q_TN;
    invPDFDirections[Q_TS] = Q_BN;
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
        } else if (LbmNodeTypeIsBB(porosity)) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          bb.Push_back(tmp);
        } else if (LbmNodeTypeIsInlet(porosity)) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          inlet.Push_back(tmp);
        } else if (LbmNodeTypeIsOutlet(porosity)) {
          tmp[0] = i;
          tmp[1] = j;
          tmp[2] = k;
          outlet.Push_back(tmp);
        } else if (LbmNodeIsObstacle(porosity)) {
          obst.Push_back(GetIndex(i,j,k));
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

void LatticeBoltzmann::create_output(const char * file, int cur)
{
  // for debug purposes
  fstream f;
  f.precision(16);
  f.open(file, ios::out);

  for(int i = 0; i < nNodes_; i++) {
    for(int  j = 0; j < n_q_; j++) {
      f << pdf(cur, i, j) << " ";
    }

    f << std::endl;
  }

  f.close();
}


void LatticeBoltzmann::prop_step()
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

          pdf(next_,x,y,z,dir) = pdf(cur_, tmp_x, tmp_y,  tmp_z, dir);
        }
      }
    }
  }

  pdfs[cur_] = pdfs[next_];
}

// Calculates macroscopic density for given element
double LatticeBoltzmann::CalcDensity(int cur, int i, int j, int k)
{
  double sum = 0;
  for (int  dir = 0; dir < n_q_; dir++) {
    sum += pdf(cur, i, j, k, dir);
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

void LatticeBoltzmann::transform_pdfs(int cur)
{
  Vector<double> tmp;
  Vector<double> result;
  tmp.Resize(n_q_);
  result.Resize(n_q_);

  for (int elem = 0; elem < nNodes_; elem++)
  {
    for (int  dir = 0; dir < n_q_; dir++)
      tmp[dir] = pdf(cur,elem,dir);
    transformation.Mult(tmp,result);
    for (int  dir = 0; dir < n_q_; dir++)
      pdfs_moments[GetPdfIndex(elem,dir)] = result[dir];
  }
}

/************************************************** 2D operators *****************************************************/

void LatticeBoltzmann::prop_coll_step2D(int cur, int next)
{
  int x, y, z = 0;

  double tmp, tmp_ux, tmp_uy, tmp_us, scale, sum;

  int index;

  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  int tmp_x,tmp_y;

#pragma omp parallel default(none)\
 private(index), \
 private(tmp_ux, tmp_uy, tmp_us, scale, sum, tmp, x, y, tmp_x, tmp_y), \
 shared(next, cur, z)
{
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  #pragma omp for collapse(2)
  for (y = 0; y < sizeY_ ; ++y) {
    for (x = 0; x < sizeX_ ; ++x) {

      index= GetIndex(x,y,z);

      // sum: macroscopic density is sum over all discrete distributions of an element
      sum = 0;
      tmp_ux = 0;
      tmp_uy = 0;
      tmp_us = 0;

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
        pdfs[dir] = pdf(cur, tmp_x, tmp_y,  z, dir); // accessed pdf value can be the old one or the neighbor's value
        sum += pdfs[dir];
        tmp_ux += microVelDirections[dir].off_x*pdfs[dir];
        tmp_uy += microVelDirections[dir].off_y*pdfs[dir];
      }

      // macroscopic scaling by design variable
      scale = scales[index];

      tmp_ux = scale * tmp_ux / sum;
      tmp_uy = scale * tmp_uy / sum;
      tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);

      tmp_ux = 3.0 * tmp_ux;
      tmp_uy = 3.0 * tmp_uy;

      // propagation and collision in one step
      for (int  dir = 0; dir < n_q_; dir++) {
        tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy ;
        // no collision on the boundaries
        if (x == 0 || y == 0 || x == sizeX_ - 1 || y == sizeY_ - 1)
          pdf(next, x, y, z, dir) = pdfs[dir];
        else
          pdf(next, x, y, z, dir) = pdfs[dir] + omega_ * (sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us) - pdfs[dir]);
      }
    }
  }
}
  return;
}


void LatticeBoltzmann::prop_coll_velinlet2D(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY, double UZ)
{
  assert(UZ == 0);

  int x, y, z = 0;
  double tmp_ux, tmp_uy, tmp_us, sum;
  double tmp;

  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];

    assert(z == 0);

    sum = CalcDensity(cur, x, y, z);

    tmp_ux = UX; // velocity at inlet is prescribed
    tmp_uy = UY;
    if (u_in.GetSize() != 0) {
      tmp_ux = u_in[i][0];
      tmp_uy = u_in[i][1];
    }

    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;

    LOG_DBG3(lattice) << "pcv: i=" << i << " tux=" << tmp_ux << " tuy=" << tmp_uy << std::endl;

    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy;
      pdf(cur, x, y, z, dir)  = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }
  return;
}

//
// Performs a bounce back step.
//
void LatticeBoltzmann::prop_coll_bounce_back2D(int cur, StdVector<StdVector<int> >& bb)
{
  int x, y, z = 0;
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = pdf(cur, x, y, z, dir);
    }
    for (int  dir = 0; dir < n_q_; dir++) {
      pdf(cur, x, y, z, GetInvDirection((Direction)dir)) = pdfs[dir];
    }
  }

  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::prop_coll_densoutlet2D(int cur, StdVector<StdVector<int> >& outlet)
{
  double tmp_ux, tmp_uy, tmp_us, sum, tmp;

  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  int x, y, z = 0;


  for(unsigned int  i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];

    assert (z == 0);

    CalcVelocitites(cur, x, y, z, tmp_ux, tmp_uy, tmp);

    sum = 1.0; // the enforced density
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;

    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy;
      pdf(cur, x, y, z, dir) =  sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


/************************************************** 3D operators *****************************************************/
void LatticeBoltzmann::prop_coll_step3D(int cur, int next)
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
          pdfs[dir] = pdf(cur, tmp_x, tmp_y, tmp_z, dir);
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
            pdf(next, x, y, z, dir) = pdfs[dir];
          else
            pdf(next, x, y, z, dir) = pdfs[dir] + omega_ * ((sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdfs[dir]);
        }

      }
    }
  }
}
  return;
}

void LatticeBoltzmann::prop_coll_velinlet3D(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY, double UZ)
{
  int x, y, z;
  double tmp_ux, tmp_uy, tmp_uz, tmp_us, sum;
  double tmp;

  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];
    z = inlet[i][2];

    sum = CalcDensity(cur, x, y, z);

    tmp_ux = UX;
    tmp_uy = UY;
    tmp_uz = UZ;

    if (!u_in.IsEmpty()) { //use parabolic profile
      tmp_ux = u_in[i][0];
      tmp_uy = u_in[i][1];
      tmp_uz = u_in[i][2];
    }
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    tmp_uz = 3.0 * tmp_uz;

    LOG_DBG3(lattice) << "pcv: i=" << i << " tux=" << tmp_ux << " tuy=" << tmp_uy << " tuz=" << tmp_uz << std::endl;


    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
      pdf(cur, x, y, z, dir)  = sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


//
// Performs a bounce back step.
//
void LatticeBoltzmann::prop_coll_bounce_back3D(int cur, StdVector<StdVector<int> >& bb)
{

  int x, y, z;
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  for(unsigned int  i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];
    z = bb[i][2];

    for (int  dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = pdf(cur, x, y, z, dir);
    }
    for (int  dir = 0; dir < n_q_; dir++) {
      pdf(cur, x, y, z, GetInvDirection((Direction)dir)) = pdfs[dir];
    }
  }
  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::prop_coll_densoutlet3D(int cur, StdVector<StdVector<int> >& outlet)
{
  double tmp_ux, tmp_uy, tmp_uz, tmp_us, sum, tmp;

  StdVector<double> pdfs;
  pdfs.Resize(n_q_);

  int x, y, z;

  for(unsigned int  i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];
    z = outlet[i][2];

    CalcVelocitites(cur, x, y, z, tmp_ux, tmp_uy, tmp_uz);

    sum = 1.0; // the enforced density
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    tmp_uz = 3.0 * tmp_uz;

    for (int  dir = 0; dir < n_q_; dir++) {
      tmp = microVelDirections[dir].off_x * tmp_ux + microVelDirections[dir].off_y * tmp_uy + microVelDirections[dir].off_z * tmp_uz;
      pdf(cur, x, y, z, dir) = sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}

} // end namespace
