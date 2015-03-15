#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Driver/basedriver.hh"
#include "Driver/staticdriver.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/formsContext.hh"
#include "Optimization/Optimization.hh"
#include "PDE/NonFEM/LatticeBoltzmann.hh"
#include "Utils/Timer.hh"
#include "PDE/basePDE.hh"

namespace CoupledField
{

using std::fstream;
using std::ios;

DECLARE_LOG(lattice)
DEFINE_LOG(lattice, "lattice")

// instantiation of the static elements
Enum<LatticeBoltzmann::Direction>         LatticeBoltzmann::directions;

LatticeBoltzmann::LatticeBoltzmann(int dim, int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency)
{
  assert(dim == 2 || dim == 3);
  // n_q_: number of discrete directions in this model, e.g. 19 for D3Q19
  if (dim == 2) {
    assert(sizeZ == 1);
    n_q_ = 9;
  }
  else
    n_q_ = 19;

  m_dim = dim;
  m_sizeX = sizeX;
  m_sizeY = sizeY;
  m_sizeZ = sizeZ;
  m_ux = ux;
  m_uy = uy;
  m_uz = uz;
  m_omega = omega;
  m_maxIter = maxIterations;
  m_maxTol = maxTolerance;
  m_writeFrequency = writeFrequency;
  m_numWriteResults = 0;

  m_nNodes = m_sizeX * m_sizeY * m_sizeZ;

  m_plot = plot;

  //matrix of the probability distributions
  LOG_DBG(lattice) << "Allocating arrays for " << m_nNodes << " PDFs (" << (sizeof(double) * m_nNodes * n_q_ * 2.0 / 1024.0 / 1024.0) << " MiB)";

  m_pdfs.Resize(2);
  m_pdfs[0].Resize(m_nNodes * n_q_);
  m_pdfs[1].Resize(m_nNodes * n_q_);

  // microDirections stores information about the 19 microscopic velocities/directions of D3Q19 model
  microDirections.Resize(n_q_);

  weights.Resize(n_q_);

  SetEnums();
  InitWeights();
  SetInvDirections();
  TestInvDirections();

  Scales.Resize(m_nNodes);
  rel.Resize(m_nNodes);
  bb.Resize(2 * m_sizeX + 2 * m_sizeY + 2 * m_sizeZ);

  m_cur  = 0;
  m_next = 1;


  SetMicroVelocities();

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

  for (int dir = 0; dir < n_q_; dir++) {
    //store current pdf values in array for better accessing
    pdfs[dir] = pdf(cur, i, j, k, dir);
    tmp_ux += microDirections[dir].off_x*pdfs[dir];
    tmp_uy += microDirections[dir].off_y*pdfs[dir];
    tmp_uz += microDirections[dir].off_z*pdfs[dir];
  }

  ux = tmp_ux / density;
  uy = tmp_uy / density;
  uz = tmp_uz / density;
}

StdVector<double>* LatticeBoltzmann::Iterate(const StdVector<double>& elements, PtrParamNode in)
{
  // this flag cannot be set in constructor, since information about optimization is not available when this constructor is called
  if (domain->GetOptimization() == NULL || domain->GetOptimization()->GetOptimizerType() == Optimization::EVALUATE_INITIAL_DESIGN ) {
    writeIntermediateResults = true;
    std::cout << "Domain has no optimization!" << std::endl;
  }
  else
    writeIntermediateResults = false;

  int it = 0;

  // counts number of written steps when not in optimization
  int count = 1;

  double res = -1.;
  double R = 1.0;
  bool steady_state = false;

  std::ofstream plot;
  if(m_plot)
    plot.open(std::string(progOpts->GetSimName() + ".lbm.dat").c_str());

  InitializePdfs();
  SetupDataStructures(elements);

  assert((int) elements.GetSize() == m_nNodes);
  for (int i = 0; i < m_nNodes; ++i) {
    Scales[i] = 1.0 - elements[i];
  }

  Timer timer;
  timer.Start();

    LOG_DBG(lattice) << "bb = " << ToString(bb);
    LOG_DBG(lattice) << "inlet = " << ToString(inlet);
    LOG_DBG(lattice) << "outlet = " << ToString(outlet);
    LOG_DBG(lattice) << "rel = " << ToString(rel);

  in->Get("converged")->SetValue("running");

  while(it < m_maxIter && !steady_state && R <= 1000)
  {
    // -- Combined propagation and collision step -------------------------
    (this->*prop_coll_step)(m_cur, m_next, m_omega);

    // -- Bounce back step ------------------------------------------------
    (this->*prop_coll_bounce_back)(m_next, bb);
    // -- Inlet condition -------------------------------------------------
    (this->*prop_coll_velinlet)(m_next, inlet, m_ux, m_uy, m_uz);
    // -- Outlet condition ------------------------------------------------
    (this->*prop_coll_densoutlet)(m_next, outlet);
    if((it == 0 || it % 100 == 0))
    {
      //Calculation of the residual
      R = 0.;

      for (int elem = 0; elem < m_nNodes; elem++) {
        //            index = k * m_sizeX * m_sizeY + j * m_sizeX + i;
        for(int dir = 0; dir < n_q_; dir++)
        {
          res = pdf(m_next, elem, dir) - pdf(m_cur, elem, dir);
          R += res * res;
        }
      }

      R = sqrt(R);

      if(R <= m_maxTol)
        steady_state = true;

      if(m_plot)
        plot << it << "\t" << R << "\n";

      in->Get("iterations")->SetValue(it);
      in->Get("residuum")->SetValue(R);
      info->ToFile(); // is not written when called too often

//      std::cout << "residual: " << R << std::endl;
    }

    m_cur  = (m_cur  + 1) % 2;
    m_next = (m_next + 1) % 2;

    it++;

    if (writeIntermediateResults) {
//      std::cout << "writing intermediate LBM results..." << std::endl;
      // check if intermediate results should write to hdf5 file
      if (it % m_writeFrequency == 0) {
        // in case of optimization we need to use StoreResults() from optimization to avoid conflicts with ...
  //      if(domain->GetOptimization() != NULL)
  //        domain->GetOptimization()->DirectStoreResults(domain->GetOptimization()->GetCurrentIteration() + 1/100000.0 * count); // <iter>.00001, ... <iter>.000301, ....
  //      else
        domain->GetDriver()->StoreResults(count,(double) count);
        count++;
      }
    }
  }

  // -- Combined propagation and collision step -------------------------
  (this->*prop_coll_step)(m_cur, m_next, m_omega);

  // -- Bounce back step ------------------------------------------------
  (this->*prop_coll_bounce_back)(m_next, bb);
  // -- Inlet condition -------------------------------------------------
  (this->*prop_coll_velinlet)(m_next, inlet, m_ux, m_uy, m_uz);
  // -- Outlet condition ------------------------------------------------
  (this->*prop_coll_densoutlet)(m_next, outlet);

  timer.Stop();

  in->Get("iterations")->SetValue(it);
  in->Get("residuum")->SetValue(R);
  in->Get("converged")->SetValue(steady_state);

  if(m_plot) {
    plot << it << "\t" << R << "\n";
    plot.flush();
  }

  double wt = timer.GetWallTime();
  double performance;
  if (m_dim == 3)
    performance = (m_sizeX - 1) * (m_sizeY - 1) * (m_sizeZ - 1) * it / wt / 1e6;
  else
    performance = (m_sizeX - 1) * (m_sizeY - 1) * it / wt / 1e6;

  in->Get("MFLUP_s")->SetValue(performance);
  in->Get("walltime")->SetValue(wt);

  if(R >= 1000)
    EXCEPTION("In LBM iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("internal LBM simulation could not converge: iterations: " << it << " residuum: " << R);

  m_numWriteResults = count;

  return &(m_pdfs[m_cur]);
}

void LatticeBoltzmann::InitializePdfs()
{
#pragma omp parallel for default(none)
  for (int elem = 0; elem < m_nNodes; elem++) {
    for (int dir = 0; dir < n_q_; dir++) {
      pdf(0, elem, dir) = weights[dir];
      pdf(1, elem, dir) = weights[dir];
    }
  }
}

//
// Checks if velocities on the boundary are really zero (no-slip b.c.)
//

void LatticeBoltzmann::CheckBoundaryVelocities(int cur, StdVector<StdVector<int> >& bb)
{
  int x, y, z;
  double ux, uy, uz;
  double eps = 1e-10;
  for(unsigned int i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];
    z = bb[i][2];
    CalcVelocitites(cur,x, y, z, ux, uy, uz);
    if (ux > eps)
      std::cout << "ux > eps. ux: " << ux << std::endl;
    if (uy > eps)
      std::cout << "uy > eps. ux: " << uy << std::endl;
    if (uz > eps)
      std::cout << "uz > eps. ux: " << uz << std::endl;
    assert(ux < eps);
    assert(uy < eps);
    assert(uz < eps);
  }
}

void LatticeBoltzmann::SetMicroVelocities()
{
  microDirections[Q_0] = MicroVelocity(0,0,0);
  microDirections[Q_E] = MicroVelocity(1,0,0);
  microDirections[Q_W] = MicroVelocity(-1,0,0);
  microDirections[Q_N] = MicroVelocity(0,1,0);
  microDirections[Q_S] = MicroVelocity(0,-1,0);
  microDirections[Q_NE] = MicroVelocity(1,1,0);
  microDirections[Q_NW] = MicroVelocity(-1,1,0);
  microDirections[Q_SE] = MicroVelocity(1,-1,0);
  microDirections[Q_SW] = MicroVelocity(-1,-1,0);
  if (m_dim == 3) {
    microDirections[Q_T] = MicroVelocity(0,0,1);
    microDirections[Q_B] = MicroVelocity(0,0,-1);
    microDirections[Q_TN] = MicroVelocity(0,1,1);
    microDirections[Q_TS] = MicroVelocity(0,-1,1);
    microDirections[Q_TE] = MicroVelocity(1,0,1);
    microDirections[Q_TW] = MicroVelocity(-1,0,1);
    microDirections[Q_BN] = MicroVelocity(0,1,-1);
    microDirections[Q_BS] = MicroVelocity(0,-1,-1);
    microDirections[Q_BE] = MicroVelocity(1,0,-1);
    microDirections[Q_BW] = MicroVelocity(-1,0,-1);
  }
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

  if (m_dim == 3) {
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
  inverseDirections.Resize(n_q_);

  inverseDirections[Q_0] = Q_0;
  inverseDirections[Q_E] = Q_W;
  inverseDirections[Q_N] = Q_S;
  inverseDirections[Q_W] = Q_E;
  inverseDirections[Q_S] = Q_N;
  inverseDirections[Q_NE] = Q_SW;
  inverseDirections[Q_SW] = Q_NE;
  inverseDirections[Q_NW] = Q_SE;
  inverseDirections[Q_SE] = Q_NW;

  if (m_dim == 3) {
    inverseDirections[Q_T] = Q_B;
    inverseDirections[Q_B] = Q_T;
    inverseDirections[Q_TN] = Q_BS;
    inverseDirections[Q_BS] = Q_TN;
    inverseDirections[Q_TS] = Q_BN;
    inverseDirections[Q_BN] = Q_TS;
    inverseDirections[Q_TE] = Q_BW;
    inverseDirections[Q_BW] = Q_TE;
    inverseDirections[Q_TW] = Q_BE;
    inverseDirections[Q_BE] = Q_TW;
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

  if (m_dim == 3) {
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

  for(int k = 0; k < m_sizeZ; k++)
  {
    for(int j = 0; j < m_sizeY; j++)
    {
      for(int i = 0; i < m_sizeX; i++)
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
        }
        ++n;
      }
    }
  }
  assert(m_nNodes == n);
}


std::string LatticeBoltzmann::ToString(const StdVector<StdVector<int> >& data)
{
  std::stringstream ss;

  for(unsigned int e = 0; e < data.GetSize(); e++)
  {
    ss << e << ": (";
    for(unsigned int d = 0; d < data[e].GetSize(); d++)
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

  for(int i = 0; i < m_nNodes; i++) {
    for(int j = 0; j < n_q_; j++) {
      f << pdf(cur, i, j) << " ";
    }

    f << std::endl;
  }

  f.close();
}


void LatticeBoltzmann::prop_step()
{

  // perform a propagation step
  int x, y, z;

  int lx = m_sizeX;
  int ly = m_sizeY;
  int lz = m_sizeZ;

  for(z = lz - 1; z >= 0; z--) {
    for(y = ly - 1; y >= 0; y--) {
      for(x = lx - 1; x > 0; x--) {
        pdf(m_cur, x, y, z, Q_E)  = pdf(m_cur, x - 1  , y     , z     , Q_E);
      }
      for(x = 0; x < lx - 1; x++) {
        pdf(m_cur, x, y, z, Q_W)  = pdf(m_cur, x + 1  , y     , z     , Q_W);
      }
    }
    for(y = ly - 1; y > 0; y--) {
      for(x = lx - 1; x > 0; x--) {
        pdf(m_cur, x, y, z, Q_NE) = pdf(m_cur, x - 1  , y - 1 , z     , Q_NE);
      }
      for(x = 0; x < lx - 1; x++) {
        pdf(m_cur, x, y, z, Q_NW) = pdf(m_cur, x + 1  , y - 1 , z     , Q_NW);
      }
    }
    for(y = 0; y < ly - 1; y++) {
      for(x = lx - 1; x > 0; x--) {
        pdf(m_cur, x, y, z, Q_SE) = pdf(m_cur, x -  1 , y + 1 , z     , Q_SE);
      }
      for(x = 0; x < lx - 1; x++) {
        pdf(m_cur, x, y, z, Q_SW) = pdf(m_cur, x + 1  , y + 1 , z     , Q_SW);
      }
    }
    for(x = 0; x < lx; x++) {
      for(y = ly - 1; y > 0; y--) {
        pdf(m_cur, x, y, z, Q_N)  = pdf(m_cur, x      , y - 1 , z     , Q_N);
      }
      for(y = 0; y < ly - 1; y++) {
        pdf(m_cur, x, y, z, Q_S)  = pdf(m_cur, x      , y + 1 , z     , Q_S);
      }
    }
  }

  if (m_dim == 2) return;

  for(z = lz - 1; z > 0; z--) {
    for(y = ly - 1; y >= 0; y--) {
      for(x = lx - 1; x >= 0; x--) {
        pdf(m_cur, x, y, z, Q_T)  = pdf(m_cur, x      , y     , z - 1 , Q_T);
      }
    }
    for(x = lx - 1; x >= 0; x--) {
      for(y = ly - 1; y > 0; y--) {
        pdf(m_cur, x, y, z, Q_TN) = pdf(m_cur, x      , y - 1 , z - 1 , Q_TN);
      }
      for(y = 0; y < ly - 1; y++) {
        pdf(m_cur, x, y, z, Q_TS) = pdf(m_cur, x      , y + 1 , z - 1 , Q_TS);
      }
    }
    for(y = ly - 1; y >= 0; y--) {
      for(x = lx - 1; x > 0; x--) {
        pdf(m_cur, x, y, z, Q_TE) = pdf(m_cur, x - 1  , y     , z - 1 , Q_TE);
      }
      for(x = 0; x < lx - 1; x++) {
        pdf(m_cur, x, y, z, Q_TW) = pdf(m_cur, x + 1  , y     , z - 1 , Q_TW);
      }
    }
  }

  for(z = 0; z < lz - 1; z++) {
    for(y = ly - 1; y >= 0; y--) {
      for(x = lx - 1; x >= 0; x--) {
        pdf(m_cur, x, y, z, Q_B)  = pdf(m_cur, x      , y     , z + 1 , Q_B);
      }
    }
    for(x = lx - 1; x >= 0; x--) {
      for(y = 0; y < ly - 1; y++) {
        pdf(m_cur, x, y, z, Q_BS) = pdf(m_cur, x      , y + 1 , z + 1 , Q_BS);
      }
      for(y = ly - 1; y > 0; y--) {
        pdf(m_cur, x, y, z, Q_BN) = pdf(m_cur, x      , y - 1 , z + 1 , Q_BN);
      }
    }
    for(y = ly - 1; y >= 0; y--) {
      for(x = lx - 1; x > 0; x--) {
        pdf(m_cur, x, y, z, Q_BE) = pdf(m_cur, x - 1  , y     , z + 1 , Q_BE);
      }
      for(x = 0; x < lx - 1; x++) {
        pdf(m_cur, x, y, z, Q_BW) = pdf(m_cur, x + 1  , y     , z + 1 , Q_BW);
      }
    }
  }
  return;
}

// Calculates macroscopic density for given element
double LatticeBoltzmann::CalcDensity(int cur, int i, int j, int k)
{
  double sum = 0;
  for (int dir = 0; dir < n_q_; dir++) {
    sum += pdf(cur, i, j, k, dir);
  }

  // debugging
  if (sum > 1.5 || sum < 1e-8) {
    std::cout << "i: " << i << " j: " << j << " k: " << j << " sum: " << sum << std::endl;
    for (int dir = 0; dir < n_q_; dir++) {
      std::cout << "dir " << dir << " pdf " <<  pdf(cur, i, j, k, dir) << std::endl;
    }
  }

  assert(sum < 1.5);
  assert(sum > 1e-8);
  return sum;
}

/************************************************** 2D operators *****************************************************/

void LatticeBoltzmann::prop_coll_step2D(int cur, int next, double omega)
{
  int x, y, z = 0;

  double tmp, tmp_ux, tmp_uy, tmp_us, scale, sum;
  double * scales  = Scales.GetPointer();

  int index;

  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  int tmp_x,tmp_y;

#pragma omp parallel default(none)\
 private(index), \
 private(tmp_ux, tmp_uy, tmp_us, scale, sum, tmp, x, y, tmp_x, tmp_y), \
 shared(next, cur, scales, omega, z)
{
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  #pragma omp for collapse(2)
  for (y = 0; y < m_sizeY ; ++y) {
    for (x = 0; x < m_sizeX ; ++x) {

      index= GetIndex(x,y,z);

      // sum: macroscopic density is sum over all discrete distributions of an element
      sum = 0;
      tmp_ux = 0;
      tmp_uy = 0;
      tmp_us = 0;

      for (int dir = 0; dir < n_q_; dir++) {

        if (OutsideDomain(x,y,inverseDirections[dir])) {
          tmp_x = x;
          tmp_y = y;
        }
        // else: standard propagation
        else {
          tmp_x = microDirections[inverseDirections[dir]].off_x + x;
          tmp_y = microDirections[inverseDirections[dir]].off_y + y;
        }

        //store current pdf values in array for better accessing
        pdfs[dir] = pdf(cur, tmp_x, tmp_y,  z, dir);
        sum += pdfs[dir];
        tmp_ux += microDirections[dir].off_x*pdfs[dir];
        tmp_uy += microDirections[dir].off_y*pdfs[dir];
      }

      // macroscopic scaling by design variable
      scale = scales[index];

      tmp_ux = scale * tmp_ux / sum;
      tmp_uy = scale * tmp_uy / sum;
      tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);

      tmp_ux = 3.0 * tmp_ux;
      tmp_uy = 3.0 * tmp_uy;

      // propagation and collision in one step
      for (int dir = 0; dir < n_q_; dir++) {
        tmp = microDirections[dir].off_x * tmp_ux + microDirections[dir].off_y * tmp_uy ;
        // no collision on the boundaries
        if (x == 0 || y == 0 || x == m_sizeX - 1 || y == m_sizeY - 1)
          pdf(next, x, y, z, dir) = pdfs[dir];
        else
          pdf(next, x, y, z, dir) = pdfs[dir] + omega * ((sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdfs[dir]);
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

  for(unsigned int i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];

    assert(z == 0);

    sum = CalcDensity(cur, x, y, z);

    tmp_ux = UX;
    tmp_uy = UY;

    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;

    LOG_DBG3(lattice) << "pcv: i=" << i << " tux=" << tmp_ux << " tuy=" << tmp_uy << std::endl;

    for (int dir = 0; dir < n_q_; dir++) {
      tmp = microDirections[dir].off_x * tmp_ux + microDirections[dir].off_y * tmp_uy;
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

  for(unsigned int i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];

    for (int dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = pdf(cur, x, y, z, dir);
    }
    for (int dir = 0; dir < n_q_; dir++) {
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


  for(unsigned int i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];

    assert (z == 0);

    CalcVelocitites(cur, x, y, z, tmp_ux, tmp_uy, tmp);

    sum = 1.0; // the enforced density
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;

    for (int dir = 0; dir < n_q_; dir++) {
      tmp = microDirections[dir].off_x * tmp_ux + microDirections[dir].off_y * tmp_uy;
      pdf(cur, x, y, z, dir) = sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


/************************************************** 3D operators *****************************************************/
void LatticeBoltzmann::prop_coll_step3D(int cur, int next, double omega)
{

  int x, y, z;
  double tmp, tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum;
  double * scales  = Scales.GetPointer();

  int index;

  int tmp_x, tmp_y, tmp_z;

#pragma omp parallel default(none)\
    private(index), \
    private(tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum, tmp, x, y,z, tmp_x, tmp_y, tmp_z), \
    shared(next, cur, scales, omega)
{
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  #pragma omp for collapse(3)
  for (z = 0; z < m_sizeZ; ++z) {
    for (y = 0; y < m_sizeY; ++y) {
      for (x = 0; x < m_sizeX; ++x) {
        index= GetIndex(x,y,z);

        // sum: macroscopic density is sum over all discrete distributions of an element
        sum = 0;
        tmp_ux = 0;
        tmp_uy = 0;
        tmp_uz = 0;
        tmp_us = 0;

        for (int dir = 0; dir < n_q_; dir++) {

          if (OutsideDomain(x,y,z,inverseDirections[dir])) { // boundary case
            tmp_x = x;
            tmp_y = y;
            tmp_z = z;
          }
          else { // standard propagation rule
            tmp_x = microDirections[inverseDirections[dir]].off_x + x;
            tmp_y = microDirections[inverseDirections[dir]].off_y + y;
            tmp_z = microDirections[inverseDirections[dir]].off_z + z;
          }

          //store current pdf values in array for better accessing
          pdfs[dir] = pdf(cur, tmp_x, tmp_y, tmp_z, dir);
          sum += pdfs[dir];
          tmp_ux += microDirections[dir].off_x*pdfs[dir];
          tmp_uy += microDirections[dir].off_y*pdfs[dir];
          tmp_uz += microDirections[dir].off_z*pdfs[dir];
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

        for (int dir = 0; dir < n_q_; dir++) {
          tmp = microDirections[dir].off_x * tmp_ux + microDirections[dir].off_y * tmp_uy + microDirections[dir].off_z * tmp_uz;
          if (x == 0 || y == 0 || x == m_sizeX - 1 || y == m_sizeY - 1 || z == 0 || z == m_sizeZ - 1)
            pdf(next, x, y, z, dir) = pdfs[dir];
          else
            pdf(next, x, y, z, dir) = pdfs[dir] + omega * ((sum * weights[dir]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdfs[dir]);
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

  for(unsigned int i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];
    z = inlet[i][2];

    sum = CalcDensity(cur, x, y, z);

    tmp_ux = UX;
    tmp_uy = UY;
    tmp_uz = UZ;
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    tmp_uz = 3.0 * tmp_uz;

    LOG_DBG3(lattice) << "pcv: i=" << i << " tux=" << tmp_ux << " tuy=" << tmp_uy << " tuz=" << tmp_uz << std::endl;


    for (int dir = 0; dir < n_q_; dir++) {
      tmp = microDirections[dir].off_x * tmp_ux + microDirections[dir].off_y * tmp_uy + microDirections[dir].off_z * tmp_uz;
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

  for(unsigned int i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];
    z = bb[i][2];

    for (int dir = 0; dir < n_q_; dir++) {
      pdfs[dir] = pdf(cur, x, y, z, dir);
    }
    for (int dir = 0; dir < n_q_; dir++) {
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

  for(unsigned int i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];
    z = outlet[i][2];

    CalcVelocitites(cur, x, y, z, tmp_ux, tmp_uy, tmp_uz);

    sum = 1.0; // the enforced density
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy + tmp_uz * tmp_uz);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;
    tmp_uz = 3.0 * tmp_uz;

    for (int dir = 0; dir < n_q_; dir++) {
      tmp = microDirections[dir].off_x * tmp_ux + microDirections[dir].off_y * tmp_uy + microDirections[dir].off_z * tmp_uz;
      pdf(cur, x, y, z, dir) = sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}

} // end namespace
