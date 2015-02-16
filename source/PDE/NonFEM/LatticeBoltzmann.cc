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
Enum<LatticeBoltzmann::Boundary>        LatticeBoltzmann::boundaries;
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

  // 8 corners + 12 edges + +6 faces + 1 map for interior nodes
  prop_maps.Resize(n_q_ + 8);
  for (int dir = 0; dir < n_q_ + 8; dir++ ){
    prop_maps[dir].Resize(n_q_ + 8);
  }

  // velocityDirections stores information about the 19 microscopic velocities/directions of D3Q19 model
  velocityDirections.Resize(n_q_);

  weights.Resize(n_q_);

  SetEnums();
  InitWeights();
  SetInvDirections();
  TestDirectionIndex();
  TestInvDirections();
  TestGetIndexbound();

  Scales.Resize(m_nNodes);
  rel.Resize(m_nNodes);
  bb.Resize(2 * m_sizeX + 2 * m_sizeY + 2 * m_sizeZ);

  m_cur  = 0;
  m_next = 1;


  //initlialize function pointers in dependence on problem's dimension
  if (dim == 2) {
    SetupTransformation2D();
    prop_coll_step = &LatticeBoltzmann::prop_coll_step2D;
    prop_coll_velinlet = &LatticeBoltzmann::prop_coll_velinlet2D;
    prop_coll_bounce_back = &LatticeBoltzmann::prop_coll_bounce_back2D;
    prop_coll_densoutlet = &LatticeBoltzmann::prop_coll_densoutlet2D;
  }
  else
  {
    SetupTransformation3D();
    prop_coll_step = &LatticeBoltzmann::prop_coll_step3D;
    prop_coll_velinlet = &LatticeBoltzmann::prop_coll_velinlet3D;
    prop_coll_bounce_back = &LatticeBoltzmann::prop_coll_bounce_back3D;
    prop_coll_densoutlet = &LatticeBoltzmann::prop_coll_densoutlet3D;
  }

  // for debugging
//  WritePropMap();
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
    tmp_ux += velocityDirections[dir].off_x*pdfs[dir];
    tmp_uy += velocityDirections[dir].off_y*pdfs[dir];
    tmp_uz += velocityDirections[dir].off_z*pdfs[dir];
  }

  ux = tmp_ux / density;
  uy = tmp_uy / density;
  uz = tmp_uz / density;
}

StdVector<double>* LatticeBoltzmann::Iterate(const StdVector<double>& elements, PtrParamNode in)
{
  int it = 0;

  // counts number of written steps
  int count = 0;

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

    // check if intermediate results should the written to hdf5 file
    if (it % m_writeFrequency == 0) {
      count++;
      domain->GetDriver()->StoreResults(count,(double)count);
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

// validates GetIndexbound function
void LatticeBoltzmann::TestGetIndexbound()
{
  assert(GetIndexBound(Q_B, Q_N, Q_W) == CORNER_BNW);
  assert(GetIndexBound(Q_B,Q_S,Q_W) == CORNER_BSW);
  assert(GetIndexBound(Q_B,Q_S,Q_E) == CORNER_BSE);
  assert(GetIndexBound(Q_B,Q_N,Q_E) == CORNER_BNE);
  assert(GetIndexBound(Q_T,Q_N,Q_W) == CORNER_TNW);
  assert(GetIndexBound(Q_T,Q_S,Q_W) == CORNER_TSW);
  assert(GetIndexBound(Q_T,Q_S,Q_E) == CORNER_TSE);
  assert(GetIndexBound(Q_T,Q_N,Q_E) == CORNER_TNE);
  assert(GetIndexBound(Q_T,Q_N) == EDGE_TN);
  assert(GetIndexBound(Q_T,Q_W) == EDGE_TW);
  assert(GetIndexBound(Q_T,Q_S) == EDGE_TS);
  assert(GetIndexBound(Q_T,Q_E) == EDGE_TE);
  assert(GetIndexBound(Q_N,Q_W) == EDGE_NW);
  assert(GetIndexBound(Q_S,Q_W) == EDGE_SW);
  assert(GetIndexBound(Q_S,Q_E) == EDGE_SE);
  assert(GetIndexBound(Q_N,Q_E) == EDGE_NE);
  assert(GetIndexBound(Q_B,Q_N) == EDGE_BN);
  assert(GetIndexBound(Q_B,Q_W) == EDGE_BW);
  assert(GetIndexBound(Q_B,Q_S) == EDGE_BS);
  assert(GetIndexBound(Q_B,Q_E) == EDGE_BE);
  assert(GetIndexBound(Q_W) == FACE_W);
  assert(GetIndexBound(Q_E) == FACE_E);
  assert(GetIndexBound(Q_T) == FACE_T);
  assert(GetIndexBound(Q_B) == FACE_B);
  assert(GetIndexBound(Q_S) == FACE_S);
  assert(GetIndexBound(Q_N) == FACE_N);
}


// validates GetIndexDir function
void LatticeBoltzmann::TestDirectionIndex()
{
  assert(GetIndexDir(Q_N, Q_E) == Q_NE);
  assert(GetIndexDir(Q_S, Q_W) == Q_SW);
  assert(GetIndexDir(Q_N, Q_W) == Q_NW);

  if (m_dim == 3) {
    assert(GetIndexDir(Q_S, Q_E) == Q_SE);
    assert(GetIndexDir(Q_T, Q_N) == Q_TN);
    assert(GetIndexDir(Q_B, Q_S) == Q_BS);
    assert(GetIndexDir(Q_T, Q_S) == Q_TS);
    assert(GetIndexDir(Q_B, Q_N) == Q_BN);
    assert(GetIndexDir(Q_T, Q_E) == Q_TE);
    assert(GetIndexDir(Q_B, Q_W) == Q_BW);
    assert(GetIndexDir(Q_T, Q_W) == Q_TW);
    assert(GetIndexDir(Q_B, Q_E) == Q_BE);
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

void LatticeBoltzmann::SetupTransformation2D()
{
  // in 2D: 8 maps (4 corners + 4 edges)
  StdVector<PropTransform> map;
  map.Resize(n_q_);
  const Direction xDirs[] = {Q_E,Q_W};
  const Direction yDirs[] = {Q_N,Q_S};
  Direction d1, d2;

  // propagation in corners
  for (int i = 0; i < 2; i++) { // E, W
    for (int j = 0; j < 2; j++) { // N, S
      // fill map to default offsets
      InitPropMap(map);
      d1 = xDirs[i];
      map[GetIndexDir(d1)] = PropTransform((d1 == Q_W ? 1: -1 ),0);
      d2 = yDirs[j];
      map[GetIndexDir(d2)] = PropTransform(0,(d2 == Q_S ? 1: -1 ));
      map[GetIndexDir(d2,d1)] = PropTransform((d1 == Q_W ? 1: -1 ),(d2 == Q_S ? 1: -1 ));
      prop_maps[GetIndexBound(Q_B,d2,d1)] = map;
    }
  }

  // propagation on horizontal edges
  for (int i = 0; i < 2; i++) { // E, W
    d1 = xDirs[i];
    // fill map to default offsets
    InitPropMap(map);
    map[GetIndexDir(d1)] = PropTransform((d1 == Q_W ? 1: -1 ),0);
    for (int j = 0; j < 2; j++) { // N, S
      d2 = yDirs[j];
      map[GetIndexDir(d2)] = PropTransform(0,(d2 == Q_S ? 1: -1 ));
      map[GetIndexDir(d2,d1)] = PropTransform((d1 == Q_W ? 1: -1 ),(d2 == Q_S ? 1: -1 ));
    }
    prop_maps[GetIndexBound(Q_B,d1)] = map;
  }

  //propagation on vertical edges
  for (int j = 0; j < 2; j++) { // N, S
    d1 = yDirs[j];
    // fill map to default offsets
    InitPropMap(map);
    map[GetIndexDir(d1)] = PropTransform(0,(d1 == Q_S ? 1: -1 ));
    for (int i = 0; i < 2; i++) { // E, W
      d2 = xDirs[i];
      map[GetIndexDir(d2)] = PropTransform((d2 == Q_W ? 1: -1 ),0);
      map[GetIndexDir(d1,d2)] = PropTransform((d2 == Q_W ? 1: -1 ),(d1 == Q_S ? 1: -1 ));
    }
    prop_maps[GetIndexBound(Q_B,d1)] = map;
  }

  InitPropMap(map);
  //propagation inside domain
  for (int i = 0; i < 2; i++) {
    d1 = xDirs[i]; // E, W
    d2 = yDirs[i]; // N, S
    map[GetIndexDir(d1)] = PropTransform((d1 == Q_W ? 1: -1 ),0);
    velocityDirections[GetIndexDir(d1)] = PropTransform((d1 == Q_E ? 1: -1 ),0);
    map[GetIndexDir(d2)] = PropTransform(0,(d2 == Q_S ? 1: -1 ));
    velocityDirections[GetIndexDir(d2)] = PropTransform(0,(d2 == Q_N ? 1: -1 ));
    for (int j = 0; j < 2; j++) {
      d2 = yDirs[j];
      map[GetIndexDir(d2,d1)] = PropTransform((d1 == Q_W ? 1: -1 ),(d2 == Q_S ? 1: -1 ));
      velocityDirections[GetIndexDir(d2,d1)] = PropTransform((d1 == Q_E ? 1: -1 ),(d2 == Q_N ? 1: -1 ));
    }
  }
  prop_maps[INTERIOR] = map;
}

void LatticeBoltzmann::SetupTransformation3D()
{
  // in 3D: 27 maps (8 corners + 12 edges + 6 faces + 1 map for interior elements)
  StdVector<PropTransform> map;
  map.Resize(n_q_);
  Direction d1, d2, d3;
  const Direction dirs1[] = {Q_T, Q_B};
  const Direction dirs2[] = {Q_N, Q_S};
  const Direction dirs3[] = {Q_E, Q_W};

  // propagation in corners
  for (int i = 0; i < 2; i++) { // T, B
    for (int j = 0; j < 2; j++) { // N, S
      for (int k = 0; k < 2; k++) { // E, W
        // fill map to default offsets
        InitPropMap(map);
        d1 = dirs1[i];
        d2 = dirs2[j];
        d3 = dirs3[k];

        // d1 can only be top or bottom
        map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q_B ? 1 : -1));
        // d2 can on only be north or south
        map[GetIndexDir(d2)] = PropTransform(0, (d2 == Q_S ? 1 : -1) , 0);
        // d3 can on only be east or west
        map[GetIndexDir(d3)] = PropTransform((d3 == Q_W ? 1 : -1), 0, 0);

        map[GetIndexDir(d1,d2)] = PropTransform(0, (d2 == Q_S ? 1 : -1) , (d1 == Q_B ? 1 : -1));
        map[GetIndexDir(d1,d3)] = PropTransform((d3 == Q_W ? 1 : -1), 0, (d1 == Q_B ? 1 : -1));
        map[GetIndexDir(d2,d3)] = PropTransform((d3 == Q_W ? 1 : -1), (d2 == Q_S ? 1 : -1), 0);

        prop_maps[GetIndexBound(d1,d2,d3)] = map;

      }
    }
  }

  // horizontal edges in x-y plane
  for (int i = 0; i < 2; i++) { // T, B
    for (int j = 0; j < 2; j++) { // N, S
      // fill map to default offsets
      InitPropMap(map);
      d1 = dirs1[i];
      d2 = dirs2[j];
      // d1 can only be top or bottom
      map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q_B ? 1 : -1));
      // d2 can on only be north or south
      map[GetIndexDir(d2)] = PropTransform(0, (d2 == Q_S ? 1 : -1) , 0);
      map[GetIndexDir(d1,d2)] = PropTransform(0, (d2 == Q_S ? 1 : -1) , (d1 == Q_B ? 1 : -1));
      for (int k = 0; k < 2; k++) { //E, W
        d3 = dirs3[k];
        map[GetIndexDir(d1,d3)] = PropTransform((d3 == Q_W ? 1 : -1), 0, (d1 == Q_B ? 1 : -1));
        map[GetIndexDir(d2,d3)] = PropTransform((d3 == Q_W ? 1 : -1), (d2 == Q_S ? 1 : -1), 0);
        // d3 can on only be east or west
        map[GetIndexDir(d3)] = PropTransform((d3 == Q_W ? 1 : -1), 0, 0);
      }
      prop_maps[GetIndexBound(d1,d2)] = map;
    }
  }

  // vertical edges in x-y plane
  for (int i = 0; i < 2; i++) { // T, B
    for (int j = 0; j < 2; j++) { // E, W
      // fill map to default offsets
      InitPropMap(map);
      d1 = dirs1[i];
      d2 = dirs3[j];
      // d1 can only be top or bottom
      map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q_B ? 1 : -1));
      // d2 can on only be east or west
      map[GetIndexDir(d2)] = PropTransform((d2 == Q_W ? 1 : -1), 0, 0);
      map[GetIndexDir(d1,d2)] = PropTransform((d2 == Q_W ? 1 : -1), 0, (d1 == Q_B ? 1 : -1));
      for (int k = 0; k < 2; k++) { //N, S
        d3 = dirs2[k];
        // d3 can on only be north or south
        map[GetIndexDir(d3)] = PropTransform(0, (d3 == Q_S ? 1 : -1) , 0);
        map[GetIndexDir(d1,d3)] = PropTransform(0,(d3 == Q_S ? 1 : -1),(d1 == Q_B ? 1 : -1));
        map[GetIndexDir(d3,d2)] = PropTransform((d2 == Q_W ? 1 : -1),(d3 == Q_S ? 1 : -1),0);
      }
      prop_maps[GetIndexBound(d1,d2)] = map;
    }
  }

  //vertical edges in x-z plane
  for (int i = 0; i < 2; i++) { // N, S
    for (int j = 0; j < 2; j++) { // E, W
      // fill map to default offsets
      InitPropMap(map);
      d1 = dirs2[i];
      d2 = dirs3[j];

      // d1 can on only be north or south
      map[GetIndexDir(d1)] = PropTransform(0, (d1 == Q_S ? 1 : -1) , 0);
      // d2 can on only be east or west
      map[GetIndexDir(d2)] = PropTransform((d2 == Q_W ? 1 : -1), 0, 0);

      map[GetIndexDir(d1,d2)] = PropTransform((d2 == Q_W ? 1 : -1),(d1 == Q_S ? 1 : -1),0);

      for (int k = 0; k < 2; k++) { // T, B
        d3 = dirs1[k];
        // d3 can only be top or bottom
        map[GetIndexDir(d3)] = PropTransform(0, 0, (d3 == Q_B ? 1 : -1));
        map[GetIndexDir(d3,d1)] = PropTransform(0,(d1 == Q_S ? 1 : -1),(d3 == Q_B ? 1 : -1));
        map[GetIndexDir(d3,d2)] = PropTransform((d2 == Q_W ? 1 : -1),0,(d3 == Q_B ? 1 : -1));
      }
      prop_maps[GetIndexBound(d1,d2)] = map;
    }
  }

  // propagation for inner nodes and information about microscopic velocities
  // fill map to default offsets
  InitPropMap(map);
  for (int i = 0; i < 2; i++) { // T, B
    d1 = dirs1[i];
    // d1 can only be top or bottom
    map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q_B ? 1 : -1));
    // offsets of microscopic velocities set by reversing offsets of propagation of interior elements
    velocityDirections[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q_T ? 1 : -1));
    for (int j = 0; j < 2; j++) { // N, S
      d2 = dirs2[j];
      // d2 can on only be north or south
      map[GetIndexDir(d2)] = PropTransform(0, (d2 == Q_S ? 1 : -1) , 0);
      velocityDirections[GetIndexDir(d2)] = PropTransform(0, (d2 == Q_N ? 1 : -1) , 0);

      map[GetIndexDir(d1,d2)] = PropTransform(0,(d2 == Q_S ? 1 : -1),(d1 == Q_B ? 1 : -1));
      velocityDirections[GetIndexDir(d1,d2)] = PropTransform(0,(d2 == Q_N ? 1 : -1),(d1 == Q_T ? 1 : -1));
      for (int k = 0; k < 2; k++) { // E, W
        d3 = dirs3[k];

        // d3 can on only be east or west
        map[GetIndexDir(d3)] = PropTransform((d3 == Q_W ? 1 : -1), 0, 0);
        velocityDirections[GetIndexDir(d3)] = PropTransform((d3 == Q_E ? 1 : -1), 0, 0);

        map[GetIndexDir(d1,d3)] = PropTransform((d3 == Q_W ? 1 : -1),0,(d1 == Q_B ? 1 : -1));
        map[GetIndexDir(d2,d3)] = PropTransform((d3 == Q_W ? 1 : -1),(d2 == Q_S ? 1 : -1),0);

        velocityDirections[GetIndexDir(d1,d3)] = PropTransform((d3 == Q_E ? 1 : -1),0,(d1 == Q_T ? 1 : -1));
        velocityDirections[GetIndexDir(d2,d3)] = PropTransform((d3 == Q_E ? 1 : -1),(d2 == Q_N ? 1 : -1),0);

      }
    }
  }
  prop_maps[INTERIOR] = map;

  // we start with the map for the interior elements and delete some entries for the corresponding cube face
  map = prop_maps[INTERIOR];
  PropTransform transform;
  int invDir;
  for (int i = 0; i < 2; i++) { // T,B
    d1 = dirs1[i];
    invDir = GetInvDirection(d1);
    transform = map[invDir];
    map = prop_maps[INTERIOR];
    for (int dir = 0; dir < n_q_; dir++) {
      if(map[dir].off_z != 0 && map[dir].off_z == transform.off_z)
        map[dir] = PropTransform();
    }
    prop_maps[GetIndexBound(d1)] = map;
  }

  for (int i = 0; i < 2; i++) { // N, S
    d2 = dirs2[i];
    invDir = GetInvDirection(d2);
    transform = map[invDir];
    map = prop_maps[INTERIOR];
    for (int dir = 0; dir < n_q_; dir++) {
      if(map[dir].off_y != 0 && map[dir].off_y == transform.off_y)
        map[dir] = PropTransform();
    }
    prop_maps[GetIndexBound(d2)] = map;
  }

  for (int i = 0; i < 2; i++) { // E, W
    d3 = dirs3[i];
    invDir = GetInvDirection(d3);
    transform = map[invDir];
    map = prop_maps[INTERIOR];
    for (int dir = 0; dir < n_q_; dir++) {
      if(map[dir].off_x != 0 && map[dir].off_x == transform.off_x)
        map[dir] = PropTransform();
    }
    prop_maps[GetIndexBound(d3)] = map;
  }
}

// for debugging purposes
void LatticeBoltzmann::WritePropMap()
{
  StdVector<Boundary> corners;
  corners.Push_back(CORNER_BSW);
  corners.Push_back(CORNER_BSE);
  corners.Push_back(CORNER_BNW);
  corners.Push_back(CORNER_BNE);

  StdVector<Direction> order;
  order.Push_back(Q_0);
  order.Push_back(Q_S);
  order.Push_back(Q_SE);
  order.Push_back(Q_SW);
  order.Push_back(Q_N);
  order.Push_back(Q_NE);
  order.Push_back(Q_NW);
  order.Push_back(Q_W);
  order.Push_back(Q_E);

  int n_corners = 4;

  if (m_dim == 3) {
    corners.Push_back(CORNER_TNW);
    corners.Push_back(CORNER_TSW);
    corners.Push_back(CORNER_TNE);
    corners.Push_back(CORNER_TSE);
    order.Clear();
    order.Push_back(Q_0);
    order.Push_back(Q_E);
    order.Push_back(Q_W);
    order.Push_back(Q_N);
    order.Push_back(Q_S);
    order.Push_back(Q_T);
    order.Push_back(Q_B);
    order.Push_back(Q_NE);
    order.Push_back(Q_SW);
    order.Push_back(Q_NW);
    order.Push_back(Q_SE);
    order.Push_back(Q_TN);
    order.Push_back(Q_BS);
    order.Push_back(Q_TS);
    order.Push_back(Q_BN);
    order.Push_back(Q_TE);
    order.Push_back(Q_BW);
    order.Push_back(Q_TW);
    order.Push_back(Q_BE);
    n_corners = 8;
  }

  fstream f;
  f.precision(16);

  // write out propagation maps for corners
  f.open("propagation_corners.txt", ios::out);
  for (int i = 0; i < n_corners; i++) {
    StdVector<PropTransform>& mymap  = prop_maps[corners[i]];
    f << "// " << boundaries.ToString(corners[i]) << std::endl;
    for (int j = 0; j < n_q_; j++) {
      f << "pdf(cur, x, y, z, " << directions.ToString(order[j]) << ") = pdf(cur, x ";
      if (mymap[order[j]].off_x == 1)
        f << " + 1";
      else if (mymap[order[j]].off_x == -1)
        f <<  " - 1";
      f << ", y";
      if (mymap[order[j]].off_y == 1)
        f << " + 1";
      else if (mymap[order[j]].off_y == -1)
        f << " - 1";
      f << ", z";
      if (mymap[order[j]].off_z == 1)
        f << " + 1";
      else if(mymap[order[j]].off_z == -1)
        f << " - 1";

      f << ", " << directions.ToString(order[j]) << ") \n";
    }
    f << std::endl;
  }

  // write out propagation maps for interior nodes
  StdVector<PropTransform>& mymap  = prop_maps[INTERIOR];
  f << "// " << boundaries.ToString(INTERIOR) << std::endl;
  for (int j = 0; j < n_q_; j++) {
    f << "pdf(cur, x, y, z, " << directions.ToString(order[j]) << ") = pdf(cur, x ";
    if (mymap[order[j]].off_x == 1)
      f << " + 1";
    else if (mymap[order[j]].off_x == -1)
      f <<  " - 1";
    f << ", y";
    if (mymap[order[j]].off_y == 1)
      f << " + 1";
    else if (mymap[order[j]].off_y == -1)
      f << " - 1";
    f << ", z";
    if (mymap[order[j]].off_z == 1)
      f << " + 1";
    else if(mymap[order[j]].off_z == -1)
      f << " - 1";

    f << ", " << directions.ToString(order[j]) << ") \n";
  }
  f << std::endl;
  f.close();

  f.open("micro_velocities.txt", ios::out);
  for (int j = 0; j < n_q_; j++) {
    f << "pdf(cur, x, y, z, " << directions.ToString(order[j]) << ") = pdf(cur, x + " << velocityDirections[order[j]].off_x << ", y +" << velocityDirections[order[j]].off_y << ", z + "
        << velocityDirections[order[j]].off_z << ", " << directions.ToString(order[j]) << ") \n";
  }
  f << std::endl;
  f.close();

  // write out propagation maps for edges
  StdVector<Boundary> edges;
  edges.Push_back(EDGE_BS);
  edges.Push_back(EDGE_BN);
  edges.Push_back(EDGE_BW);
  edges.Push_back(EDGE_BE);

  int n_edges = 4;
  if (m_dim == 3) {
    edges.Push_back(EDGE_TN);
    edges.Push_back(EDGE_TW);
    edges.Push_back(EDGE_TS);
    edges.Push_back(EDGE_TE);
    edges.Push_back(EDGE_NW);
    edges.Push_back(EDGE_SW);
    edges.Push_back(EDGE_SE);
    edges.Push_back(EDGE_NE);
    n_edges = 12;
  }
  f.open("propagation_edges.txt", ios::out);
  for (int i = 0; i < n_edges; i++) {
    StdVector<PropTransform>& mymap  = prop_maps[edges[i]];
    f << "// " << boundaries.ToString(edges[i]) << std::endl;
    for (int j = 0; j < n_q_; j++) {
      f << "pdf(cur, x, y, z, " << directions.ToString(order[j]) << ") = pdf(cur, x";
      if (mymap[order[j]].off_x == 1)
        f << " + 1";
      else if (mymap[order[j]].off_x == -1)
        f <<  " - 1";
      f << ", y";
      if (mymap[order[j]].off_y == 1)
        f << " + 1";
      else if (mymap[order[j]].off_y == -1)
        f << " - 1";
      f << ", z";
      if (mymap[order[j]].off_z == 1)
        f << " + 1";
      else if(mymap[order[j]].off_z == -1)
        f << " - 1";

      f << ", " << directions.ToString(order[j]) << ") \n";
    }
    f << std::endl;
  }
  f.close();

  if (m_dim == 3)
  {
    // write out propagation maps for faces
    const Boundary faces[] = {FACE_T, FACE_B, FACE_S, FACE_N, FACE_W, FACE_E};
    f.open("propagation_faces.txt", ios::out);
    for (int i = 0; i < 6; i++) {
      StdVector<PropTransform>& mymap  = prop_maps[faces[i]];
      f << "// " << boundaries.ToString(faces[i]) << std::endl;
      for (int j = 0; j < n_q_; j++) {
        f << "pdf(cur, x, y, z, " << directions.ToString(order[j]) << ") = pdf(cur, x";
        if (mymap[order[j]].off_x == 1)
          f << " + 1";
        else if (mymap[order[j]].off_x == -1)
          f <<  " - 1";
        f << ", y";
        if (mymap[order[j]].off_y == 1)
          f << " + 1";
        else if (mymap[order[j]].off_y == -1)
          f << " - 1";
        f << ", z";
        if (mymap[order[j]].off_z == 1)
          f << " + 1";
        else if(mymap[order[j]].off_z == -1)
          f << " - 1";

        f << ", " << directions.ToString(order[j]) << ") \n";
      }
      f << std::endl;
    }
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

  boundaries.SetName("Corners and edges of a cube");
  boundaries.Add(INTERIOR,"interior");
  boundaries.Add(CORNER_BNW, "left bottom back corner");
  boundaries.Add(CORNER_BSW, "left bottom front corner");
  boundaries.Add(CORNER_BSE, "right bottom front corner");
  boundaries.Add(CORNER_BNE, "right bottom back corner");
  boundaries.Add(CORNER_TNW, "left top back corner");
  boundaries.Add(CORNER_TSW, "left top front corner");
  boundaries.Add(CORNER_TSE, "right top front corner");
  boundaries.Add(CORNER_TNE, "right top back corner");
  boundaries.Add(EDGE_TN ,"top back edge");
  boundaries.Add(EDGE_TW ,"top left edge");
  boundaries.Add(EDGE_TS ,"top front edge");
  boundaries.Add(EDGE_TE ,"top right edge");
  boundaries.Add(EDGE_NW ,"left back edge");
  boundaries.Add(EDGE_SW ,"left front edge");
  boundaries.Add(EDGE_SE ,"right front edge");
  boundaries.Add(EDGE_NE ,"right back edge");
  boundaries.Add(EDGE_BN ,"bottom back edge");
  boundaries.Add(EDGE_BW ,"bottom left edge");
  boundaries.Add(EDGE_BS ,"bottom front edge");
  boundaries.Add(EDGE_BE ,"bottom right edge");
  boundaries.Add(FACE_W,"left face");
  boundaries.Add(FACE_E,"right face");
  boundaries.Add(FACE_T,"top face");
  boundaries.Add(FACE_B,"bottom face");
  boundaries.Add(FACE_S,"front face");
  boundaries.Add(FACE_N,"back face");
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
  directionsInv.Resize(n_q_);

  directionsInv[Q_0] = Q_0;
  directionsInv[Q_E] = Q_W;
  directionsInv[Q_N] = Q_S;
  directionsInv[Q_W] = Q_E;
  directionsInv[Q_S] = Q_N;
  directionsInv[Q_NE] = Q_SW;
  directionsInv[Q_SW] = Q_NE;
  directionsInv[Q_NW] = Q_SE;
  directionsInv[Q_SE] = Q_NW;

  if (m_dim == 3) {
    directionsInv[Q_T] = Q_B;
    directionsInv[Q_B] = Q_T;
    directionsInv[Q_TN] = Q_BS;
    directionsInv[Q_BS] = Q_TN;
    directionsInv[Q_TS] = Q_BN;
    directionsInv[Q_BN] = Q_TS;
    directionsInv[Q_TE] = Q_BW;
    directionsInv[Q_BW] = Q_TE;
    directionsInv[Q_TW] = Q_BE;
    directionsInv[Q_BE] = Q_TW;
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

//fills given propagation map with given default offset values (0 in all directions
inline void LatticeBoltzmann::InitPropMap(StdVector<PropTransform>& map)
{
  map.Init(PropTransform());
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

  // perform a propagation step
  int x, y, z = 0;
  const Boundary corners[] = {CORNER_BNW, CORNER_BSW, CORNER_BNE, CORNER_BSE};

  StdVector<PropTransform> *map1, *map2;
  PropTransform *transform1, *transform2;
  Boundary bound;

  int n_corners = 4;

  // ---------------------------------------------------------------------
  // Propagation for nodes in the corners....
  // ---------------------------------------------------------------------
  // propagation maps for corners stored in prop_maps from 1 to 4
  for (int i = 0; i < n_corners; i++) {
    bound = corners[i];
    switch(bound)
    {
      case CORNER_BNW: x = 0; y = m_sizeY - 1 ; break;
      case CORNER_BSW: x = 0; y = 0; break;
      case CORNER_BNE: x = m_sizeX -1; y = m_sizeY - 1; break;
      case CORNER_BSE: x = m_sizeX - 1; y = 0; break;
      default:
        EXCEPTION("Trying to propagate into a corner that does not exist!");
    }
    map1 = &prop_maps[bound];

    for (int dir = 0; dir < n_q_; dir++) {
      transform1 = &(*map1)[dir];
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform1->off_x, y + transform1->off_y, z, dir);
    }
  }

  // Propagation along the boundaries of the X dimension.
  // lower and upper edge
  map1 = &prop_maps[EDGE_BS];
  map2 = &prop_maps[EDGE_BN];

  for (int dir = 0; dir < n_q_; dir++) {
    transform1 = &(*map1)[dir];
    transform2 = &(*map2)[dir];
    for (x = 1; x < m_sizeX - 1; ++x) {
      //BS
      y = 0;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform1->off_x, y + transform1->off_y, z, dir);

      //BN
      y = m_sizeY - 1;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform2->off_x, y + transform2->off_y, z, dir);
    }
  }

  // Propagation along the boundaries of the Y dimension.
  // edges: left and right edge
  map1 = &prop_maps[EDGE_BW];
  map2 = &prop_maps[EDGE_BE];

  for (int dir = 0; dir < n_q_; dir++) {
    transform1 = &(*map1)[dir];
    transform2 = &(*map2)[dir];
    for (y = 1; y < m_sizeY - 1; ++y) {
      //BW: left edge
      x = 0;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform1->off_x, y + transform1->off_y, z, dir);

      //BE: right edge
      x = m_sizeX - 1;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform2->off_x, y + transform2->off_y, z, dir);
    }
  }

  //--------------------- propagation over inner elements--------------------------------------------------------------------

  const StdVector<PropTransform>& map_interior = prop_maps[INTERIOR];
  double tmp, tmp_ux, tmp_uy, tmp_us, scale, sum;
  double * scales  = Scales.GetPointer();

  int index;

#pragma omp parallel default(none)\
    private(index), \
    private(tmp_ux, tmp_uy, tmp_us, scale, sum, tmp, x, y), \
    shared(z, next, cur, map_interior, scales, omega)
{
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  #pragma omp for collapse(2)
  for (y = 1; y < m_sizeY - 1; ++y) {
    for (x = 1; x < m_sizeX - 1; ++x) {

      index= GetIndex(x,y,z);

      // sum: macroscopic density is sum over all discrete distributions of an element
      sum = 0;
      tmp_ux = 0;
      tmp_uy = 0;
      tmp_us = 0;

      for (int dir = 0; dir < n_q_; dir++) {
        //store current pdf values in array for better accessing
        pdfs[dir] = pdf(cur, x + map_interior[dir].off_x, y + map_interior[dir].off_y,  z, dir);
        sum += pdfs[dir];
        tmp_ux += velocityDirections[dir].off_x*pdfs[dir];
        tmp_uy += velocityDirections[dir].off_y*pdfs[dir];
      }

      // macroscopic scaling by design variable
      scale = scales[index];

      tmp_ux = scale * tmp_ux / sum;
      tmp_uy = scale * tmp_uy / sum;
      tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);

      tmp_ux = 3.0 * tmp_ux;
      tmp_uy = 3.0 * tmp_uy;

      for (int dir = 0; dir < n_q_; dir++) {
        tmp = velocityDirections[dir].off_x * tmp_ux + velocityDirections[dir].off_y * tmp_uy ;
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
      tmp = velocityDirections[dir].off_x * tmp_ux + velocityDirections[dir].off_y * tmp_uy;
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
      tmp = velocityDirections[dir].off_x * tmp_ux + velocityDirections[dir].off_y * tmp_uy;
      pdf(cur, x, y, z, dir) = sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


/************************************************** 3D operators *****************************************************/
void LatticeBoltzmann::prop_coll_step3D(int cur, int next, double omega)
{

  // perform a propagation step
  int x, y, z;
  const Boundary corners[] = {CORNER_BNW, CORNER_BSW, CORNER_BNE, CORNER_BSE, CORNER_TNW, CORNER_TSW, CORNER_TNE, CORNER_TSE};

  StdVector<PropTransform> *map1, *map2;
  PropTransform *transform1, *transform2;
  Boundary bound;

  // only necessary in 3D
  StdVector<PropTransform> *map3, *map4;
  PropTransform *transform3, *transform4;

  int n_corners = 8;

  // ---------------------------------------------------------------------
  // Propagation for nodes in the corners....
  // ---------------------------------------------------------------------
  // propagation maps for corners stored in prop_maps from 1 to 8
  for (int i = 0; i < n_corners; i++) {
    bound = corners[i];
    switch(bound)
    {
      case CORNER_BNW: x = 0; y = m_sizeY - 1 ; z = 0; break;
      case CORNER_BSW: x = 0; y = 0; z = 0; break;
      case CORNER_BNE: x = m_sizeX -1; y = m_sizeY - 1; z = 0; break;
      case CORNER_BSE: x = m_sizeX - 1, y = 0, z = 0; break;
      case CORNER_TNW: x = 0; y = m_sizeY - 1; z = m_sizeZ - 1; break;
      case CORNER_TSW: x = 0; y = 0; z = m_sizeZ - 1; break;
      case CORNER_TNE: x = m_sizeX - 1; y = m_sizeY - 1; z = m_sizeZ - 1; break;
      case CORNER_TSE: x = m_sizeX - 1; y = 0; z = m_sizeZ - 1; break;
      default:
        EXCEPTION("Trying to propagate into a corner that does not exist!");
    }
    map1 = &prop_maps[bound];

    for (int dir = 0; dir < n_q_; dir++) {
      transform1 = &(*map1)[dir];
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform1->off_x, y + transform1->off_y, z + transform1->off_z, dir);
    }
  }

  // Propagation along the boundaries of the X dimension.
  // edges: BS, BN, TS, TN
  map1 = &prop_maps[EDGE_BS];
  map2 = &prop_maps[EDGE_BN];

  map3 = &prop_maps[EDGE_TS];
  map4 = &prop_maps[EDGE_TN];

  for (int dir = 0; dir < n_q_; dir++) {
    transform1 = &(*map1)[dir];
    transform2 = &(*map2)[dir];
    for (x = 1; x < m_sizeX - 1; ++x) {
      //BS: bottom front
      y = 0; z = 0;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform1->off_x, y + transform1->off_y, z + transform1->off_z, dir);

      //BN: bottom back
      y = m_sizeY - 1; z = 0;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform2->off_x, y + transform2->off_y, z + transform2->off_z, dir);

      transform3 = &(*map3)[dir];
      transform4 = &(*map4)[dir];
      //TS: top front
      y = 0; z = m_sizeZ - 1;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform3->off_x, y + transform3->off_y, z + transform3->off_z, dir);

      //TN: top back
      y = m_sizeY - 1; z = m_sizeZ - 1;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform4->off_x, y + transform4->off_y, z + transform4->off_z, dir);
    }
  }

  // Propagation along the boundaries of the Y dimension.
  // edges: TW, TE, BW, BE
  //required in 2D and 3D
  map1 = &prop_maps[EDGE_BW];
  map2 = &prop_maps[EDGE_BE];

  //additional in 3D
  map3 = &prop_maps[EDGE_TW];
  map4 = &prop_maps[EDGE_TE];

  for (int dir = 0; dir < n_q_; dir++) {
    transform1 = &(*map1)[dir];
    transform2 = &(*map2)[dir];
    for (y = 1; y < m_sizeY - 1; ++y) {
      //BW: bottom left
      x = 0; z = 0;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform1->off_x, y + transform1->off_y, z + transform1->off_z, dir);

      //BE: bottom right
      x = m_sizeX - 1; z = 0;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform2->off_x, y + transform2->off_y, z + transform2->off_z, dir);

      transform3 = &(*map3)[dir];
      transform4 = &(*map4)[dir];
      //TW: top left
      x = 0; z = m_sizeZ - 1;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform3->off_x, y + transform3->off_y, z + transform3->off_z, dir);

      //TE: top right
      x = m_sizeX - 1; z = m_sizeZ - 1;
      pdf(next, x, y, z, dir)  = pdf(cur, x + transform4->off_x, y + transform4->off_y, z + transform4->off_z, dir);
    }
  }

  // Propagation along the boundaries of the Z dimension.
  // edges: SW, NW, SE, NE
  map1 = &prop_maps[EDGE_SW];
  map2 = &prop_maps[EDGE_NW];
  map3 = &prop_maps[EDGE_SE];
  map4 = &prop_maps[EDGE_NE];
  for (z = 1; z < m_sizeZ - 1; ++z) {
    for (int i = 0; i < n_q_; i++) {
      transform1 = &(*map1)[i];
      transform2 = &(*map2)[i];
      transform3 = &(*map3)[i];
      transform4 = &(*map4)[i];

      // SW: left front
      x = 0; y = 0;
      pdf(next, x, y, z, i)  = pdf(cur, x + transform1->off_x, y + transform1->off_y, z + transform1->off_z, i);

      // NW: left back
      x = 0; y = m_sizeY - 1;
      pdf(next, x, y, z, i)  = pdf(cur, x + transform2->off_x, y + transform2->off_y, z + transform2->off_z, i);

      // SE: right front
      x = m_sizeX - 1; y = 0;
      pdf(next, x, y, z, i)  = pdf(cur, x + transform3->off_x, y + transform3->off_y, z + transform3->off_z, i);

      // NE: right back
      x = m_sizeX - 1; y = m_sizeY - 1;
      pdf(next, x, y, z, i)  = pdf(cur, x + transform4->off_x, y + transform4->off_y, z + transform4->off_z, i);
    }
  }

  //propagation along top and bottom face
  map1 = &prop_maps[FACE_T];
  map2 = &prop_maps[FACE_B];
  for (int dir = 0; dir < n_q_; dir++) {
    transform1 = &(*map1)[dir];
    transform2 = &(*map2)[dir];
    for (y = 1; y < m_sizeY - 1; y++) {
      for (x = 1; x < m_sizeX - 1; x++) {

        z = m_sizeZ - 1; // top face
        pdf(next, x, y, z, dir) = pdf(cur, x + transform1->off_x, y + transform1->off_y, z + transform1->off_z, dir);

        z = 0; // bottom face
        pdf(next, x, y, z, dir) = pdf(cur, x + transform2->off_x, y + transform2->off_y, z + transform2->off_z, dir);

      }
    }
  }

  //propagation along front and back face
  map1 = &prop_maps[FACE_N];
  map2 = &prop_maps[FACE_S];
  for (int dir = 0; dir < n_q_; dir++) {
    transform1 = &(*map1)[dir];
    transform2 = &(*map2)[dir];
    for (z = 1; z < m_sizeZ - 1; z++) {
      for (x = 1; x < m_sizeX - 1; x++) {

        y = m_sizeY - 1; // back face
        pdf(next, x, y, z, dir) = pdf(cur, x + transform1->off_x, y + transform1->off_y, z + transform1->off_z, dir);

        y = 0; // front face
        pdf(next, x, y, z, dir) = pdf(cur, x + transform2->off_x, y + transform2->off_y, z + transform2->off_z, dir);
      }
    }
  }

  //propagation along left and right face
  map1 = &prop_maps[FACE_W];
  map2 = &prop_maps[FACE_E];
  for (int dir = 0; dir < n_q_; dir++) {
    transform1 = &(*map1)[dir];
    transform2 = &(*map2)[dir];
    for (z = 1; z < m_sizeZ - 1; z++) {
      for (y = 1; y < m_sizeY - 1; y++) {
        x = 0; // left face
        pdf(next, x, y, z, dir) = pdf(cur, x + transform1->off_x, y + transform1->off_y, z + transform1->off_z, dir);

        x = m_sizeX - 1; // right face
        pdf(next, x, y, z, dir) = pdf(cur, x + transform2->off_x, y + transform2->off_y, z + transform2->off_z, dir);
      }
    }
  }


  //--------------------- propagation over inner elements--------------------------------------------------------------------

  const StdVector<PropTransform>& map_interior = prop_maps[INTERIOR];
  double tmp, tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum;
  double * scales  = Scales.GetPointer();

  int index;

#pragma omp parallel default(none)\
    private(index), \
    private(tmp_ux, tmp_uy, tmp_uz, tmp_us, scale, sum, tmp, x, y,z), \
    shared(next, cur, map_interior, scales, omega)
{
  StdVector<double> pdfs;
  pdfs.Resize(n_q_);
  #pragma omp for collapse(3)
  for (z = 1; z < m_sizeZ - 1; ++z) {
    for (y = 1; y < m_sizeY - 1; ++y) {
      for (x = 1; x < m_sizeX - 1; ++x) {
        index= GetIndex(x,y,z);

        // sum: macroscopic density is sum over all discrete distributions of an element
        sum = 0;
        tmp_ux = 0;
        tmp_uy = 0;
        tmp_uz = 0;
        tmp_us = 0;

        for (int dir = 0; dir < n_q_; dir++) {
          //store current pdf values in array for better accessing
          pdfs[dir] = pdf(cur, x + map_interior[dir].off_x, y + map_interior[dir].off_y,  z + map_interior[dir].off_z, dir);
          sum += pdfs[dir];
          tmp_ux += velocityDirections[dir].off_x*pdfs[dir];
          tmp_uy += velocityDirections[dir].off_y*pdfs[dir];
          tmp_uz += velocityDirections[dir].off_z*pdfs[dir];
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
          tmp = velocityDirections[dir].off_x * tmp_ux + velocityDirections[dir].off_y * tmp_uy + velocityDirections[dir].off_z * tmp_uz;
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
      tmp = velocityDirections[dir].off_x * tmp_ux + velocityDirections[dir].off_y * tmp_uy + velocityDirections[dir].off_z * tmp_uz;
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
      tmp = velocityDirections[dir].off_x * tmp_ux + velocityDirections[dir].off_y * tmp_uy + velocityDirections[dir].off_z * tmp_uz;
      pdf(cur, x, y, z, dir) = sum * weights[dir] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}

} // end namespace
