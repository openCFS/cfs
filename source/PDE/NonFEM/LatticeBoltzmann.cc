#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Driver/basedriver.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/NonFEM/LatticeBoltzmann.hh"
#include "Utils/Timer.hh"


namespace CoupledField
{

using std::fstream;
using std::ios;

DECLARE_LOG(lattice)
DEFINE_LOG(lattice, "lattice")

// instantiation of the static elements
Enum<LatticeBoltzmann::Boundary>        LatticeBoltzmann::boundaries;
Enum<LatticeBoltzmann::Direction>        LatticeBoltzmann::directions;

 LatticeBoltzmann::LatticeBoltzmann(int sizeX, int sizeY, double ux, double uy, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency)
 {
   m_sizeX = sizeX;
   m_sizeY = sizeY;
   m_ux = ux;
   m_uy = uy;
   m_omega = omega;
   m_maxIter = maxIterations;
   m_maxTol = maxTolerance;

   m_nNodes = m_sizeX * m_sizeY;

   m_plot = plot;

   m_writeFrequency = writeFrequency;
   m_numWriteResults = 0;

   LOG_DBG(lattice) << "Write frequency for hdf5 file: " << m_writeFrequency;

   //matrix of the probability distributions
   LOG_DBG(lattice) << "Allocating arrays for " << m_nNodes << " PDFs (" << (sizeof(double) * m_nNodes * N_DIR * 2.0 / 1024.0 / 1024.0) << " MiB)";

   m_pdfs.Resize(2);
   m_pdfs[0].Resize(m_nNodes * N_DIR);
   m_pdfs[1].Resize(m_nNodes * N_DIR);

   prop_maps.Resize(N_DIR);
   for (int i = 0; i < N_DIR; i++ ){
     prop_maps[i].Resize(N_DIR);
   }

   SetEnums();
   InitWeights();
//   TestDirectionIndex();
//   testInvDirections();
   SetupTransformation();

   // for debugging
//   WritePropMap();

   Scales.Resize(m_nNodes);
   rel.Resize(m_nNodes);
   bb.Resize(2 * m_sizeX + 2 * m_sizeY);

   m_cur  = 0;
   m_next = 1;
 }

 LatticeBoltzmann::~LatticeBoltzmann()
 {
 }

StdVector<double>* LatticeBoltzmann::Iterate(const StdVector<double>& elements, PtrParamNode in)
{
  int index;
  int it = 0;
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
  for (int i = 0; i < m_nNodes; ++i)
    Scales[i] = 1.0 - elements[i];

  Timer timer;
  timer.Start();

  // create_output("node.dat", m_cur);

  LOG_DBG(lattice) << "bb = " << ToString(bb);
  LOG_DBG(lattice) << "inlet = " << ToString(inlet);
  LOG_DBG(lattice) << "outlet = " << ToString(outlet);
  LOG_DBG(lattice) << "rel = " << ToString(rel);

  in->Get("converged")->SetValue("running");

//  domain->GetDriver()->StoreResults(0,0.0);

  while(it < m_maxIter && !steady_state && R <= 1000)
  {
    if (it % m_writeFrequency == 0) {
      domain->GetDriver()->StoreResults(count,(double)count);
      count++;
    }
    // -- Combined propagation and collision step -------------------------
    prop_coll_step(m_cur, m_next, m_omega);

    // -- Bounce back step ------------------------------------------------
    prop_coll_bounce_back(m_next, bb);

    // -- Inlet condition -------------------------------------------------
    prop_coll_velinlet(m_next, inlet, m_ux, m_uy);

    // -- Outlet condition ------------------------------------------------
    prop_coll_densoutlet(m_next, outlet);

    if((it == 0 || it % 100 == 0))
    {
      //Calculation of the residual
      R = 0.;

      // TODO: parallel computation of residual, take care of non_sing and
      //       summation of residual.
      // TG ORIG: for(int i = 0; i < m_sizeX; i++) {
      // TG ORIG:   for(int j = 0; j < m_sizeY; j++) {
      for(int j = 0; j < m_sizeY; j++) {
        for(int i = 0; i < m_sizeX; i++) {
          index = j * m_sizeX + i;

          for(int k = 0; k < N_DIR; k++)
          {
            res = PDF_IDX(m_next, index, k) - PDF_IDX(m_cur, index, k);
            R += res * res;
          }
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
    }
    m_cur  = (m_cur  + 1) % 2;
    m_next = (m_next + 1) % 2;

    it++;
  }
  timer.Stop();

  in->Get("iterations")->SetValue(it);
  in->Get("residuum")->SetValue(R);
  in->Get("converged")->SetValue(steady_state);

  if(m_plot) {
    plot << it << "\t" << R << "\n";
    plot.flush();
  }

  double wt = timer.GetWallTime();
  double performance = (m_sizeX - 1) * (m_sizeY - 1) * it / wt / 1e6;
  in->Get("MFLUP_s")->SetValue(performance);
  in->Get("walltime")->SetValue(wt);

  if(R >= 1000) 
    EXCEPTION("In LBM iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("internal LBM simulation could not converge: iterations: " << it << " residuum: " << R);

  m_numWriteResults = count;
  return &(m_pdfs[m_cur]);
}

int LatticeBoltzmann::GetNumWriteResults()
{
  return m_numWriteResults;
}


void LatticeBoltzmann::InitializePdfs()
{
  int i;

  #ifdef _OPENMP
    // Perform NUMA-correct placement.
    #pragma omp parallel for default(none), private(i), collapse(2)
  #endif
  for (int y = 0; y < m_sizeY; ++y) {
    for (int x = 0; x < m_sizeX; ++x) {

      i = y * m_sizeX + x;

      PDF_IDX(0, i, 0) = 4. / 9.;
      PDF_IDX(1, i, 0) = 4. / 9.;

      for(int j = 1; j < 5; j++) {
        PDF_IDX(0, i, j) = 1. / 9.;
        PDF_IDX(1, i, j) = 1. / 9.;
      }

      for(int j = 5; j < 9; j++) {
        PDF_IDX(0, i, j) = 1. / 36.;
        PDF_IDX(1, i, j) = 1. / 36.;
      }
    }
  }
}

inline int LatticeBoltzmann::GetIndexDir(Direction dir1)
{
//  assert(directions.IsValid(dir1));
  return dir1;
}


inline int LatticeBoltzmann::GetIndexDir(Direction dir1, Direction dir2)
{
//  assert(directions.IsValid(dir1));
  assert(dir2 != 0);
//  assert(directions.IsValid(dir1));
//  assert(directions.IsValid(dir2));
  assert((dir1 == Q9_S) || (dir1 == Q9_N));
  assert((dir2 == Q9_E) || (dir2 == Q9_W));
  return 2*dir1+1.5*dir2-0.5*dir1*dir2+0.5;
}

StdVector<double> LatticeBoltzmann::GetPdfs() {
  return m_pdfs[m_cur];
}

// validates GetIndexDir function
void LatticeBoltzmann::TestDirectionIndex()
{
  assert(GetIndexDir(Q9_0) == 0);
  assert(GetIndexDir(Q9_E) == 1);
  assert(GetIndexDir(Q9_N) == 2);
  assert(GetIndexDir(Q9_W) == 3);
  assert(GetIndexDir(Q9_S) == 4);
  assert(GetIndexDir(Q9_N,Q9_E) == 5);
  assert(GetIndexDir(Q9_N,Q9_W) == 6);
  assert(GetIndexDir(Q9_S,Q9_W) == 7);
  assert(GetIndexDir(Q9_S,Q9_E) == 8);
}

void LatticeBoltzmann::SetupTransformation()
{
  // in 2D: 8 maps (4 corners + 4 edges)
  StdVector<PropTransform> map;
  map.Resize(N_DIR);
  const Direction xDirs[] = {Q9_E,Q9_W};
  const Direction yDirs[] = {Q9_N,Q9_S};
  Direction d1, d2;

//  // propagation in corners
  for (int i = 0; i < 2; i++) { // E, W
    for (int j = 0; j < 2; j++) { // N, S
      // fill map to default offsets
      InitPropMap(map);
      d1 = xDirs[i];
      map[GetIndexDir(d1)] = PropTransform((d1 == Q9_W ? 1: -1 ),0);
      d2 = yDirs[j];
      map[GetIndexDir(d2)] = PropTransform(0,(d2 == Q9_S ? 1: -1 ));
      map[GetIndexDir(d2,d1)] = PropTransform((d1 == Q9_W ? 1: -1 ),(d2 == Q9_S ? 1: -1 ));
      prop_maps[GetIndexDir(d2,d1)] = map;
    }
  }

  // propagation on horizontal edges
  for (int i = 0; i < 2; i++) { // E, W
    d1 = xDirs[i];
    // fill map to default offsets
    InitPropMap(map);
    map[GetIndexDir(d1)] = PropTransform((d1 == Q9_W ? 1: -1 ),0);
    for (int j = 0; j < 2; j++) { // N, S
      d2 = yDirs[j];
      map[GetIndexDir(d2)] = PropTransform(0,(d2 == Q9_S ? 1: -1 ));
      map[GetIndexDir(d2,d1)] = PropTransform((d1 == Q9_W ? 1: -1 ),(d2 == Q9_S ? 1: -1 ));
    }
    prop_maps[GetIndexDir(d1)] = map;
  }

  //propagation on vertical edges
  for (int j = 0; j < 2; j++) { // N, S
    d1 = yDirs[j];
    // fill map to default offsets
    InitPropMap(map);
    map[GetIndexDir(d1)] = PropTransform(0,(d1 == Q9_S ? 1: -1 ));
    for (int i = 0; i < 2; i++) { // E, W
      d2 = xDirs[i];
      map[GetIndexDir(d2)] = PropTransform((d2 == Q9_W ? 1: -1 ),0);
      map[GetIndexDir(d1,d2)] = PropTransform((d2 == Q9_W ? 1: -1 ),(d1 == Q9_S ? 1: -1 ));
    }
    prop_maps[GetIndexDir(d1)] = map;
  }

  InitPropMap(map);
  //propagation inside domain
  for (int i = 0; i < 2; i++) {
    d1 = xDirs[i]; // E, W
    d2 = yDirs[i]; // N, S
    map[GetIndexDir(d1)] = PropTransform((d1 == Q9_W ? 1: -1 ),0);
    map[GetIndexDir(d2)] = PropTransform(0,(d2 == Q9_S ? 1: -1 ));
    for (int j = 0; j < 2; j++) {
      d2 = yDirs[j];
      map[GetIndexDir(d2,d1)] = PropTransform((d1 == Q9_W ? 1: -1 ),(d2 == Q9_S ? 1: -1 ));
    }
  }
  prop_maps[0] = map;

}

// for debugging purposes
void LatticeBoltzmann::WritePropMap()
{
  // order of output for comparison with old implementation
  Direction order[9];
  order[0] = Q9_0;
  order[1] = Q9_S;
  order[2] = Q9_SE;
  order[3] = Q9_SW;
  order[4] = Q9_N;
  order[5] = Q9_NE;
  order[6] = Q9_NW;
  order[7] = Q9_W;
  order[8] = Q9_E;

  fstream f;
  f.precision(16);
  f.open("propagation_map.txt", ios::out);
  const Boundary allBoundaries[] = {C, E, N, W, S, NE, NW, SW, SE};
  for (int i = 0; i < N_DIR; i++) {
    StdVector<PropTransform>& mymap  = prop_maps[i];
    f << boundaries.ToString(allBoundaries[i]) << std::endl;
    f << "---------------------------------------" << std::endl;
    for (int i = 0; i < N_DIR; i++) {
      f << directions.ToString(order[i]) << " \t x+" << mymap[order[i]].off_x << ", y+" << mymap[order[i]].off_y << std::endl;
    }
    f << std::endl;
  }
  f.close();
}

void LatticeBoltzmann::SetEnums()
{
  boundaries.SetName("Corners and edges in 2D");
  boundaries.Add(C,"interior");
  boundaries.Add(E,"right edge");
  boundaries.Add(N,"top edge");
  boundaries.Add(W,"left edge");
  boundaries.Add(S,"bottom edge");
  boundaries.Add(NE,"top right corner");
  boundaries.Add(NW,"top left corner");
  boundaries.Add(SW,"bottom left corner");
  boundaries.Add(SE,"bottom right corner");

  directions.SetName("Q9 directions");
  directions.Add(Q9_0,"C");
  directions.Add(Q9_E,"E");
  directions.Add(Q9_N,"N");
  directions.Add(Q9_W,"W");
  directions.Add(Q9_S,"S");
  directions.Add(Q9_NE,"NE");
  directions.Add(Q9_NW,"NW");
  directions.Add(Q9_SW,"SW");
  directions.Add(Q9_SE,"SE");
}

void LatticeBoltzmann::InitWeights()
{
  weights.Resize(N_DIR);
  weights[Q9_0] = 4.0 / 9.0;
  weights[Q9_E] = 1.0 /  9.0;
  weights[Q9_N] = 1.0 /  9.0;
  weights[Q9_W] = 1.0 /  9.0;
  weights[Q9_S] = 1.0 /  9.0;
  weights[Q9_NE] = 1.0 / 36.0;
  weights[Q9_NW] = 1.0 / 36.0;
  weights[Q9_SW] = 1.0 / 36.0;
  weights[Q9_SE] = 1.0 / 36.0;
}

// depends on numbering of directions
inline int LatticeBoltzmann::GetInvDirection(Direction dir)
{
  assert(directions.IsValid(dir));
  if (dir == Q9_0)
    return dir;
  if (dir == Q9_E || dir == Q9_N || dir == Q9_NE || dir == Q9_NW)
    return dir + 2;
  else
    return dir - 2;
}

void LatticeBoltzmann::TestInvDirections()
{
  assert(GetInvDirection(Q9_0) == Q9_0);
  assert(GetInvDirection(Q9_E) == Q9_W);
  assert(GetInvDirection(Q9_N) == Q9_S);
  assert(GetInvDirection(Q9_W) == Q9_E);
  assert(GetInvDirection(Q9_S) == Q9_N);
  assert(GetInvDirection(Q9_NE) == Q9_SW);
  assert(GetInvDirection(Q9_NW) == Q9_SE);
  assert(GetInvDirection(Q9_SW) == Q9_NE);
  assert(GetInvDirection(Q9_SE) == Q9_NW);
}

void LatticeBoltzmann::SetupDataStructures(const StdVector<double>& elements)
{
  rel.Clear();
  bb.Clear();
  inlet.Clear();
  outlet.Clear();

  StdVector<int> tmp(2);

  int n = 0;
  double porosity;

  for(int j = 0; j < m_sizeY; j++)
  {
    for(int i = 0; i < m_sizeX; i++)
    {
      porosity = elements[n];

      if (LbmNodeTypeIsFluid(porosity)) {
        tmp[0] = i;
        tmp[1] = j;
        rel.Push_back(tmp);
      } else if (LbmNodeTypeIsBB(porosity)) {
        tmp[0] = i;
        tmp[1] = j;
        bb.Push_back(tmp);
      } else if (LbmNodeTypeIsInlet(porosity)) {
        tmp[0] = i;
        tmp[1] = j;
        inlet.Push_back(tmp);
      } else if (LbmNodeTypeIsOutlet(porosity)) {
        tmp[0] = i;
        tmp[1] = j;
        outlet.Push_back(tmp);
      }
      ++n;
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

inline void LatticeBoltzmann::InitPropMap(StdVector<PropTransform>& map)
{
  for (int dir = 0; dir < N_DIR; dir++) {
    map[dir] = PropTransform();
  }
}

void LatticeBoltzmann::create_output(const char * file, int cur)
{
  // for debug purposes
  fstream f;
  f.precision(16);
  f.open(file, ios::out);

  for(int i = 0; i < m_nNodes; i++) {
    for(int j = 0; j < N_DIR; j++) {
      f << PDF_IDX(cur, i, j) << " ";
    }

    f << std::endl;
  }

  f.close();
}


void LatticeBoltzmann::prop_step()
{
  // perform a propagation step
  int x, y;

  int lx = m_sizeX;
  int ly = m_sizeY;

  int cur = m_cur;


  for(y = ly - 1; y >= 0; y--) {
    for(x = lx - 1; x > 0; x--) {
      PDF(cur, x, y, Q9_E) = PDF(cur, x - 1, y, Q9_E);
    }

    for(x = 0; x < lx - 1; x++) {
      PDF(cur, x, y, Q9_W) = PDF(cur, x + 1, y, Q9_W);
    }
  }

  for(y = ly - 1; y > 0; y--) {
    for(x = lx - 1; x > 0; x--) {
      PDF(cur, x, y, Q9_NE) = PDF(cur, x - 1, y - 1, Q9_NE);
    }

    for(x = 0; x < lx - 1; x++) {
      PDF(cur, x, y, Q9_NW) = PDF(cur, x + 1, y - 1, Q9_NW);
    }
  }

  for(y = 0; y < ly - 1; y++) {
    for(x = lx - 1; x > 0; x--) {
      PDF(cur, x, y, Q9_SE) = PDF(cur, x - 1, y + 1, Q9_SE);
    }

    for(x = 0; x < lx - 1; x++) {
      PDF(cur, x, y, Q9_SW) = PDF(cur, x + 1, y + 1, Q9_SW);
    }
  }

  for(x = 0; x < lx; x++) {
    for(y = ly - 1; y > 0; y--) {
      PDF(cur, x, y, Q9_N) = PDF(cur, x, y - 1, Q9_N);
    }

    for(y = 0; y < ly - 1; y++) {
      PDF(cur, x, y, Q9_S) = PDF(cur, x, y + 1, Q9_S);
    }
  }

  return;
}


void LatticeBoltzmann::prop_coll_step(int m_cur, int m_next, double omega)
{
  // perform a propagation step
  int x, y, index;

//  Direction dir;(Direction)

  const int lx = m_sizeX;
  const int ly = m_sizeY;

  double tmp_ux, tmp_uy, tmp_us, scale, sum;
  double tmp;

  // ---------------------------------------------------------------------
  // Propagation for nodes in the corners....
  // ---------------------------------------------------------------------

  const Boundary corners[] = {NE, NW, SW, SE};
  const Boundary vEdges[] = {W, E};
  const Boundary hEdges[] = {S, N};
  StdVector<PropTransform> map1, map2;

  // Propagation in the corners
  for (int i = 0; i < 4; i++) {
    switch (corners[i])
    { case NE: x = lx - 1; y = ly - 1; break;
      case NW: x =  0; y = ly -1; break;
      case SW: x = 0; y = 0; break;
      case SE: x = lx -1; y = 0; break;
      default:
        EXCEPTION("Trying to propagate into a corner that does not exist!");
    }

    map1 = prop_maps[corners[i]];
    assert(map1.GetSize() == N_DIR);
    for (int j = 0; j < N_DIR; j++) {
//      dir = allDirections[j];
      PDF(m_next, x, y, j)  = PDF(m_cur, x+map1[j].off_x,     y+map1[j].off_y,     j);
    }
  }

  // Propagation along the boundaries of the X dimension.
  map1 = prop_maps[hEdges[0]];
  map2 = prop_maps[hEdges[1]];
  for (int x = 1; x < lx - 1; ++x) {
    for (int i = 0; i < N_DIR; i++) {
//      dir = allDirections[i];
      y = 0;
      PDF(m_next, x, y, i)  = PDF(m_cur, x+map1[i].off_x,     y+map1[i].off_y,     i);
      y = ly -1;
      PDF(m_next, x, y, i)  = PDF(m_cur, x+map2[i].off_x,     y+map2[i].off_y,     i);
    }
  }

  // Propagation along the boundaries of the Y dimension.
  map1 = prop_maps[vEdges[0]];
  map2 = prop_maps[vEdges[1]];
  for (int y = 1; y < ly - 1; ++y) {
    for (int i = 0; i < N_DIR; i++) {
//      dir = allDirections[i];
      x = 0;
      PDF(m_next, x, y, i)  = PDF(m_cur, x+map1[i].off_x,     y+map1[i].off_y,     i);
      x = lx - 1;
      PDF(m_next, x, y, i)  = PDF(m_cur, x+map2[i].off_x,     y+map2[i].off_y,     i);
    }
  }

  double * scales  = Scales.GetPointer();


  #ifdef _OPENMP
    #pragma omp parallel for collapse(2) \
      default(none), \
      private(index, pdf_0, pdf_w, pdf_s, pdf_e, pdf_n, pdf_sw, pdf_se, pdf_ne, pdf_nw), \
      private(tmp_ux, tmp_uy, tmp_us, scale, sum, tmp, x, y), \
      shared(m_cur, m_next, omega, scales)
  #endif

  StdVector<PropTransform> map = prop_maps[0];
  StdVector<double> pdfs;
  pdfs.Resize(N_DIR);

  for (y = 1; y < ly - 1; ++y) {
    for (x = 1; x < lx - 1; ++x) {

      index = y * lx + x;

      sum = 0;
      tmp_ux = 0;
      tmp_uy = 0;
      tmp_us = 0;

//      for (int i = 0; i < N_DIR; i++) {
////        dir = allDirections[i];
//        //store current pdf values in array for better accessing
//        pdfs[i] = PDF(m_cur, x+map[i].off_x,     y+map[i].off_y,     i);
//        // add up all distribution functions of one element
//        sum += pdfs[i];
//        if (i != Q9_0) {
//          if (i != Q9_N && i != Q9_S) {
//            tmp_ux += map[GetInvDirection((Direction)i)].off_x*pdfs[i];
//          }
//          if (i != Q9_E && i != Q9_W)
//            tmp_uy += map[GetInvDirection((Direction)i)].off_y*pdfs[i];
//        }
//      }
//      pdfs[Q9_0]  = PDF(m_cur, x,     y,     Q9_0);
//      pdfs[Q9_E]  = PDF(m_cur, x - 1, y,     Q9_E);
//      pdfs[Q9_N]  = PDF(m_cur, x,     y - 1, Q9_N);
//      pdfs[Q9_W]  = PDF(m_cur, x + 1, y,     Q9_W);
//      pdfs[Q9_S]  = PDF(m_cur, x,     y + 1, Q9_S);
//      pdfs[Q9_NE] = PDF(m_cur, x - 1, y - 1, Q9_NE);
//      pdfs[Q9_NW] = PDF(m_cur, x + 1, y - 1, Q9_NW);
//      pdfs[Q9_SW] = PDF(m_cur, x + 1, y + 1, Q9_SW);
//      pdfs[Q9_SE] = PDF(m_cur, x - 1, y + 1, Q9_SE);

      sum = 0;
      for (int i = 0; i < N_DIR; i++) {
        pdfs[i] = PDF(m_cur, x+map[i].off_x,     y+map[i].off_y,     i);
        sum += pdfs[i];
        tmp_ux += (-1) * map[i].off_x*pdfs[i];
        tmp_uy += (-1) * map[i].off_y*pdfs[i];
      }

//      if (sum > 2) {
//        std::cout << "sum = " << sum << " x = " << x << " y = " << y <<  std::endl;
//      }
//      if (tmp_ux < 1) {
//        std::cout << "tmp_ux = " << tmp_ux << " x = " << x << " y = " << y << std::endl;
//      }
//      if (tmp_uy < 1) {
//        std::cout << "tmp_uy = " << tmp_uy << " x = " << x << " y = " << y << std::endl;
//      }

//      sum = pdfs[Q9_0] + pdfs[Q9_E] + pdfs[Q9_N] + pdfs[Q9_W] + pdfs[Q9_S] + pdfs[Q9_NE] + pdfs[Q9_NW] + pdfs[Q9_SW] + pdfs[Q9_SE];

      // macroscopic scaling by design variable
      scale = scales[index];

      tmp_ux = scale * tmp_ux / sum;
      tmp_uy = scale * tmp_uy / sum;
//      tmp_ux = scale * (pdfs[Q9_E] + pdfs[Q9_NE] + pdfs[Q9_SE] - pdfs[Q9_W] - pdfs[Q9_NW] - pdfs[Q9_SW]) / sum;
//      tmp_uy = scale * (pdfs[Q9_N] + pdfs[Q9_NE] + pdfs[Q9_NW] - pdfs[Q9_S] - pdfs[Q9_SW] - pdfs[Q9_SE]) / sum;
      tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);

      tmp_ux = 3.0 * tmp_ux;
      tmp_uy = 3.0 * tmp_uy;

      // propagate and collide for inner elements
      for (int i = 0; i < N_DIR; i++) {
//        dir = allDirections[i];
        tmp = map[GetInvDirection((Direction)i)].off_x * tmp_ux + map[GetInvDirection((Direction)i)].off_y * tmp_uy;
        PDF(m_next, x, y, i) = pdfs[i] + omega * ((sum * weights[i]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdfs[i]);
      }
    }
  }
  return;
}

void LatticeBoltzmann::prop_coll_velinlet(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY)
{

  int x, y;
  double tmp_ux, tmp_uy, tmp_us, sum;
  double tmp;

  StdVector<double> pdfs;
  pdfs.Resize(N_DIR);
  StdVector<PropTransform> map = prop_maps[0];

  for(unsigned int i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];

    for (int j = 0; j < N_DIR; j++) {
      pdfs[j] = PDF(cur, x, y, j);
    }

    sum = 0;
    for (int j = 0; j < N_DIR; j++)
      sum += pdfs[j];

    tmp_ux = UX;
    tmp_uy = UY;
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;

    LOG_DBG3(lattice) << "pcv: i=" << i << " tux=" << tmp_ux << " tuy=" << tmp_uy << std::endl;


    for (int i = 0; i < N_DIR; i++) {
      tmp = map[GetInvDirection((Direction)i)].off_x * tmp_ux + map[GetInvDirection((Direction)i)].off_y * tmp_uy;
      PDF(cur, x, y, i)  = sum * weights[i]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


//
// Performs a bounce back step.
//
void LatticeBoltzmann::prop_coll_bounce_back(int cur, StdVector<StdVector<int> >& bb)
{
  int x;
  int y;

  StdVector<double> pdfs;
  pdfs.Resize(N_DIR);

  for(unsigned int i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];

    for (int j = 0; j < N_DIR; j++) {
      pdfs[j] = PDF(cur, x, y, j);
    }

    for (int j = 0; j < N_DIR; j++) {
      PDF(cur, x, y, GetInvDirection((Direction)j)) = pdfs[j] ;
    }
  }

  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::prop_coll_densoutlet(int cur, StdVector<StdVector<int> >& outlet)
{

  double tmp_ux, tmp_uy, tmp_us, sum, tmp;

  StdVector<PropTransform> map = prop_maps[0];
  StdVector<double> pdfs;
  pdfs.Resize(N_DIR);

  int x;
  int y;

  for(unsigned int i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];

    sum = 0;
    tmp_ux = 0;
    tmp_uy = 0;

    for (int i = 0; i < N_DIR; i++) {
      //store current pdf values in array for better accessing
      pdfs[i] = PDF(cur, x, y, i);
      // add up all distribution functions of one element
      sum += pdfs[i];
      if (i != Q9_0) {
        if (i != Q9_N && i != Q9_S) {
          tmp_ux += map[GetInvDirection((Direction)i)].off_x*pdfs[i];
        }
        if (i != Q9_E && i != Q9_W)
          tmp_uy += map[GetInvDirection((Direction)i)].off_y*pdfs[i];
      }
    }

    tmp_ux = tmp_ux / sum;
    tmp_uy = tmp_uy / sum;

    sum = 1.0; // the enforced density
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;


    for (int j = 0; j < N_DIR; j++) {
      tmp = map[GetInvDirection((Direction)j)].off_x * tmp_ux + map[GetInvDirection((Direction)j)].off_y * tmp_uy;
      PDF(cur, x, y, j) = sum * weights[j] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
    }
  }

  return;
}


} // end namespace
