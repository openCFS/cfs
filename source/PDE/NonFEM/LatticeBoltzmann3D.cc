#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "PDE/NonFEM/LatticeBoltzmann3D.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{

using std::fstream;
using std::ios;

DECLARE_LOG(lattice)
DEFINE_LOG(lattice, "lattice")

// instantiation of the static elements
Enum<LatticeBoltzmann3D::Boundary>        LatticeBoltzmann3D::boundaries;
Enum<LatticeBoltzmann3D::Direction>         LatticeBoltzmann3D::directions;
//Enum<LatticeBoltzmann3D::Corner>            LatticeBoltzmann3D::corners;
//Enum<LatticeBoltzmann3D::Edge>              LatticeBoltzmann3D::edges;

 LatticeBoltzmann3D::LatticeBoltzmann3D(int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, double omega, int maxIterations, double maxTolerance, bool plot)
 {
   m_sizeX = sizeX;
   m_sizeY = sizeY;
   m_sizeZ = sizeZ;

   m_ux = ux;
   m_uy = uy;
   m_uz = uz;
   m_omega = omega;
   m_maxIter = maxIterations;
   m_maxTol = maxTolerance;

   m_nNodes = m_sizeX * m_sizeY * m_sizeZ;

   m_plot = plot;

   //matrix of the probability distributions
   LOG_DBG(lattice) << "Allocating arrays for " << m_nNodes << " PDFs (" << (sizeof(double) * m_nNodes * N_DIR * 2.0 / 1024.0 / 1024.0) << " MiB)";

   // N_DIR: number of discrete directions in this model, e.g. 19 for D3Q19
   m_pdfs.Resize(2);
   m_pdfs[0].Resize(m_nNodes * N_DIR);
   m_pdfs[1].Resize(m_nNodes * N_DIR);

   // 8 corners + 12 edges + 1 map for interior nodes
   prop_maps.Resize(N_DIR + 2);
   for (int i = 0; i < N_DIR +2; i++ ){
     prop_maps[i].Resize(N_DIR + 2);
   }

   SetEnums();
   InitWeights();
   SetInvDirections();
   TestDirectionIndex();
   TestInvDirections();
   TestGetIndexbound();
   SetupTransformation();

   // for debugging
   WritePropMap();

   Scales.Resize(m_nNodes);
   rel.Resize(m_nNodes);
   bb.Resize(2 * m_sizeX + 2 * m_sizeY + 2 * m_sizeZ);

   m_cur  = 0;
   m_next = 1;
 }

 LatticeBoltzmann3D::~LatticeBoltzmann3D()
 {
 }

StdVector<double>* LatticeBoltzmann3D::Iterate(const StdVector<double>& elements, PtrParamNode in)
{
  int index;
  int it = 0;

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

  while(it < m_maxIter && !steady_state && R <= 1000)
  {
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
      for (int k = 0;  k < m_sizeZ; k++) {
        for(int j = 0; j < m_sizeY; j++) {
          for(int i = 0; i < m_sizeX; i++) {
            index = k * m_sizeX * m_sizeY + j * m_sizeX + i;

            for(int k = 0; k < N_DIR; k++)
            {
              res = PDF_IDX(m_next, index, k) - PDF_IDX(m_cur, index, k);
              R += res * res;
            }
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
  double performance = (m_sizeX - 1) * (m_sizeY - 1) * (m_sizeZ - 1) * it / wt / 1e6;
  in->Get("MFLUP_s")->SetValue(performance);
  in->Get("walltime")->SetValue(wt);

  if(R >= 1000) 
    EXCEPTION("In LBM iteration " << it << " residuum " << R << " too large ... abort");

  if(!steady_state)
    EXCEPTION("internal LBM simulation could not converge: iterations: " << it << " residuum: " << R);

  return &(m_pdfs[m_cur]);
}


void LatticeBoltzmann3D::InitializePdfs()
{
  int i;

  #ifdef _OPENMP
    // Perform NUMA-correct placement.
    #pragma omp parallel for default(none), private(i), collapse(2)
  #endif
  double weight;
  for (int y = 0; y < m_sizeY; ++y) {
    for (int x = 0; x < m_sizeX; ++x) {

      i = y * m_sizeX + x;

      // weight for center
      weight = weights[0];
      assert(weight == 1./3.);
      PDF_IDX(0, i, 0) = weight;
      PDF_IDX(1, i, 0) = weight;

      // weight for direct neighbors
      weight = weights[1];
      assert(weight == 1./18.);
      for(int j = 1; j < 7; j++) {
        PDF_IDX(0, i, j) = weight;
        PDF_IDX(1, i, j) = weight;
      }

      // weight for all other neighbors
      weight = weights[N_DIR - 1];
      assert(weight == 1./36.);
      for(int j = 7; j < N_DIR; j++) {
        PDF_IDX(0, i, j) = weight;
        PDF_IDX(1, i, j) = weight;
      }
    }
  }
}

//returns associated integer value of velocity direction two given principal directions
inline int LatticeBoltzmann3D::GetIndexDir(Direction dir1)
{
  assert(directions.IsValid(dir1));
  return dir1;
}


// first direction must be: S, N, T or B
// second direction can only be: E, W, N, S
inline int LatticeBoltzmann3D::GetIndexDir(Direction dir1, Direction dir2)
{
  assert(directions.IsValid(dir1));
  assert(dir2 != 0);
  assert(directions.IsValid(dir2));
  assert(dir1 != dir2);
  assert((dir1 == Q19_S) || (dir1 == Q19_N) || (dir1 == Q19_T) || (dir1 == Q19_B));
  assert((dir2 == Q19_E) || (dir2 == Q19_W) || (dir2 == Q19_N) || (dir2 == Q19_S));
  // this formula is dependent on indexing of directions
  if (dir1 == Q19_N && dir2 == Q19_E)
    return Q19_NE;
  if (dir1 ==  Q19_N && dir2 ==  Q19_W)
    return  Q19_NW;
  if (dir1 == Q19_S  && dir2 == Q19_W)
      return  Q19_SW;
  if (dir1 ==  Q19_S && dir2 == Q19_E)
      return Q19_SE ;
  if (dir1 == Q19_T  && dir2 == Q19_N)
      return  Q19_TN;
  if (dir1 == Q19_B  && dir2 ==  Q19_S)
      return  Q19_BS;
  if (dir1 == Q19_T  && dir2 == Q19_S)
      return Q19_TS ;
  if (dir1 == Q19_B  && dir2 ==  Q19_N)
      return  Q19_BN;
  if (dir1 ==  Q19_T && dir2 == Q19_E)
    return  Q19_TE;
  if (dir1 == Q19_B  && dir2 == Q19_W)
    return  Q19_BW;
  if (dir1 == Q19_T  && dir2 ==  Q19_W)
    return  Q19_TW;
  if (dir1 == Q19_B  && dir2 ==  Q19_E)
    return  Q19_BE;

  WARN("I should not be here!");
  return -1;
}

//returns associated integer value of boundary enum(edges)
// dir1 can only be top, bottom, north or south
// dir2 can only be north, south, west or east
inline int LatticeBoltzmann3D::GetIndexBound(Direction dir1, Direction dir2)
{
  assert(directions.IsValid(dir1));
  assert(dir2 != 0);
  assert(directions.IsValid(dir2));
  if (dir1 == Q19_T && dir2 == Q19_N)
    return TN;
  if (dir1 == Q19_T && dir2 == Q19_W)
    return TW;
  if (dir1 == Q19_T && dir2 == Q19_S)
    return TS;
  if (dir1 == Q19_T && dir2 == Q19_E)
    return TE;
  if (dir1 == Q19_S && dir2 == Q19_W)
      return SW;
  if (dir1 == Q19_S && dir2 == Q19_E)
    return SE;
  if (dir1 == Q19_N && dir2 == Q19_W)
      return NW;
  if (dir1 == Q19_N && dir2 == Q19_E)
    return NE;
  if (dir1 == Q19_B && dir2 == Q19_N)
    return BN;
  if (dir1 == Q19_B && dir2 == Q19_W)
    return BW;
  if (dir1 == Q19_B && dir2 == Q19_S)
    return BS;
  if (dir1 == Q19_B && dir2 == Q19_E)
    return BE;
  return -1;
}

//returns associated integer value of boundary enum(corners)
// dir1 can only be top or bottom
// dir2 can only be north or south
// dir3 can only be east or west
inline int LatticeBoltzmann3D::GetIndexBound(Direction dir1, Direction dir2, Direction dir3)
{
  assert(directions.IsValid(dir1));
  assert(dir2 != 0);
  assert(directions.IsValid(dir2));
  assert(dir3 != 0);
  assert(directions.IsValid(dir3));
  assert(dir1 == Q19_T || dir1 == Q19_B);
  assert(dir2 == Q19_N || dir2 == Q19_S);
  assert(dir3 == Q19_E || dir3 == Q19_W);

  if (dir1 == Q19_B && dir2 == Q19_N && dir3 == Q19_W )
    return BNW;
  if (dir1 == Q19_B && dir2 == Q19_S && dir3 == Q19_W)
      return BSW;
  if (dir1 == Q19_B && dir2 == Q19_N && dir3 == Q19_E)
      return BNE;
  if (dir1 == Q19_B && dir2 == Q19_S && dir3 == Q19_E)
      return BSE;
  if (dir1 == Q19_T && dir2 == Q19_N && dir3 == Q19_W )
    return TNW;
  if (dir1 == Q19_T && dir2 == Q19_S && dir3 == Q19_W)
    return TSW;
  if (dir1 == Q19_T && dir2 == Q19_N && dir3 == Q19_E)
    return TNE;
  if (dir1 == Q19_T && dir2 == Q19_S && dir3 == Q19_E)
    return TSE;
  return -1;
}

// validates GetIndexbound function
void LatticeBoltzmann3D::TestGetIndexbound()
{
  assert(GetIndexBound(Q19_B, Q19_N, Q19_W) == 1);
  assert(GetIndexBound(Q19_B,Q19_S,Q19_W) == 2);
  assert(GetIndexBound(Q19_B,Q19_S,Q19_E) == 3);
  assert(GetIndexBound(Q19_B,Q19_N,Q19_E) == 4);
  assert(GetIndexBound(Q19_T,Q19_N,Q19_W) == 5);
  assert(GetIndexBound(Q19_T,Q19_S,Q19_W) == 6);
  assert(GetIndexBound(Q19_T,Q19_S,Q19_E) == 7);
  assert(GetIndexBound(Q19_T,Q19_N,Q19_E) == 8);
  assert(GetIndexBound(Q19_T,Q19_N) == 9);
  assert(GetIndexBound(Q19_T,Q19_W) == 10);
  assert(GetIndexBound(Q19_T,Q19_S) == 11);
  assert(GetIndexBound(Q19_T,Q19_E) == 12);
  assert(GetIndexBound(Q19_N,Q19_W) == 13);
  assert(GetIndexBound(Q19_S,Q19_W) == 14);
  assert(GetIndexBound(Q19_S,Q19_E) == 15);
  assert(GetIndexBound(Q19_N,Q19_E) == 16);
  assert(GetIndexBound(Q19_B,Q19_N) == 17);
  assert(GetIndexBound(Q19_B,Q19_W) == 18);
  assert(GetIndexBound(Q19_B,Q19_S) == 19);
  assert(GetIndexBound(Q19_B,Q19_E) == 20);
}


// validates GetIndexDir function
void LatticeBoltzmann3D::TestDirectionIndex()
{
  assert(GetIndexDir(Q19_0) == 0);
  assert(GetIndexDir(Q19_E) == 1);
  assert(GetIndexDir(Q19_W) == 2);
  assert(GetIndexDir(Q19_N) == 3);
  assert(GetIndexDir(Q19_S) == 4);
  assert(GetIndexDir(Q19_T) == 5);
  assert(GetIndexDir(Q19_B) == 6);

  assert(GetIndexDir(Q19_N, Q19_E) == 7);
  assert(GetIndexDir(Q19_S, Q19_W) == 8);
  assert(GetIndexDir(Q19_N, Q19_W) == 9);
  assert(GetIndexDir(Q19_S, Q19_E) == 10);
  assert(GetIndexDir(Q19_T, Q19_N) == 11);
  assert(GetIndexDir(Q19_B, Q19_S) == 12);
  assert(GetIndexDir(Q19_T, Q19_S) == 13);
  assert(GetIndexDir(Q19_B, Q19_N) == 14);
  assert(GetIndexDir(Q19_T, Q19_E) == 15);
  assert(GetIndexDir(Q19_B, Q19_W) == 16);
  assert(GetIndexDir(Q19_T, Q19_W) == 17);
  assert(GetIndexDir(Q19_B, Q19_E) == 18);
}

void LatticeBoltzmann3D::SetupTransformation()
{
  // in 3D: 20 maps (8 corners + 12 edges)
  StdVector<PropTransform> map;
  map.Resize(N_DIR);
//  const Direction xDirs[] = {Q19_E,Q19_W};
//  const Direction yDirs[] = {Q19_N,Q19_S};
  Direction d1, d2, d3;
  const Direction dirs1[] = {Q19_T, Q19_B};
  const Direction dirs2[] = {Q19_N, Q19_S};
  const Direction dirs3[] = {Q19_E, Q19_W};

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
        map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q19_B ? 1 : -1));
        // d2 can on only be north or south
        map[GetIndexDir(d2)] = PropTransform(0, (d2 == Q19_S ? 1 : -1) , 0);
        // d3 can on only be east or west
        map[GetIndexDir(d3)] = PropTransform((d3 == Q19_W ? 1 : -1), 0, 0);

        map[GetIndexDir(d1,d2)] = PropTransform(0, (d2 == Q19_S ? 1 : -1) , (d1 == Q19_B ? 1 : -1));
        map[GetIndexDir(d1,d3)] = PropTransform((d3 == Q19_W ? 1 : -1), 0, (d1 == Q19_B ? 1 : -1));
        map[GetIndexDir(d2,d3)] = PropTransform((d3 == Q19_W ? 1 : -1), (d2 == Q19_S ? 1 : -1), 0);

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
      map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q19_B ? 1 : -1));
      // d2 can on only be north or south
      map[GetIndexDir(d2)] = PropTransform(0, (d2 == Q19_S ? 1 : -1) , 0);
      map[GetIndexDir(d1,d2)] = PropTransform(0, (d2 == Q19_S ? 1 : -1) , (d1 == Q19_B ? 1 : -1));
      for (int k = 0; k < 2; k++) { //E, W
        d3 = dirs3[k];
        map[GetIndexDir(d1,d3)] = PropTransform((d3 == Q19_W ? 1 : -1), 0, (d1 == Q19_B ? 1 : -1));
        map[GetIndexDir(d2,d3)] = PropTransform((d3 == Q19_W ? 1 : -1), (d2 == Q19_S ? 1 : -1), 0);
        // d3 can on only be east or west
        map[GetIndexDir(d3)] = PropTransform((d3 == Q19_W ? 1 : -1), 0, 0);
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
      map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q19_B ? 1 : -1));
      // d2 can on only be east or west
      map[GetIndexDir(d2)] = PropTransform((d2 == Q19_W ? 1 : -1), 0, 0);
      map[GetIndexDir(d1,d2)] = PropTransform((d2 == Q19_W ? 1 : -1), 0, (d1 == Q19_B ? 1 : -1));
      for (int k = 0; k < 2; k++) { //N, S
        d3 = dirs2[k];
        // d3 can on only be north or south
        map[GetIndexDir(d3)] = PropTransform(0, (d3 == Q19_S ? 1 : -1) , 0);
        map[GetIndexDir(d1,d3)] = PropTransform(0,(d3 == Q19_S ? 1 : -1),(d1 == Q19_B ? 1 : -1));
        map[GetIndexDir(d3,d2)] = PropTransform((d2 == Q19_W ? 1 : -1),(d3 == Q19_S ? 1 : -1),0);
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
      map[GetIndexDir(d1)] = PropTransform(0, (d1 == Q19_S ? 1 : -1) , 0);
      // d2 can on only be east or west
      map[GetIndexDir(d2)] = PropTransform((d2 == Q19_W ? 1 : -1), 0, 0);

      map[GetIndexDir(d1,d2)] = PropTransform((d2 == Q19_W ? 1 : -1),(d1 == Q19_S ? 1 : -1),0);

      for (int k = 0; k < 2; k++) { // T, B
        d3 = dirs1[k];
        // d3 can only be top or bottom
        map[GetIndexDir(d3)] = PropTransform(0, 0, (d3 == Q19_B ? 1 : -1));
        map[GetIndexDir(d3,d1)] = PropTransform(0,(d1 == Q19_S ? 1 : -1),(d3 == Q19_B ? 1 : -1));
        map[GetIndexDir(d3,d2)] = PropTransform((d2 == Q19_W ? 1 : -1),0,(d3 == Q19_B ? 1 : -1));
      }
      prop_maps[GetIndexBound(d1,d2)] = map;
    }
  }

  // propagation for inner nodes
  // fill map to default offsets
  InitPropMap(map);
  for (int i = 0; i < 2; i++) { // T, B
    d1 = dirs1[i];
    // d1 can only be top or bottom
    map[GetIndexDir(d1)] = PropTransform(0, 0, (d1 == Q19_B ? 1 : -1));
    for (int j = 0; j < 2; j++) { // N, S
      d2 = dirs2[j];
      // d2 can on only be north or south
      map[GetIndexDir(d2)] = PropTransform(0, (d2 == Q19_S ? 1 : -1) , 0);

      map[GetIndexDir(d1,d2)] = PropTransform(0,(d2 == Q19_S ? 1 : -1),(d1 == Q19_B ? 1 : -1));
      for (int k = 0; k < 2; k++) { // E, W
        d3 = dirs3[k];

        // d3 can on only be east or west
        map[GetIndexDir(d3)] = PropTransform((d3 == Q19_W ? 1 : -1), 0, 0);

        map[GetIndexDir(d1,d3)] = PropTransform((d3 == Q19_W ? 1 : -1),0,(d1 == Q19_B ? 1 : -1));
        map[GetIndexDir(d2,d3)] = PropTransform((d3 == Q19_W ? 1 : -1),(d2 == Q19_S ? 1 : -1),0);
      }
    }
  }
  prop_maps[C] = map;
}

// for debugging purposes
void LatticeBoltzmann3D::WritePropMap()
{
  // write out propagation maps for corners
  const Boundary corners[] = {BNW, BSW, BNE, BSE, TNW, TSW, TNE, TSE};
  const Direction order[] = {Q19_0, Q19_E, Q19_N, Q19_W, Q19_S, Q19_T, Q19_B, Q19_NE, Q19_NW, Q19_SW, Q19_SE, Q19_TN, Q19_TS, Q19_BS, Q19_BN, Q19_TE, Q19_BW, Q19_TW, Q19_BE};

  fstream f;
  f.precision(16);
  f.open("propagation_corners.txt", ios::out);
  for (int i = 0; i < 8; i++) {
    StdVector<PropTransform>& mymap  = prop_maps[corners[i]];
    f << boundaries.ToString(corners[i]) << std::endl;
    f << "---------------------------------------" << std::endl;
    for (int j = 0; j < N_DIR; j++) {
      f << "PDF(cur, x, y, z, " << directions.ToString(order[j]) << ") = PDF(cur, x + " << mymap[order[j]].off_x << ", y +" << mymap[order[j]].off_y << ", z + "
          << mymap[order[j]].off_z << ", " << directions.ToString(order[j]) << ") \n";
    }
    f << std::endl;
  }
  f.close();

  // write out propagation maps for edges
  const Boundary edges[] = {TN, TW, TS, TE, NW, SW, SE, NE, BN, BW, BS, BE};
  f.open("propagation_edges.txt", ios::out);
  for (int i = 0; i < 12; i++) {
    StdVector<PropTransform>& mymap  = prop_maps[edges[i]];
    f << boundaries.ToString(edges[i]) << std::endl;
    for (int j = 0; j < N_DIR; j++) {
      f << "PDF(cur, x, y, z, " << directions.ToString(order[j]) << ") = PDF(cur, x + " << mymap[order[j]].off_x << ", y +" << mymap[order[j]].off_y << ", z + "
          << mymap[order[j]].off_z << ", " << directions.ToString(order[j]) << ") \n";
    }
    f << std::endl;
  }
  // write out propagation maps for interior nodes
  StdVector<PropTransform>& mymap  = prop_maps[C];
  f << boundaries.ToString(C) << std::endl;
  for (int j = 0; j < N_DIR; j++) {
    f << "PDF(cur, x, y, z, " << directions.ToString(order[j]) << ") = PDF(cur, x + " << mymap[order[j]].off_x << ", y +" << mymap[order[j]].off_y << ", z + "
        << mymap[order[j]].off_z << ", " << directions.ToString(order[j]) << ") \n";
  }
  f << std::endl;
  f.close();
}

void LatticeBoltzmann3D::SetEnums()
{
  directions.SetName("Q19 directions");
  directions.Add(Q19_0,"C");
  directions.Add(Q19_E,"E");
  directions.Add(Q19_N,"N");
  directions.Add(Q19_W,"W");
  directions.Add(Q19_S,"S");
  directions.Add(Q19_T,"T");
  directions.Add(Q19_B,"B");
  directions.Add(Q19_NE,"NE");
  directions.Add(Q19_NW,"NW");
  directions.Add(Q19_SW,"SW");
  directions.Add(Q19_SE,"SE");
  directions.Add(Q19_TN,"TN");
  directions.Add(Q19_BS,"BS");
  directions.Add(Q19_TS,"TS");
  directions.Add(Q19_BN,"BN");
  directions.Add(Q19_TE,"TE");
  directions.Add(Q19_BW,"BW");
  directions.Add(Q19_TW,"TW");
  directions.Add(Q19_BE,"BE");

//  corners.SetName("Corners of a cube");
//  corners.Add(B_NW, "left bottom back");
//  corners.Add(B_SW, "left bottom front");
//  corners.Add(B_SE, "right bottom front");
//  corners.Add(B_NE, "right bottom back");
//  corners.Add(T_NW, "left top back");
//  corners.Add(T_SW, "left top front");
//  corners.Add(T_SE, "right top front");
//  corners.Add(T_NE, "right top back");
//
//  edges.SetName("Edges of a cube");
//  edges.Add(TN ,"top back");
//  edges.Add(TW ,"top left");
//  edges.Add(TS ,"top front");
//  edges.Add(TE ,"top right");
//  edges.Add(NW ,"left back");
//  edges.Add(SW ,"left front");
//  edges.Add(SE ,"right front");
//  edges.Add(NE ,"right back");
//  edges.Add(BN ,"bottom back");
//  edges.Add(BW ,"bottom left");
//  edges.Add(BS ,"bottom front");
//  edges.Add(BE ,"bottom right");

  boundaries.SetName("Corners and edges of a cube");
  boundaries.Add(C,"interior nodes");
  boundaries.Add(BNW, "left bottom back");
  boundaries.Add(BSW, "left bottom front");
  boundaries.Add(BSE, "right bottom front");
  boundaries.Add(BNE, "right bottom back");
  boundaries.Add(TNW, "left top back");
  boundaries.Add(TSW, "left top front");
  boundaries.Add(TSE, "right top front");
  boundaries.Add(TNE, "right top back");
  boundaries.Add(TN ,"top back");
  boundaries.Add(TW ,"top left");
  boundaries.Add(TS ,"top front");
  boundaries.Add(TE ,"top right");
  boundaries.Add(NW ,"left back");
  boundaries.Add(SW ,"left front");
  boundaries.Add(SE ,"right front");
  boundaries.Add(NE ,"right back");
  boundaries.Add(BN ,"bottom back");
  boundaries.Add(BW ,"bottom left");
  boundaries.Add(BS ,"bottom front");
  boundaries.Add(BE ,"bottom right");
}

void LatticeBoltzmann3D::InitWeights()
{
  weights.Resize(N_DIR);
  weights[Q19_0] = 1. / 3.;
  weights[Q19_E] = 1. / 18.;
  weights[Q19_N] = 1. / 18.;
  weights[Q19_W] = 1. / 18.;
  weights[Q19_S] = 1. / 18.;
  weights[Q19_T] = 1. / 18.;
  weights[Q19_B] = 1. / 18.;
  weights[Q19_NE] = 1. / 36.;
  weights[Q19_NW] = 1. / 36.;
  weights[Q19_SW] = 1. / 36.;
  weights[Q19_SE] = 1. / 36.;
  weights[Q19_TN] = 1. / 36.;
  weights[Q19_BS] = 1. / 36.;
  weights[Q19_TS] = 1. / 36.;
  weights[Q19_BN] = 1. / 36.;
  weights[Q19_TE] = 1. / 36.;
  weights[Q19_BW] = 1. / 36.;
  weights[Q19_TW] = 1. / 36.;
  weights[Q19_BE] = 1. / 36.;
}

void LatticeBoltzmann3D::SetInvDirections()
{
  directionsInv.Resize(N_DIR);
  directionsInv[Q19_0] = Q19_0;
  directionsInv[Q19_E] = Q19_W;
  directionsInv[Q19_N] = Q19_S;
  directionsInv[Q19_W] = Q19_E;
  directionsInv[Q19_S] = Q19_N;
  directionsInv[Q19_T] = Q19_B;
  directionsInv[Q19_B] = Q19_T;
  directionsInv[Q19_NE] = Q19_SW;
  directionsInv[Q19_NW] = Q19_SE;
  directionsInv[Q19_SW] = Q19_NE;
  directionsInv[Q19_SE] = Q19_NW;
  directionsInv[Q19_TN] = Q19_BS;
  directionsInv[Q19_TS] = Q19_BN;
  directionsInv[Q19_BS] = Q19_TN;
  directionsInv[Q19_BN] = Q19_TS;
  directionsInv[Q19_TE] = Q19_BW;
  directionsInv[Q19_BW] = Q19_TE;
  directionsInv[Q19_TW] = Q19_BE;
  directionsInv[Q19_BE] = Q19_TW;
}
// depends on numbering of directions
inline int LatticeBoltzmann3D::GetInvDirection(Direction dir)
{
//  assert(directions.IsValid(dir));
  return directionsInv[dir];
}

void LatticeBoltzmann3D::TestInvDirections()
{
  assert(GetInvDirection(Q19_0) == Q19_0);
  assert(GetInvDirection(Q19_E) == Q19_W);
  assert(GetInvDirection(Q19_N) == Q19_S);
  assert(GetInvDirection(Q19_W) == Q19_E);
  assert(GetInvDirection(Q19_S) == Q19_N);
  assert(GetInvDirection(Q19_T) == Q19_B);
  assert(GetInvDirection(Q19_B) == Q19_T);
  assert(GetInvDirection(Q19_NE) == Q19_SW);
  assert(GetInvDirection(Q19_NW) == Q19_SE);
  assert(GetInvDirection(Q19_SW) == Q19_NE);
  assert(GetInvDirection(Q19_SE) == Q19_NW);

  assert(GetInvDirection(Q19_TN) == Q19_BS);
  assert(GetInvDirection(Q19_TS) == Q19_BN);
  assert(GetInvDirection(Q19_BS) == Q19_TN);
  assert(GetInvDirection(Q19_BN) == Q19_TS);
  assert(GetInvDirection(Q19_TE) == Q19_BW);
  assert(GetInvDirection(Q19_BW) == Q19_TE);
  assert(GetInvDirection(Q19_TW) == Q19_BE);
  assert(GetInvDirection(Q19_BE) == Q19_TW);
}

void LatticeBoltzmann3D::SetupDataStructures(const StdVector<double>& elements)
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


std::string LatticeBoltzmann3D::ToString(const StdVector<StdVector<int> >& data)
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
inline void LatticeBoltzmann3D::InitPropMap(StdVector<PropTransform>& map)
{
  for (int dir = 0; dir < N_DIR; dir++) {
    map[dir] = PropTransform();
  }
}

void LatticeBoltzmann3D::create_output(const char * file, int cur)
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


void LatticeBoltzmann3D::prop_step()
{

  // TODO implement
  // perform a propagation step
//  int x, y, z;
//
//  int lx = m_sizeX;
//  int ly = m_sizeY;
//  int lz = m_sizeZ;
//
//  int cur = m_cur;
//
//
//  for(y = ly - 1; y >= 0; y--) {
//    for(x = lx - 1; x > 0; x--) {
//      PDF(cur, x, y, Q19_E) = PDF(cur, x - 1, y, Q19_E);
//    }
//
//    for(x = 0; x < lx - 1; x++) {
//      PDF(cur, x, y, Q19_W) = PDF(cur, x + 1, y, Q19_W);
//    }
//  }
//
//  for(y = ly - 1; y > 0; y--) {
//    for(x = lx - 1; x > 0; x--) {
//      PDF(cur, x, y, Q19_NE) = PDF(cur, x - 1, y - 1, Q19_NE);
//    }
//
//    for(x = 0; x < lx - 1; x++) {
//      PDF(cur, x, y, Q19_NW) = PDF(cur, x + 1, y - 1, Q19_NW);
//    }
//  }
//
//  for(y = 0; y < ly - 1; y++) {
//    for(x = lx - 1; x > 0; x--) {
//      PDF(cur, x, y, Q19_SE) = PDF(cur, x - 1, y + 1, Q19_SE);
//    }
//
//    for(x = 0; x < lx - 1; x++) {
//      PDF(cur, x, y, Q19_SW) = PDF(cur, x + 1, y + 1, Q19_SW);
//    }
//  }
//
//  for(x = 0; x < lx; x++) {
//    for(y = ly - 1; y > 0; y--) {
//      PDF(cur, x, y, Q19_N) = PDF(cur, x, y - 1, Q19_N);
//    }
//
//    for(y = 0; y < ly - 1; y++) {
//      PDF(cur, x, y, Q19_S) = PDF(cur, x, y + 1, Q19_S);
//    }
//  }

  return;
}


void LatticeBoltzmann3D::prop_coll_step(int m_cur, int m_next, double omega)
{
  // TODO implement
//  // perform a propagation step
//  int x, y, index;
//
//  const int lx = m_sizeX;
//  const int ly = m_sizeY;
//
//  double tmp_ux, tmp_uy, tmp_us, scale, sum;
//  double tmp;
//
//  // ---------------------------------------------------------------------
//  // Propagation for nodes in the corners....
//  // ---------------------------------------------------------------------
//
//  const Boundary corners[] = {NE, NW, SW, SE};
//  const Boundary vEdges[] = {W, E};
//  const Boundary hEdges[] = {S, N};
//  StdVector<PropTransform>* map1, *map2;
//
//  // Propagation in the corners
//  for (int i = 0; i < 4; i++) {
//    switch (corners[i])
//    { case NE: x = lx - 1; y = ly - 1; break;
//      case NW: x =  0; y = ly -1; break;
//      case SW: x = 0; y = 0; break;
//      case SE: x = lx -1; y = 0; break;
//      default:
//        EXCEPTION("Trying to propagate into a corner that does not exist!");
//    }
//
//    map1 = &prop_maps[corners[i]];
//    assert(map1->GetSize() == N_DIR);
//    for (int j = 0; j < N_DIR; j++) {
//      PDF(m_next, x, y, j)  = PDF(m_cur, x+ (*map1)[j].off_x,     y+(*map1)[j].off_y,     j);
//    }
//  }
//
//  // Propagation along the boundaries of the X dimension.
//  map1 = &prop_maps[hEdges[0]];
//  map2 = &prop_maps[hEdges[1]];
//  for (int x = 1; x < lx - 1; ++x) {
//    for (int i = 0; i < N_DIR; i++) {
//      y = 0;
//      PDF(m_next, x, y, i)  = PDF(m_cur, x+(*map1)[i].off_x,     y+(*map1)[i].off_y,     i);
//      y = ly -1;
//      PDF(m_next, x, y, i)  = PDF(m_cur, x+(*map2)[i].off_x,     y+(*map2)[i].off_y,     i);
//    }
//  }
//
//  // Propagation along the boundaries of the Y dimension.
//  map1 = &prop_maps[vEdges[0]];
//  map2 = &prop_maps[vEdges[1]];
//  for (int y = 1; y < ly - 1; ++y) {
//    for (int i = 0; i < N_DIR; i++) {
//      x = 0;
//      PDF(m_next, x, y, i)  = PDF(m_cur, x+(*map1)[i].off_x,     y+(*map1)[i].off_y,     i);
//      x = lx - 1;
//      PDF(m_next, x, y, i)  = PDF(m_cur, x+(*map2)[i].off_x,     y+(*map2)[i].off_y,     i);
//    }
//  }
//
//  double * scales  = Scales.GetPointer();
//
//
//  #ifdef _OPENMP
//    #pragma omp parallel for collapse(2) \
//      default(none), \
//      private(index, pdf_0, pdf_w, pdf_s, pdf_e, pdf_n, pdf_sw, pdf_se, pdf_ne, pdf_nw), \
//      private(tmp_ux, tmp_uy, tmp_us, scale, sum, tmp, x, y), \
//      shared(m_cur, m_next, omega, scales)
//  #endif
//
//  const StdVector<PropTransform>& map = prop_maps[0];
//  StdVector<double> pdfs;
//  pdfs.Resize(N_DIR);
//
//  for (y = 1; y < ly - 1; ++y) {
//    for (x = 1; x < lx - 1; ++x) {
//
//      index = y * lx + x;
//
//      sum = 0;
//      tmp_ux = 0;
//      tmp_uy = 0;
//      tmp_us = 0;
//
//      for (int i = 0; i < N_DIR; i++) {
//        //store current pdf values in array for better accessing
//        pdfs[i] = PDF(m_cur, x+map[i].off_x,     y+map[i].off_y,     i);
//        // add up all distribution functions of one element
//        sum += pdfs[i];
//        if (i != Q19_0) {
//          if (i != Q19_N && i != Q19_S) {
//            tmp_ux += map[GetInvDirection((Direction)i)].off_x*pdfs[i];
//          }
//          if (i != Q19_E && i != Q19_W)
//            tmp_uy += map[GetInvDirection((Direction)i)].off_y*pdfs[i];
//        }
//      }
//
//      // macroscopic scaling by design variable
//      scale = scales[index];
//
//      tmp_ux = scale * tmp_ux / sum;
//      tmp_uy = scale * tmp_uy / sum;
//      tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
//
//      tmp_ux = 3.0 * tmp_ux;
//      tmp_uy = 3.0 * tmp_uy;
//
//      // propagate and collide for inner elements
//      for (int i = 0; i < N_DIR; i++) {
//        tmp = map[GetInvDirection((Direction)i)].off_x * tmp_ux + map[GetInvDirection((Direction)i)].off_y * tmp_uy;
//        PDF(m_next, x, y, i) = pdfs[i] + omega * ((sum * weights[i]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdfs[i]);
//      }
//    }
//  }
  return;
}

void LatticeBoltzmann3D::prop_coll_velinlet(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY)
{
  // TODO Implement!
//  int x, y;
//  double tmp_ux, tmp_uy, tmp_us, sum;
//  double tmp;
//
//  StdVector<double> pdfs;
//  pdfs.Resize(N_DIR);
//  const StdVector<PropTransform>& map = prop_maps[0];
//
//  for(unsigned int i = 0; i < inlet.GetSize(); i++) {
//    x = inlet[i][0];
//    y = inlet[i][1];
//
//    for (int j = 0; j < N_DIR; j++) {
//      pdfs[j] = PDF(cur, x, y, j);
//    }
//
//    sum = 0;
//    for (int j = 0; j < N_DIR; j++)
//      sum += pdfs[j];
//
//    tmp_ux = UX;
//    tmp_uy = UY;
//    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
//    tmp_ux = 3.0 * tmp_ux;
//    tmp_uy = 3.0 * tmp_uy;
//
//    LOG_DBG3(lattice) << "pcv: i=" << i << " tux=" << tmp_ux << " tuy=" << tmp_uy << std::endl;
//
//
//    for (int i = 0; i < N_DIR; i++) {
//      tmp = map[GetInvDirection((Direction)i)].off_x * tmp_ux + map[GetInvDirection((Direction)i)].off_y * tmp_uy;
//      PDF(cur, x, y, i)  = sum * weights[i]  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
//    }
//  }
//
  return;
}
//
//
////
//// Performs a bounce back step.
////
void LatticeBoltzmann3D::prop_coll_bounce_back(int cur, StdVector<StdVector<int> >& bb)
{
//  int x;
//  int y;
//
//  StdVector<double> pdfs;
//  pdfs.Resize(N_DIR);
//
//  for(unsigned int i = 0; i < bb.GetSize(); i++) {
//    x = bb[i][0];
//    y = bb[i][1];
//
//    for (int j = 0; j < N_DIR; j++) {
//      pdfs[j] = PDF(cur, x, y, j);
//    }
//
//    for (int j = 0; j < N_DIR; j++) {
//      PDF(cur, x, y, GetInvDirection((Direction)j)) = pdfs[j] ;
//    }
//  }

  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann3D::prop_coll_densoutlet(int cur, StdVector<StdVector<int> >& outlet)
{
  // TODO Implement!
//  double tmp_ux, tmp_uy, tmp_us, sum, tmp;
//
//  const StdVector<PropTransform>& map = prop_maps[0];
//  StdVector<double> pdfs;
//  pdfs.Resize(N_DIR);
//
//  int x;
//  int y;
//
//  for(unsigned int i = 0; i < outlet.GetSize(); i++) {
//    x = outlet[i][0];
//    y = outlet[i][1];
//
//    sum = 0;
//    tmp_ux = 0;
//    tmp_uy = 0;
//
//    for (int i = 0; i < N_DIR; i++) {
//      //store current pdf values in array for better accessing
//      pdfs[i] = PDF(cur, x, y, i);
//      // add up all distribution functions of one element
//      sum += pdfs[i];
//      if (i != Q19_0) {
//        if (i != Q19_N && i != Q19_S) {
//          tmp_ux += map[GetInvDirection((Direction)i)].off_x*pdfs[i];
//        }
//        if (i != Q19_E && i != Q19_W)
//          tmp_uy += map[GetInvDirection((Direction)i)].off_y*pdfs[i];
//      }
//    }
//
//    tmp_ux = tmp_ux / sum;
//    tmp_uy = tmp_uy / sum;
//
//    sum = 1.0; // the enforced density
//    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
//    tmp_ux = 3.0 * tmp_ux;
//    tmp_uy = 3.0 * tmp_uy;
//
//
//    for (int j = 0; j < N_DIR; j++) {
//      tmp = map[GetInvDirection((Direction)j)].off_x * tmp_ux + map[GetInvDirection((Direction)j)].off_y * tmp_uy;
//      PDF(cur, x, y, j) = sum * weights[j] * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
//    }
//  }

  return;
}


} // end namespace
