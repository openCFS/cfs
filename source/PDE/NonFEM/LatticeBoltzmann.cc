#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "PDE/NonFEM/LatticeBoltzmann.hh"

namespace CoupledField
{

using std::fstream;
using std::ios;

// #define __PDFS__ lattice.Pdfs


DECLARE_LOG(lattice)
DEFINE_LOG(lattice, "lattice")


 LatticeBoltzmann::LatticeBoltzmann(int sizeX, int sizeY, double _ux, double _uy, double omega, int maxIterations, double maxTolerance, const StdVector<double>& elements)
 {
   m_sizeX = sizeX;
   m_sizeY = sizeY;
   m_penalty = 3.0; // penalty;
   m_ux = _ux;
   m_uy = _uy;
   m_density = 1.0; //density;
   m_omega = omega;
   m_maxIter = maxIterations;
   m_maxTol = maxTolerance;
   non_sing_length = 99;


   m_nNodes = m_sizeX * m_sizeY;

   //quadrature weights of LBM
   t.Resize(9);
   set_t();

   //matrix of the probability distributions
   LOG_DBG(lattice) << "Allocating arrays for " << m_nNodes << " PDFs (" << (sizeof(double) * m_nNodes * 9 * 2.0 / 1024.0 / 1024.0) << " MiB)";

   m_pdfs[0] = (double *)malloc(sizeof(double) * m_nNodes * 9);
   m_pdfs[1] = (double *)malloc(sizeof(double) * m_nNodes * 9);

   assert(m_pdfs[0] != NULL);
   assert(m_pdfs[1] != NULL);

   m_lattice.Pdfs[0] = m_pdfs[0];
   m_lattice.Pdfs[1] = m_pdfs[1];

   m_lattice.Scales = (double *)malloc(sizeof(double) * m_nNodes);
   assert(m_lattice.Scales != NULL);

   for (int i = 0; i < m_nNodes; ++i) {
     m_lattice.Scales[i] = 1.0 - pow(por[i], m_penalty);
   }

   m_lattice.SizeX = m_sizeX;
   m_lattice.SizeY = m_sizeY;

   m_cur  = 0;
   m_next = 1;

   InitializePdfs();

   // // indices of the fluid m_nodes
   // rel = new StdVector<StdVector<int> >(1, StdVector<int>(2));

   // // indices of the bounce-back m_nodes
   // bb = new StdVector<StdVector<int> >(1, StdVector<int>(2));
   // // indices of the inlet m_nodes

   // inlet = new StdVector<StdVector<int> >(1, StdVector<int>(2));
   // // indices of the outlet m_nodes
   // outlet = new StdVector<StdVector<int> >(1, StdVector<int>(2));

   // set_data();

   ux.Resize(m_nNodes);
   uy.Resize(m_nNodes);
   dloc.Resize(m_nNodes);

 }

 LatticeBoltzmann::~LatticeBoltzmann()
 {
   free(m_pdfs[0]);
   free(m_pdfs[1]);
 }


void LatticeBoltzmann::Iterate(double * output)
{
    int index;
    int it = 0;
    int count;

    double res;
    double R;
    double refnorm;
    bool nonsteady = true;

    create_output("node.dat", m_cur);
    create_inlet("inlet.dat");

    std::cout << "bb: ";
    for(unsigned int i = 0; i < bb.GetSize(); i++)
      std::cout << bb[i][0] << ":" << bb[i][1] << " ";

    std::cout << "\ninlet: ";
    for(unsigned int i = 0; i < inlet.GetSize(); i++)
      std::cout << inlet[i][0] << ":" << inlet[i][1] << " ";

    std::cout << "\noutlet: ";
    for(unsigned int i = 0; i < outlet.GetSize(); i++)
      std::cout << outlet[i][0] << ":" << outlet[i][1] << " ";

    std::cout << "\nrel: ";
    for(unsigned int i = 0; i < rel.GetSize(); i++)
      std::cout << rel[i][0] << ":" << rel[i][1] << " ";

    std::cout << "\n";
    std::cout.flush();

    while(it < m_maxIter && nonsteady)
    {
      // -- Combined propagation and collision step -------------------------
      prop_coll_step(m_lattice, m_cur, m_next, m_omega);

      // -- Bounce back step ------------------------------------------------
      prop_coll_bounce_back(m_next, bb);

      create_output(std::string("before_inlet_node_" + boost::lexical_cast<std::string>(it) + ".dat").c_str(), m_next);

      // -- Inlet condition -------------------------------------------------
      prop_coll_velinlet(m_next, inlet, m_lbmCase.InletVelocities);

      create_output(std::string("after_inlet_node_" + boost::lexical_cast<std::string>(it) + ".dat").c_str(), m_next);

      // -- Outlet condition ------------------------------------------------
      prop_coll_densoutlet(m_next, outlet, m_density);

      if((it == 0 || it % 100 == 0))
      {
        //Calculation of the residual
        R = 0.;
        count = 0;

        // TODO: parallel computation of residual, take care of non_sing and
        //       summation of residual.
        // TG ORIG: for(int i = 0; i < m_sizeX; i++) {
        // TG ORIG:   for(int j = 0; j < m_sizeY; j++) {
        for(int j = 0; j < m_sizeY; j++) {
          for(int i = 0; i < m_sizeX; i++) {
            index = j * m_sizeX + i;

            for(int k = 0; k < 9; k++)
            {
              res = PDF_IDX(m_next, index, k) - PDF_IDX(m_cur, index, k);
              R += res * res;
              count++;
            }
          }
        }

        R = sqrt(R);

        if(it == 0) {
          refnorm = 1.;

        } else {
          if(R / refnorm <= m_maxTol) {
            nonsteady = false;

          } else if(R >= 1000) {
            std::stringstream ss;
            ss << "In LBM iteration " << it << " residuum " << R << " too large -> abort";
            info->Get("LBM")->Get(ParamNode::WARNING)->SetValue(ss.str());
            nonsteady = false;
          }
        }
      }

      m_cur  = (m_cur  + 1) % 2;
      m_next = (m_next + 1) % 2;

      it++;

    }

    // DumpFluidNodes("fluid-nodes.dat", it, m_cur);

    create_output("node_steady.dat", m_cur);
    // output number of iterations; Residual
    output[0] = --it;
    output[1] = R / refnorm;
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

void LatticeBoltzmann::set_t()
{
  // sets quadrature weights
  t[0] = 4. / 9.;

  for(int i = 1; i < 5; i++) {
    t[i] = 1. / 9.;
  }

  for(int i = 5; i < 9; i++) {
    t[i] = 1. / 36.;
  }
}


void LatticeBoltzmann::SetupDataStructures(const StdVector<double>& elements)
{
  StdVector<int> tmp(2);

  int n = 0;
  double porosity;

  // TODO: maybe interchange loops as x is the fast dimension.

  // TODO: this only works for two dimensions...

  // TODO: preallocate the StdVectors and just set the indices

  for(int i = 0; i < m_sizeX; i++) {
    for(int j = 0; j < m_sizeY; j++) {

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


void LatticeBoltzmann::create_output(const char * file, int cur)
{
  // for debug purposes
  fstream f;
  f.precision(16);
  f.open(file, ios::out);

  for(int i = 0; i < m_nNodes; i++) {
    for(int j = 0; j < 9; j++) {
      f << PDF_IDX(cur, i, j) << " ";
    }

    f << std::endl;
  }

  f.close();
}

void LatticeBoltzmann::create_outlet(const char * file) {
  // for debug purposes
  fstream f;
  f.open(file, ios::out);
  f.precision(16);

  for(unsigned int i = 0; i < outlet.GetSize(); i++) {
    for(unsigned int j = 0; j < 2; j++) {
      f << outlet[i][j] << " ";
    }

    f << std::endl;
  }

  f.close();
}

void LatticeBoltzmann::create_inlet(const char * file)
{
  // for debug purposes
  fstream f;
  f.open(file, ios::out);
  f.precision(16);

  for(unsigned int i = 0; i < inlet.GetSize(); i++) {
    for(unsigned int j = 0; j < 2; j++) {
      f << inlet[i][j] << " ";
    }

    f << std::endl;
  }

  f.close();
}


void LatticeBoltzmann::create_bb(const char * file)
{
  // for debug purposes
  fstream f;
  f.open(file, ios::out);
  f.precision(16);

  for(unsigned int i = 0; i < bb.GetSize(); i++) {
    for(unsigned int j = 0; j < 2; j++) {
      f << bb[i][j] << " ";
    }

    f << std::endl;
  }

  f.close();
}




void LatticeBoltzmann::prop_step(Lattice& lattice, int cur)
{
  // perform a propagation step
  int x, y;

  int lx = lattice.SizeX;
  int ly = lattice.SizeY;

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



void LatticeBoltzmann::prop_coll_step(Lattice& lattice, int m_cur, int m_next, double omega)
{
  // perform a propagation step
  int x, y, index;

  int lx = lattice.SizeX;
  int ly = lattice.SizeY;

  double tmp_ux, tmp_uy, tmp_us, scale, sum;
  double pdf_0, pdf_w, pdf_s, pdf_e, pdf_n, pdf_sw, pdf_se, pdf_ne, pdf_nw;
  double tmp;


  assert(lattice.SizeX > 0);
  assert(lattice.SizeY > 0);
  assert(lattice.Pdfs[0] != NULL);
  assert(lattice.Pdfs[1] != NULL);

  // {{{ Propagation at the boundaries

  // ---------------------------------------------------------------------
  // Propagation for nodes in the corners....
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  // Bottom Left
  x = 0; y = 0;
  PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);

  PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x,     y + 1, Q9_S);
  PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x,     y,     Q9_SE);
  PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x + 1, y + 1, Q9_SW);

  PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x,     y,     Q9_N);
  PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x,     y,     Q9_NE);
  PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x,     y,     Q9_NW);

  PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x + 1, y,     Q9_W);
  PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x,     y,     Q9_E);

  // ---------------------------------------------------------------------
  // Bottom Right
  x = lx - 1; y = 0;
  PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);

  PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x,     y + 1, Q9_S);
  PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x - 1, y + 1,  Q9_SE);
  PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x,     y,     Q9_SW);

  PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x,     y,     Q9_N);
  PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x,     y,     Q9_NE);
  PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x,     y,     Q9_NW);

  PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x,     y,     Q9_W);
  PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x - 1, y,     Q9_E);

  // ---------------------------------------------------------------------
  // Top Left
  x = 0; y = ly - 1;
  PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);

  PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x,     y,     Q9_S);
  PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x,     y,     Q9_SE);
  PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x,     y,     Q9_SW);

  PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x,     y - 1, Q9_N);
  PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x,     y,     Q9_NE);
  PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x + 1, y - 1, Q9_NW);

  PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x + 1, y,     Q9_W);
  PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x,     y,     Q9_E);

  // ---------------------------------------------------------------------
  // Top Right
  x = lx - 1; y = ly - 1;
  PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);

  PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x,     y,     Q9_S);
  PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x,     y,     Q9_SE);
  PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x,     y,     Q9_SW);

  PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x,     y - 1, Q9_N);
  PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x - 1, y - 1, Q9_NE);
  PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x,     y,     Q9_NW);

  PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x,     y,     Q9_W);
  PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x - 1, y,     Q9_E);


  // Propagation along the boundaries of the X dimension.

  for (x = 1; x < lx - 1; ++x) {
    // Bottom
    y = 0;
    PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);
    PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x,     y + 1,      Q9_S);
    PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x - 1, y + 1,      Q9_SE);
    PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x + 1, y + 1,      Q9_SW);

    PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x, y, Q9_N);
    PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x, y, Q9_NE);
    PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x, y, Q9_NW);

    PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x - 1,     y,      Q9_E);
    PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x + 1,     y,      Q9_W);

    // Top
    y = ly - 1;
    PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);
    PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x, y,      Q9_S);
    PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x, y,      Q9_SE);
    PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x, y,      Q9_SW);

    PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x,     y - 1, Q9_N);
    PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x - 1, y - 1, Q9_NE);
    PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x + 1, y - 1, Q9_NW);

    PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x - 1,     y,      Q9_E);
    PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x + 1,     y,      Q9_W);
  }

  // Propagation along the boundaries of the Y dimension.


  for (y = 1; y < ly - 1; ++y) {
    // Left
    x = 0;
    PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);
    PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x,     y + 1, Q9_S);
    PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x,     y,     Q9_SE);
    PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x + 1, y + 1, Q9_SW);

    PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x,     y - 1, Q9_N);
    // ASSERT(PDF(m_next, 0, 1, Q9_N) < 0.0);

    PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x,     y,     Q9_NE);
    PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x + 1, y - 1, Q9_NW);

    PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x,     y,     Q9_E);
    PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x + 1, y,     Q9_W);

    // Rigth
    x = lx - 1;
    PDF(m_next, x, y, Q9_0)  = PDF(m_cur, x,     y,     Q9_0);
    PDF(m_next, x, y, Q9_S)  = PDF(m_cur, x,     y + 1, Q9_S);
    PDF(m_next, x, y, Q9_SE) = PDF(m_cur, x - 1, y + 1, Q9_SE);
    PDF(m_next, x, y, Q9_SW) = PDF(m_cur, x,     y,     Q9_SW);

    PDF(m_next, x, y, Q9_N)  = PDF(m_cur, x,     y - 1, Q9_N);
    PDF(m_next, x, y, Q9_NE) = PDF(m_cur, x - 1, y - 1, Q9_NE);
    PDF(m_next, x, y, Q9_NW) = PDF(m_cur, x,     y,     Q9_NW);

    PDF(m_next, x, y, Q9_E)  = PDF(m_cur, x - 1, y,     Q9_E);
    PDF(m_next, x, y, Q9_W)  = PDF(m_cur, x,     y,     Q9_W);
  }

  // }}}

  double * scales  = lattice.Scales;


  #ifdef _OPENMP
    #pragma omp parallel for collapse(2) \
      default(none), \
      private(index, pdf_0, pdf_w, pdf_s, pdf_e, pdf_n, pdf_sw, pdf_se, pdf_ne, pdf_nw), \
      private(tmp_ux, tmp_uy, tmp_us, scale, sum, tmp, x, y), \
      shared(lx, ly, m_cur, m_next, omega, scales, lattice)
  #endif
  for (y = 1; y < ly - 1; ++y) {
    for (x = 1; x < lx - 1; ++x) {

      index = y * lx + x;
      // ++index;

      // TODO: look at the generated assembly code for addressing the arrays.
      //       If it's too bad we should use a more sophisticated referenceing
      //       scheme.
      pdf_0  = PDF(m_cur, x,     y,     Q9_0);
      pdf_e  = PDF(m_cur, x - 1, y,     Q9_E);
      pdf_n  = PDF(m_cur, x,     y - 1, Q9_N);
      pdf_w  = PDF(m_cur, x + 1, y,     Q9_W);
      pdf_s  = PDF(m_cur, x,     y + 1, Q9_S);
      pdf_ne = PDF(m_cur, x - 1, y - 1, Q9_NE);
      pdf_nw = PDF(m_cur, x + 1, y - 1, Q9_NW);
      pdf_sw = PDF(m_cur, x + 1, y + 1, Q9_SW);
      pdf_se = PDF(m_cur, x - 1, y + 1, Q9_SE);

      sum = pdf_0 + pdf_e + pdf_n + pdf_w + pdf_s + pdf_ne + pdf_nw + pdf_sw + pdf_se;

      // macroscopic scaling by design variable
      // scale = 1.0 - pow(por[index], penal);
      scale = scales[index];

      tmp_ux = scale * (pdf_e + pdf_ne + pdf_se - pdf_w - pdf_nw - pdf_sw) / sum;
      tmp_uy = scale * (pdf_n + pdf_ne + pdf_nw - pdf_s - pdf_sw - pdf_se) / sum;

      tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
      tmp_ux = 3.0 * tmp_ux;
      tmp_uy = 3.0 * tmp_uy;

      PDF(m_next, x, y, Q9_0)  = pdf_0  + omega * ((sum * T_Q9_0  * (1.0 - tmp_us)) - pdf_0);

      tmp = tmp_ux;
      PDF(m_next, x, y, Q9_E)  = pdf_e  + omega * ((sum * T_Q9_E  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_e);

      tmp = tmp_uy;
      PDF(m_next, x, y, Q9_N)  = pdf_n  + omega * ((sum * T_Q9_N  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_n);

      tmp = -tmp_ux;
      PDF(m_next, x, y, Q9_W)  = pdf_w  + omega * ((sum * T_Q9_W  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_w);

      tmp = -tmp_uy;
      PDF(m_next, x, y, Q9_S)  = pdf_s  + omega * ((sum * T_Q9_S  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_s);

      tmp = tmp_ux + tmp_uy;
      PDF(m_next, x, y, Q9_NE) = pdf_ne + omega * ((sum * T_Q9_NE * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_ne);

      tmp = - tmp_ux + tmp_uy;
      PDF(m_next, x, y, Q9_NW) = pdf_nw + omega * ((sum * T_Q9_NW * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_nw);

      tmp = -tmp_ux - tmp_uy;
      PDF(m_next, x, y, Q9_SW) = pdf_sw + omega * ((sum * T_Q9_SW * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_sw);

      tmp = tmp_ux - tmp_uy;
      PDF(m_next, x, y, Q9_SE) = pdf_se + omega * ((sum * T_Q9_SE * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us)) - pdf_se);
    }
  }

  return;
}

void LatticeBoltzmann::prop_coll_velinlet(int cur, StdVector<StdVector<int> >& inlet, StdVector<double>& velocities)
{

  int x, y;
  double tmp_ux, tmp_uy, tmp_us, sum;
  double tmp;

  double pdf_0, pdf_w, pdf_s, pdf_e, pdf_n, pdf_sw, pdf_se, pdf_ne, pdf_nw;

  for(unsigned int i = 0; i < inlet.GetSize(); i++) {
    x = inlet[i][0];
    y = inlet[i][1];

    pdf_0  = PDF(cur, x, y, Q9_0);
    pdf_e  = PDF(cur, x, y, Q9_E);
    pdf_n  = PDF(cur, x, y, Q9_N);
    pdf_w  = PDF(cur, x, y, Q9_W);
    pdf_s  = PDF(cur, x, y, Q9_S);
    pdf_ne = PDF(cur, x, y, Q9_NE);
    pdf_nw = PDF(cur, x, y, Q9_NW);
    pdf_sw = PDF(cur, x, y, Q9_SW);
    pdf_se = PDF(cur, x, y, Q9_SE);

    sum = pdf_0 + pdf_e + pdf_n + pdf_w + pdf_s + pdf_ne + pdf_nw + pdf_sw + pdf_se;

    tmp_ux = velocities[i * 3 + 0];
    tmp_uy = velocities[i * 3 + 1];
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;

    std::cout << "prop_coll_velinlet new: i=" << i << " ux_org=" << velocities[i * 3 + 0] << " uy_org=" << velocities[i * 3 + 1] << " tux=" << tmp_ux << " tuy=" << tmp_uy << std::endl;

    PDF(cur, x, y, Q9_0)  = sum * T_Q9_0  * (1.0 - tmp_us);

    tmp = tmp_ux;
    PDF(cur, x, y, Q9_E)  = sum * T_Q9_E  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = tmp_uy;
    PDF(cur, x, y, Q9_N)  = sum * T_Q9_N  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = -tmp_ux;
    PDF(cur, x, y, Q9_W)  = sum * T_Q9_W  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = -tmp_uy;
    PDF(cur, x, y, Q9_S)  = sum * T_Q9_S  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = tmp_ux + tmp_uy;
    PDF(cur, x, y, Q9_NE) = sum * T_Q9_NE * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = - tmp_ux + tmp_uy;
    PDF(cur, x, y, Q9_NW) = sum * T_Q9_NW * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = -tmp_ux - tmp_uy;
    PDF(cur, x, y, Q9_SW) = sum * T_Q9_SW * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = tmp_ux - tmp_uy;
    PDF(cur, x, y, Q9_SE) = sum * T_Q9_SE * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

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

  double pdf_w, pdf_s, pdf_e, pdf_n, pdf_sw, pdf_se, pdf_ne, pdf_nw;

  for(unsigned int i = 0; i < bb.GetSize(); i++) {
    x = bb[i][0];
    y = bb[i][1];

    pdf_e  = PDF(cur, x, y, Q9_E);
    pdf_n  = PDF(cur, x, y, Q9_N);
    pdf_w  = PDF(cur, x, y, Q9_W);
    pdf_s  = PDF(cur, x, y, Q9_S);
    pdf_ne = PDF(cur, x, y, Q9_NE);
    pdf_nw = PDF(cur, x, y, Q9_NW);
    pdf_sw = PDF(cur, x, y, Q9_SW);
    pdf_se = PDF(cur, x, y, Q9_SE);

    PDF(cur, x, y, Q9_INV_E)  = pdf_e ;
    PDF(cur, x, y, Q9_INV_N)  = pdf_n ;
    PDF(cur, x, y, Q9_INV_W)  = pdf_w ;
    PDF(cur, x, y, Q9_INV_S)  = pdf_s ;
    PDF(cur, x, y, Q9_INV_NE) = pdf_ne;
    PDF(cur, x, y, Q9_INV_NW) = pdf_nw;
    PDF(cur, x, y, Q9_INV_SW) = pdf_sw;
    PDF(cur, x, y, Q9_INV_SE) = pdf_se;
  }

  return;
}

//
// Density outlet condition.
//
void LatticeBoltzmann::prop_coll_densoutlet(int cur, StdVector<StdVector<int> >& outlet, double density)
{

  double tmp_ux, tmp_uy, tmp_us, sum, tmp;
  double pdf_0, pdf_w, pdf_s, pdf_e, pdf_n, pdf_sw, pdf_se, pdf_ne, pdf_nw;

  int x;
  int y;

  for(unsigned int i = 0; i < outlet.GetSize(); i++) {
    x = outlet[i][0];
    y = outlet[i][1];

    pdf_0  = PDF(cur, x, y, Q9_0);
    pdf_e  = PDF(cur, x, y, Q9_E);
    pdf_n  = PDF(cur, x, y, Q9_N);
    pdf_w  = PDF(cur, x, y, Q9_W);
    pdf_s  = PDF(cur, x, y, Q9_S);
    pdf_ne = PDF(cur, x, y, Q9_NE);
    pdf_nw = PDF(cur, x, y, Q9_NW);
    pdf_sw = PDF(cur, x, y, Q9_SW);
    pdf_se = PDF(cur, x, y, Q9_SE);


    sum = pdf_0 + pdf_e + pdf_n + pdf_w + pdf_s + pdf_ne + pdf_nw + pdf_sw + pdf_se;

    tmp_ux = (pdf_e + pdf_ne + pdf_se - pdf_w - pdf_nw - pdf_sw) / sum;
    tmp_uy = (pdf_n + pdf_ne + pdf_nw - pdf_s - pdf_sw - pdf_se) / sum;

    sum = density;
    tmp_us = 1.5 * (tmp_ux * tmp_ux + tmp_uy * tmp_uy);
    tmp_ux = 3.0 * tmp_ux;
    tmp_uy = 3.0 * tmp_uy;


    PDF(cur, x, y, Q9_0)  = sum * T_Q9_0  * (1.0 - tmp_us);

    tmp = tmp_ux;
    PDF(cur, x, y, Q9_E)  = sum * T_Q9_E  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = tmp_uy;
    PDF(cur, x, y, Q9_N)  = sum * T_Q9_N  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = -tmp_ux;
    PDF(cur, x, y, Q9_W)  = sum * T_Q9_W  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = -tmp_uy;
    PDF(cur, x, y, Q9_S)  = sum * T_Q9_S  * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = tmp_ux + tmp_uy;
    PDF(cur, x, y, Q9_NE) = sum * T_Q9_NE * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = - tmp_ux + tmp_uy;
    PDF(cur, x, y, Q9_NW) = sum * T_Q9_NW * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = -tmp_ux - tmp_uy;
    PDF(cur, x, y, Q9_SW) = sum * T_Q9_SW * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);

    tmp = tmp_ux - tmp_uy;
    PDF(cur, x, y, Q9_SE) = sum * T_Q9_SE * (1.0 + tmp + 0.5 * tmp * tmp - tmp_us);
  }

  return;
}


} // end namespace
