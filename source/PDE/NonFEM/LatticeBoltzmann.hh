#ifndef LATTICEBOLTZMANN_HH_
#define LATTICEBOLTZMANN_HH_


#include <set>
#include "Utils/StdVector.hh"


/** This file is an extract of the LBM code of Markus Wittmann (RRZE) based on the code of Thomas Guess (AM2) based on the code of Georg Pingen (USA).
The original code from Markus can be found at

svn+ssh://eamc080.eam.uni-erlangen.de/home/svn_repo/repository/catalyst/OptiLBM/src

This version makes it simpler and extends to 3D (Fabian Wein) */

namespace CoupledField
{
  // Values in the matrix which are smaller or equal than
  // following value are treated as zeros.
  #define NON_ZERO_TOLERANCE  1.0e-20

  #define Q9_0  0
  #define Q9_E  1
  #define Q9_N  2
  #define Q9_W  3
  #define Q9_S  4
  #define Q9_NE 5
  #define Q9_NW 6
  #define Q9_SW 7
  #define Q9_SE 8


  #define Q9_INV_0  Q9_0
  #define Q9_INV_E  Q9_W
  #define Q9_INV_N  Q9_S
  #define Q9_INV_W  Q9_E
  #define Q9_INV_S  Q9_N
  #define Q9_INV_NE Q9_SW
  #define Q9_INV_NW Q9_SE
  #define Q9_INV_SW Q9_NE
  #define Q9_INV_SE Q9_NW


  #define T_Q9_0   (4.0 /  9.0)
  #define T_Q9_E   (1.0 /  9.0)
  #define T_Q9_N   (1.0 /  9.0)
  #define T_Q9_W   (1.0 /  9.0)
  #define T_Q9_S   (1.0 /  9.0)
  #define T_Q9_NE  (1.0 / 36.0)
  #define T_Q9_NW  (1.0 / 36.0)
  #define T_Q9_SW  (1.0 / 36.0)
  #define T_Q9_SE  (1.0 / 36.0)

  #ifndef __LX__
    #define __LX__  m_sizeX
  #endif

  #ifndef __LY__
    #define __LY__  m_sizeY
  #endif

  #ifndef __PDFS__
    #define __PDFS__ m_pdfs
  #endif

    // -- ARRAYS OF STRUCTURES LAYOUT --

    // Address PDFs via Coordinates and Direction.
    // Assuming array of structures layout.
    #define PDF(_arrayIndex, _x, _y, _direction)  __PDFS__[_arrayIndex][ ((_y) * __LX__ + (_x)) * 9 + (_direction)]

    // Address PDFs via Index and Direction only.
    // Assuming array of structures layout.
    #define PDF_IDX(_arrayIndex, _index, _direction)  __PDFS__[_arrayIndex][(_index) * 9 + (_direction)]

    // -- STRUCTURES OF ARRAYS LAYOUT --

    // Address PDFs via Coordinates and Direction.
    // Assuming structures of array layout.
    // #define PDF(_arrayIndex, _x, _y, _direction)  __PDFS__[_arrayIndex][ ((__LX__) * (__LY__)) * (_direction) +  (_y) * __LX__ + (_x) ]

    // Address PDFs via Index and Direction only.
    // Assuming structures of array layout.
    // #define PDF_IDX(_arrayIndex, _index, _direction)  __PDFS__[_arrayIndex][ ((__LX__) * (__LY__)) * (_direction) + (_index) ]

    // Fluid nodes have values in the range of [0.0; 1.0].
    #define LBM_NODE_TYPE_BB      (-1.0)
    #define LBM_NODE_TYPE_INLET   (-2.0)
    #define LBM_NODE_TYPE_OUTLET  (-3.0)

  class LatticeBoltzmann
  {
public:


    LatticeBoltzmann(int sizeX, int sizeY, double ux, double uy, double omega, int maxIterations, double maxTolerance, bool plot);

    ~LatticeBoltzmann();

    /** Performs all the LBM iterations until a steady-state with a given tolerance is reached.
     * @param info stores current and final info there. Saves under way
     * @return pdfs will be subject to coll_step() called from LatticeBoltzmannPDE */
    StdVector<double>* Iterate(const StdVector<double>& elements, PtrParamNode info);

    /*** performs a single propagation step on the current array. Called only by LatticeBoltzmannPDE to prepare for the adjoint calculation */
    void prop_step();

  private:

    void SetupDataStructures(const StdVector<double>& elements);


    /** Sets distribution functions to initial value. */
    void InitializePdfs();

    /** debug information */
    std::string ToString(const StdVector<StdVector<int> >& data);


    inline bool LbmNodeTypeIsFluid(double value)
    {
      return (value >= 0.0 && value <= 1.0);
    }

    inline bool LbmNodeTypeIsBB(double value)
    {
      return (value > (LBM_NODE_TYPE_BB - 0.5) && value < (LBM_NODE_TYPE_BB + 0.5));
    }

    inline bool LbmNodeTypeIsInlet(double value)
    {
      return (value > (LBM_NODE_TYPE_INLET - 0.5) && value < (LBM_NODE_TYPE_INLET + 0.5));
    }

    inline bool LbmNodeTypeIsOutlet(double value)
    {
      return (value > (LBM_NODE_TYPE_OUTLET - 0.5) && value < (LBM_NODE_TYPE_OUTLET + 0.5));
    }

    void create_output(const char * file, int cur);


    void prop_coll_step(int m_cur, int m_next, double omega);

    void prop_coll_velinlet(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY);

    void prop_coll_bounce_back(int cur, StdVector<StdVector<int> >& bb);

    void prop_coll_densoutlet(int cur, StdVector<StdVector<int> >&outlet);

    int m_sizeX;
    int m_sizeY;
    int m_nNodes;     // Number of total nodes (fluid + obstacle)
    double m_ux;      // Inlet x velocity
    double m_uy;      // Inlet y velocity
    double m_omega;
    int m_maxIter;
    double m_maxTol;
    /** plot the residuum over lbm iterations */
    bool m_plot;

    StdVector<double> Scales;

    StdVector< StdVector<double> > m_pdfs;

    //double * m_pdfs[2];
    int m_cur;
    int m_next;

    StdVector<StdVector<int> > inlet;
    StdVector<StdVector<int> > outlet;
    StdVector<StdVector<int> > bb;
    StdVector<StdVector<int> > rel; // indices of the fluid m_nodes
  }; // end LatticeBoltzmann



} // namespace



#endif /* LATTICEBOLTZMANN_HH_ */
