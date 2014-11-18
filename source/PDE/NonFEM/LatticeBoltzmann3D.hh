#ifndef LATTICEBOLTZMANN3D_HH_
#define LATTICEBOLTZMANN3D_HH_


#include <set>
#include "Utils/StdVector.hh"
#include "General/Enum.hh"
/** This file is an extract of the LBM code of Markus Wittmann (RRZE) based on the code of Thomas Guess (AM2) based on the code of Georg Pingen (USA).
The original code from Markus can be found at

svn+ssh://eamc080.eam.uni-erlangen.de/home/svn_repo/repository/catalyst/OptiLBM/src

This version makes it simpler and extends to 3D (Fabian Wein) */

namespace CoupledField
{
  // Values in the matrix which are smaller or equal than
  // following value are treated as zeros.
  #define NON_ZERO_TOLERANCE  1.0e-20

  // Number of velocity directions per element, e.g. 19 for D3Q19 model
  #define N_DIR 19

  #ifndef __LX__
    #define __LX__  m_sizeX
  #endif

  #ifndef __LY__
    #define __LY__  m_sizeY
  #endif

  #ifndef __LZ__
    #define __LZ__  m_sizeZ
  #endif


  #ifndef __PDFS__
    #define __PDFS__ m_pdfs
  #endif

    // -- ARRAYS OF STRUCTURES LAYOUT --

    // Address PDFs via Coordinates and Direction.
    // Assuming array of structures layout.
    #define PDF(_arrayIndex, _x, _y, _z,  _direction)  __PDFS__[_arrayIndex][ ((_z) * __LX__ * __LY__  + (_y) * __LX__ + (_x)) * N_DIR + (_direction)]

    // Address PDFs via Index and Direction only.
    // Assuming array of structures layout.
    #define PDF_IDX(_arrayIndex, _index, _direction)  __PDFS__[_arrayIndex][(_index) * N_DIR + (_direction)]

    // -- STRUCTURES OF ARRAYS LAYOUT --

    // Address PDFs via Coordinates and Direction.
    // Assuming structures of array layout.
//     #define PDF(_arrayIndex, _x, _y, _direction)  __PDFS__[_arrayIndex][ ((__LX__) * (__LY__)) * (_direction) +  (_y) * __LX__ + (_x) ]

    // Address PDFs via Index and Direction only.
    // Assuming structures of array layout.
//     #define PDF_IDX(_arrayIndex, _index, _direction)  __PDFS__[_arrayIndex][ ((__LX__) * (__LY__)) * (_direction) + (_index) ]

    // Fluid nodes have values in the range of [0.0; 1.0].
    #define LBM_NODE_TYPE_BB      (-1.0)
    #define LBM_NODE_TYPE_INLET   (-2.0)
    #define LBM_NODE_TYPE_OUTLET  (-3.0)

  class LatticeBoltzmann3D
  {
  public:


    LatticeBoltzmann3D(int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, double omega, int maxIterations, double maxTolerance, bool plot);

    ~LatticeBoltzmann3D();

    /** Performs all the LBM iterations until a steady-state with a given tolerance is reached.
     * @param info stores current and final info there. Saves under way
     * @return pdfs will be subject to coll_step() called from LatticeBoltzmannPDE */
    StdVector<double>* Iterate(const StdVector<double>& elements, PtrParamNode info);

    /*** performs a single propagation step on the current array. Called only by LatticeBoltzmannPDE to prepare for the adjoint calculation */
    void prop_step();

  private:
    struct PropTransform{
      int off_x;
      int off_y;
      int off_z;
      PropTransform(int offx, int offy, int offz): off_x(offx), off_y(offy), off_z(offz){}
      PropTransform(): off_x(0),off_y(0), off_z(0){}
    };

    // Boundary C corresponds to interior
//    typedef enum {C=0, E=1, N=2, W=3, S=4, NE=5, NW=6, SW=7, SE=8} Boundary;
    // Corner C stands for center element
//    typedef enum {C = 0, B_NW = 1, B_SW = 2, B_SE = 3, B_NE = 4, T_NW = 5, T_SW = 6, T_SE = 7, T_NE = 8} Corner;
//    typedef enum {TN = 0, TW = 1, TS = 2, TE = 3, NW = 4, SW = 5, SE = 6, NE = 7, BN = 8, BW = 9, BS = 10, BE = 11} Edge;
    // corners: B_NW to T_NE, edges: TN to BE
    typedef enum {C = 0, BNW = 1, BSW = 2, BSE = 3, BNE = 4, TNW = 5, TSW = 6, TSE = 7, TNE = 8, TN = 9,
                  TW = 10, TS = 11, TE = 12, NW = 13, SW = 14, SE = 15, NE = 16, BN = 17, BW = 18, BS = 19, BE = 20} Boundary;
    typedef enum {Q19_0 = 0, Q19_E = 1, Q19_W = 2, Q19_N = 3, Q19_S = 4, Q19_T = 5, Q19_B = 6,
                  Q19_NE = 7, Q19_SW = 8, Q19_NW = 9, Q19_SE = 10,
                  Q19_TN = 11, Q19_BS = 12, Q19_TS = 13, Q19_BN = 14,
                  Q19_TE = 15, Q19_BW = 16, Q19_TW = 17, Q19_BE = 18} Direction;
//    typedef enum{ Q19_INV_0 = Q19_0, Q19_INV_E = Q19_W, Q19_INV_W = Q19_E, Q19_INV_N = Q19_S, Q19_INV_S = Q19_N, Q19_INV_T = Q19_B, Q19_INV_B = Q19_T,
//                  Q19_INV_NE = Q19_SW, Q19_INV_SW = Q19_NE, Q19_INV_NW = Q19_SE, Q19_INV_SE = Q19_NW,
//                  Q19_INV_TN = Q19_BS, Q19_INV_BS = Q19_TN, Q19_INV_TS = Q19_BS, Q19_INV_BN = Q19_TS,
//                  Q19_INV_TE = Q19_BW, Q19_INV_BW = Q19_TE, Q19_INV_TW = Q19_BE, Q19_INV_BE = Q19_TW} DirectionInv;

    static Enum<Direction> directions;
    static Enum<Boundary> boundaries;
//    static Enum<Corner> corners;
//    static Enum<Edge> edges;
//    static Enum<DirectionInv> directionsInv;

    void SetupDataStructures(const StdVector<double>& elements);

    /** Sets distribution functions to initial value. */
    void InitializePdfs();

    /**
     * returns associated integer value of velocity direction two given principal directions
     * e.g. getIndexDir(S,E) returns Q19_SE
     */
    int GetIndexDir(Direction dir1);
    int GetIndexDir(Direction dir1, Direction dir2);

    /**
     * returns associated integer value of boundary of a cube
     * e.g. getIndexBound(T,N,E) returns right top back corner
     *
     * dir1 can only be top, bottom, north or south
     * dir2 can only be north, south, west or east
     */
    int GetIndexBound(Direction dir1, Direction dir2);

    /**
      dir1 can only be top or bottom
      dir2 can only be north or south
      dir3 can only be east or west
    */
    int GetIndexBound(Direction dir1, Direction dir2, Direction dir3);

    // validates GetIndexbound function
    void TestGetIndexbound();

    // validates GetIndexDir function
    void TestDirectionIndex();

    /**
     * Sets transformation map to handle propagation at boundaries (corners and edges)
     */
    void SetupTransformation();

    /** dump propagation maps */
    void WritePropMap();

    /** set enumerations for directions and boundaries*/
    void SetEnums();

    /** Set LBM weights for D3Q19 model*/
    void InitWeights();

    /** Set lookup table for inverse directions */
    void SetInvDirections();

    /** get inverse direction of D2Q9 direction */
    int GetInvDirection(Direction dir);

    void TestInvDirections();



    /** debug information */
    std::string ToString(const StdVector<StdVector<int> >& data);

    /** fills given propagation map with given default offset values (0 in all directions)*/
    void InitPropMap(StdVector<PropTransform>& map);

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
    int m_sizeZ;
    int m_nNodes;     // Number of total nodes (fluid + obstacle)
    double m_ux;      // Inlet x velocity
    double m_uy;      // Inlet y velocity
    double m_uz;      // Inlet z velocity
    double m_omega;
    int m_maxIter;
    double m_maxTol;
    /** plot the residuum over lbm iterations */
    bool m_plot;

    StdVector<double> Scales;

    StdVector<double> weights;

    StdVector< StdVector<double> > m_pdfs;

//    // contains propagation rules for all corners of the domain
//    StdVector< StdVector<PropTransform> > prop_corners;
//    // contains propagation rules for all edges of the domain
//    StdVector< StdVector<PropTransform> > prop_edges;
    // contains propagation rules for all types of elements (edge, corner, interior)
    StdVector< StdVector<PropTransform> > prop_maps;
    // lookup table to get inverse directions
    StdVector<Direction> directionsInv;
    //double * m_pdfs[2];
    int m_cur;
    int m_next;

    StdVector<StdVector<int> > inlet;
    StdVector<StdVector<int> > outlet;
    StdVector<StdVector<int> > bb;
    StdVector<StdVector<int> > rel; // indices of the fluid m_nodes
  }; // end LatticeBoltzmann3D



} // namespace



#endif /* LATTICEBOLTZMANN_HH_ */
