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
  // Fluid nodes have values in the range of [0.0; 1.0].
#define LBM_NODE_TYPE_BB      (-1.0)
#define LBM_NODE_TYPE_INLET   (-2.0)
#define LBM_NODE_TYPE_OUTLET  (-3.0)

  class LatticeBoltzmannBase
  {
    public:

      typedef enum {C = 0, BNW = 1, BSW = 2, BSE = 3, BNE = 4, TNW = 5, TSW = 6, TSE = 7, TNE = 8, TN = 9,
              TW = 10, TS = 11, TE = 12, NW = 13, SW = 14, SE = 15, NE = 16, BN = 17, BW = 18, BS = 19, BE = 20, W, E, T, B, S, N} Boundary;
      typedef enum {Q_0 = 0, Q_E = 1, Q_W = 2, Q_N = 3, Q_S = 4, Q_T = 5, Q_B = 6,
                Q_NE = 7, Q_SW = 8, Q_NW = 9, Q_SE = 10,
                Q_TN = 11, Q_BS = 12, Q_TS = 13, Q_BN = 14,
                Q_TE = 15, Q_BW = 16, Q_TW = 17, Q_BE = 18} Direction;
  };

  class LatticeBoltzmann3D
  {
    public:


      LatticeBoltzmann3D(int dim, int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency);

      ~LatticeBoltzmann3D();

      /** Performs all the LBM iterations until a steady-state with a given tolerance is reached.
       * @param info stores current and final info there. Saves under way
       * @return pdfs will be subject to coll_step() called from LatticeBoltzmannPDE */
      StdVector<double>* Iterate(const StdVector<double>& elements, PtrParamNode info);

      /*** performs a single propagation step on the current array. Called only by LatticeBoltzmannPDE to prepare for the adjoint calculation */
      void prop_step();

      /** Returns a copy of current pdf array for calculations of macroscopic values in LatticeBoltzmannPDE during the Iterate() function
       *  @retun copy of pdfs
       */
      inline StdVector<double> GetPdfs() {return m_pdfs[m_next];}

      /**
       * returns number of simulations results we have already written out. We need this to know which number the StoreResults() for the converged solution in staticDriver gets
       */
      inline int GetNumWriteResults() { return m_numWriteResults; }

    private:
      struct PropTransform{
        int off_x;
        int off_y;
        int off_z;
        PropTransform(int offx, int offy, int offz): off_x(offx), off_y(offy), off_z(offz){}
        PropTransform(): off_x(0),off_y(0), off_z(0){}
      };

      // Corner C stands for interior elements
      // BNW, ..., TNE are corners
      // TN, ..., BE are edges
      // LEFT, ..., BACK are faces
      typedef enum {C = 0, BNW = 1, BSW = 2, BSE = 3, BNE = 4, TNW = 5, TSW = 6, TSE = 7, TNE = 8, TN = 9,
        TW = 10, TS = 11, TE = 12, NW = 13, SW = 14, SE = 15, NE = 16, BN = 17, BW = 18, BS = 19, BE = 20, W, E, T, B, S, N} Boundary;
        typedef enum {Q_0 = 0, Q_E = 1, Q_W = 2, Q_N = 3, Q_S = 4, Q_T = 5, Q_B = 6,
          Q_NE = 7, Q_SW = 8, Q_NW = 9, Q_SE = 10,
          Q_TN = 11, Q_BS = 12, Q_TS = 13, Q_BN = 14,
          Q_TE = 15, Q_BW = 16, Q_TW = 17, Q_BE = 18} Direction;

          static Enum<Direction> directions;
          static Enum<Boundary> boundaries;

          void SetupDataStructures(const StdVector<double>& elements);

          /** Sets distribution functions to initial value. */
          void InitializePdfs();

          /**
           * returns associated integer value of velocity direction two given principal directions
           * e.g. getIndexDir(S,E) returns Q_SE
           */
          int GetIndexDir(Direction dir1);
          int GetIndexDir(Direction dir1, Direction dir2);

          /**
           * returns associated integer value of boundary of a cube
           * e.g. getIndexBound(T) returns the index of the top face
           *
           */
          int GetIndexBound(Direction dir1);


          /**
           * returns associated integer value of boundary of a cube
           * e.g. getIndexBound(T,N,E) returns the index of the right top back corner
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

          // make sure that no-slip b.c. is imposed on the boundaries (bounce-back nodes)
          void CheckBoundaryVelocities(int cur, StdVector<StdVector<int> >& bb);

          /**
           * Calculates x, y and z velocity for element with coordinate (i,j,k)
           */
          void CalcVelocitites(int cur, int i, int j, int k, double& ux, double& uy, double& uz);

          /**
           * Calculates macroscopic density for given element
           */
          double CalcDensity(int cur, int i, int j, int k);

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

          // calculates array index for given grid coordinate in 3D
          inline int GetIndex(int x, int y, int z) const
          {
            return z * m_sizeX * m_sizeY + y * m_sizeX + x;
          }

          // calculates array index for given grid coordinate in 2D
          inline int GetIndex(int x, int y) const
          {
            return y * m_sizeX + x;
          }

          // calculates index of neighbor element we are propagating to (depends on the propagation direction)
          // @index: index of the current element
          //    int CalcPropIndex(int index, int dir);

          //--------------- array of structures----------------------------------
          inline double& pdf(int cur, int x, int y, int z, int direction)
          {
            return m_pdfs[cur][direction + n_q_ * GetIndex(x, y, z)];
          }

          inline const double pdf(int cur, int x, int y, int z, int direction) const
          {
            return m_pdfs[cur][direction + n_q_ * GetIndex(x, y, z)];
          }

          inline double& pdf(int cur, int elem, int direction)
          {
            return m_pdfs[cur][direction + n_q_ * elem];
          }

          inline const double pdf(int cur, int elem, int direction) const
          {
            return m_pdfs[cur][direction + n_q_ * elem];
          }

          //---------------------------------------------------------------------

          //--------------- structure of arrays----------------------------------
          //    inline const double pdf(int cur, int x, int y, int z, int direction) const
          //    {
          //      return m_pdfs[cur][m_sizeX * m_sizeY * m_sizeZ * direction + index(x, y, z)];
          //    }

          //    inline double& pdf(int cur, int x, int y, int z, int direction)
          //    {
          //      return m_pdfs[cur][m_sizeX * m_sizeY * m_sizeZ * direction + index(x,y,z)];
          //    }
          //
          //    inline double& pdf(int cur, int index, int direction)
          //    {
          //      return m_pdfs[cur][m_sizeX * m_sizeY * m_sizeZ * direction + index];
          //    }

          //    inline const double pdf(int cur, int index, int direction) const
          //    {
          //      return m_pdfs[cur][m_sizeX * m_sizeY * m_sizeZ * direction + index];
          //    }
          //---------------------------------------------------------------------




          void create_output(const char * file, int cur);

          void prop_coll_step(int cur, int next, double omega);

          void prop_coll_velinlet(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY, double UZ);

          void prop_coll_bounce_back(int cur, StdVector<StdVector<int> >& bb);

          void prop_coll_densoutlet(int cur, StdVector<StdVector<int> >&outlet);

          int m_dim;
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
          int m_writeFrequency;
          // counts how many intermediate steps we have already written to hdf5 file
          int m_numWriteResults;

          // number of discrete velocities in LBM model, e.g. 9 for D2Q19 or 19 for D3Q19
          int n_q_;

          StdVector<double> Scales;

          StdVector<double> weights;

          StdVector< StdVector<double> > m_pdfs;

          // contains propagation rules for all types of elements (edge, corner, interior)
          StdVector< StdVector<PropTransform> > prop_maps;
          // stores microscopic velocities (directions) of D3Q19 model: e.g. for Q_N: e_N = (0,1,0)
          StdVector<PropTransform> velocityDirections;
          // lookup table to get inverse directions
          StdVector<Direction> directionsInv;
          //double * m_pdfs[2];
          int m_cur;
          int m_next;


          StdVector<StdVector<int> > inlet;
          StdVector<StdVector<int> > outlet;
          StdVector<StdVector<int> > bb;
          StdVector<StdVector<int> > rel; // indices of the fluid m_nodes

          ResultHandler* rh = NULL;
  }; // end LatticeBoltzmann3D



} // namespace



#endif /* LATTICEBOLTZMANN_HH_ */
