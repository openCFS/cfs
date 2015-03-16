#ifndef LATTICEBOLTZMANN_HH_
#define LATTICEBOLTZMANN_HH_


#include <set>
#include "Utils/StdVector.hh"
#include "General/Enum.hh"
/** This file is an extract of the LBM code of Markus Wittmann (RRZE) based on the code of Thomas Guess (AM2) based on the code of Georg Pingen (USA).
The original code from Markus can be found at

svn+ssh://eamc080.eam.uni-erlangen.de/home/svn_repo/repository/catalyst/OptiLBM/src

This version makes it simpler and extends to 3D (Fabian Wein) */

namespace CoupledField
{
#define LBM_NODE_TYPE_BB      (-1.0)
#define LBM_NODE_TYPE_INLET   (-2.0)
#define LBM_NODE_TYPE_OUTLET  (-3.0)

  class LatticeBoltzmannBase
  {
    public:

        // In 2D: 9 microscopic directions
        // In 3D: 19 microscopic directions
        typedef enum {Q_0=0, Q_E=1, Q_N=2, Q_W=3, Q_S=4, Q_NE=5, Q_NW=6, Q_SW=7, Q_SE=8,          // 2D
          Q_T = 9, Q_B = 10, Q_TN = 11, Q_BS = 12, Q_TS = 13, Q_BN = 14,                          // 3D
          Q_TE = 15, Q_BW = 16, Q_TW = 17, Q_BE = 18} Direction;                                  // 3D
  };

  class LatticeBoltzmann: LatticeBoltzmannBase
  {
    public:
      struct LatticeVector{
        int off_x;
        int off_y;
        int off_z;
        LatticeVector(int offx, int offy, int offz): off_x(offx), off_y(offy), off_z(offz){}
        LatticeVector(int offx, int offy): off_x(offx), off_y(offy), off_z(0){}
        LatticeVector(): off_x(1),off_y(0), off_z(0){}
      };

      LatticeBoltzmann(int dim, int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency);

      ~LatticeBoltzmann();

      /** Performs all the LBM iterations until a steady-state with a given tolerance is reached.
       * @param info stores current and final info there. Saves under way
       * @return pdfs will be subject to coll_step() called from LatticeBoltzmannPDE */
      StdVector<double>* Iterate(const StdVector<double>& elements, PtrParamNode info);

      /*** performs a single propagation step on the current array. Called only by LatticeBoltzmannPDE to prepare for the adjoint calculation */
      void prop_step();

      /** Returns a copy of current pdf array for calculations of macroscopic values in LatticeBoltzmannPDE during the Iterate() function
       *  @retun copy of pdfs
       */
      inline StdVector<double> GetPdfs() {return m_pdfs[m_cur];}

      /**
       * returns number of simulations results we have already written out. We need this to know which number the StoreResults() for the converged solution in staticDriver gets
       */
      inline int GetNumWriteResults() { return m_numWriteResults; }

      inline StdVector<LatticeVector>* GetVelocityDirections() { return &microDirections; }
      inline StdVector<Direction>* GetInverseDirections() { return &inverseDirections; }

    private:

          static Enum<Direction> directions;

          void SetupDataStructures(const StdVector<double>& elements);

          /**
           * Sets distribution functions to initial value, independently from dimension of problem (no if-statement necessary)
           */
          void InitializePdfs();

          /**
           * Calculates x, y and z velocity for element with coordinate (i,j,k)
           */
          void CalcVelocitites(int cur, int i, int j, int k, double& ux, double& uy, double& uz);

          /**
           * Calculates macroscopic density for given element
           */
          double CalcDensity(int cur, int i, int j, int k);

          /** set enumerations for directions and boundaries*/
          void SetEnums();

          /** Set LBM weights for DxQy model*/
          void InitWeights();

          /** Set lookup table for inverse directions */
          void SetInvDirections();

          /** Set basis vectors of DmQn model */
          void SetMicroVelocities();

          /** get inverse direction of D2Q9 direction */
          // depends on numbering of directions
          inline int GetInvDirection(Direction dir)
          {
            assert(directions.IsValid(dir));
            return inverseDirections[dir];
          }

          void TestInvDirections();

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

          // calculates array index for given grid coordinate in 3D
          inline int GetIndex(int x, int y, int z) const
          {
            return z * m_sizeX * m_sizeY + y * m_sizeX + x;
          }

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

          void create_output(const char * file, int cur);

          inline bool OutsideDomain(int x, int y, int z, int dir)
          {
            LatticeVector tmp = microDirections[dir];
            int tmp_x = x + tmp.off_x;
            int tmp_y = y + tmp.off_y;
            int tmp_z = z + tmp.off_z;
            return (tmp_x < 0 || tmp_x >= m_sizeX || tmp_y < 0 || tmp_y >= m_sizeY || tmp_z < 0 || tmp_z >= m_sizeZ) ;
          }

          /**
           * LBM operators in 2D
           */

          void prop_coll_step2D(int cur, int next, double omega);

          void prop_coll_velinlet2D(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY, double UZ);

          void prop_coll_bounce_back2D(int cur, StdVector<StdVector<int> >& bb);

          void prop_coll_densoutlet2D(int cur, StdVector<StdVector<int> >&outlet);


          /**
           * LBM operators in 3D
           */
          void prop_coll_step3D(int cur, int next, double omega);

          void prop_coll_velinlet3D(int cur, StdVector<StdVector<int> >& inlet, double UX, double UY, double UZ);

          void prop_coll_bounce_back3D(int cur, StdVector<StdVector<int> >& bb);

          void prop_coll_densoutlet3D(int cur, StdVector<StdVector<int> >&outlet);

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
          // indicates whether LBM simulation results should be written to hdf5 file every m_writeFrequency'th step or not
          // set to true, if no optimization is done (e.g. design is prescribed by density file)
          bool writeIntermediateResults;
          int m_writeFrequency;
          // counts how many intermediate steps we have already written to hdf5 file
          int m_numWriteResults;

          // number of microscopic velocities in LBM model, e.g. 9 for D2Q19 or 19 for D3Q19
          int n_q_;

          StdVector<double> Scales;

          StdVector<double> weights;

          StdVector< StdVector<double> > m_pdfs;

          // stores microscopic velocities (directions) of D3Q19 model: e.g. for Q_N: e_N = (0,1,0)
          StdVector<LatticeVector> microDirections;
          // lookup table to get inverse directions
          StdVector<Direction> inverseDirections;
          int m_cur;
          int m_next;


          StdVector<StdVector<int> > inlet;
          StdVector<StdVector<int> > outlet;
          StdVector<StdVector<int> > bb;
          StdVector<StdVector<int> > rel; // indices of the fluid m_nodes

          ResultHandler* rh = NULL;

          // function pointers to LBM operators (propagation, collision); use these to avoid many if-statements to distinguish 2D from 3D case
          void (LatticeBoltzmann::*prop_coll_step)(int, int, double);
          void (LatticeBoltzmann::*prop_coll_velinlet)(int, StdVector<StdVector<int> >&, double, double, double);
          void (LatticeBoltzmann::*prop_coll_bounce_back)(int, StdVector<StdVector<int> >&);
          void (LatticeBoltzmann::*prop_coll_densoutlet)(int, StdVector<StdVector<int> >&);

  }; // end LatticeBoltzmann


} // namespace



#endif /* LATTICEBOLTZMANN_HH_ */
