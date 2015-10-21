#ifndef LATTICEBOLTZMANN_HH_
#define LATTICEBOLTZMANN_HH_


#include <set>
#include "Utils/StdVector.hh"
#include "MatVec/Matrix.hh"
#include "General/Enum.hh"
#include "DataInOut/ResultHandler.hh"
/** This file is an extract of the LBM code of Markus Wittmann (RRZE) based on the code of Thomas Guess (AM2) based on the code of Georg Pingen (USA).
The original code from Markus can be found at

svn+ssh://eamc080.eam.uni-erlangen.de/home/svn_repo/repository/catalyst/OptiLBM/src

An generalization and extension to 3D of the SRT model was added.
The MRT model was added in correspondance to the descriptions made in the paper of Geng Liu et al. (2014)
*/

namespace CoupledField
{
#define LBM_NODE_TYPE_BB      (-1.0)
#define LBM_NODE_TYPE_INLET   (-2.0)
#define LBM_NODE_TYPE_OUTLET  (-3.0)
#define LBM_NODE_TYPE_OBSTACLE  (1.0)

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
      struct PDFDirectionVector{
        int off_x;
        int off_y;
        int off_z;
        PDFDirectionVector(int offx, int offy, int offz): off_x(offx), off_y(offy), off_z(offz){}
        PDFDirectionVector(int offx, int offy): off_x(offx), off_y(offy), off_z(0){}
        PDFDirectionVector(): off_x(1),off_y(0), off_z(0){}
      };

      LatticeBoltzmann(int dim, int sizeX, int sizeY, int sizeZ, double ux, double uy, double uz, StdVector< StdVector<double> > uin, double omega, int maxIterations, double maxTolerance, bool plot, int writeFrequency, bool srt, double omega_e, double omega_eps, double omega_q, double alpha_max);

      ~LatticeBoltzmann();

      /** Performs all the LBM iterations until a steady-state with a given tolerance is reached.
       * @param info stores current and final info there. Saves under way
       * @param niter stores how many iterations til convergence for one simulation call is needed
       * @return pdfs will be subject to coll_step() called from LatticeBoltzmannPDE */
      StdVector<double>* Iterate(const StdVector<double>& elements, PtrParamNode info);

      /** Performs all the adjoint LBM iterations until a steady-state with a given tolerance is reached.
       * @param info stores current and final info there */
      StdVector<double>* IterateAdjoint(PtrParamNode info);

      /*** performs a single propagation step on the current array. Called only by LatticeBoltzmannPDE to prepare for the adjoint calculation */
      void Prop_step();

      /** Returns a copy of current pdf array for calculations of macroscopic values in LatticeBoltzmannPDE during the Iterate() function
       *  @return copy of pdfs
       */
      inline StdVector<double>& GetPdfs() {return pdfs_[cur_];}

      /** Returns a copy of current pdf array for calculations of macroscopic values in LatticeBoltzmannPDE during the Iterate() function
       *  @return copy of pdfs
       */
      inline StdVector<double>& GetAdjPdfs() {return adjMoments_[adjCur_];}

      /** Returns overall dissipation calculated with current pdfs */
      inline double GetDissipation()
      {
        Vector<double> pdfs(n_q_);
        double dissipation = 0.0;
        for (int elem = 0; elem < nNodes_; elem++) {
          for (int dir = 0; dir < n_q_; dir++)
            pdfs[dir] = PDF(cur_,elem,dir);

          Vector<double> moms(n_q_);
          transformation.Mult(pdfs,moms);
          Vector<double> eqMoms;
          CalcEquilMoments(moms,eqMoms);

          Vector<double> f1,f2;
          CalcDarcyForce(moms,elem,f1,f2);
          dissipation += CalcDissipation(moms,eqMoms,f1[3],f1[5]);
        }
        return dissipation;
      }

      /**
       * returns number of simulations results we have already written out. We need this to know which number the StoreResults() for the converged solution in staticDriver gets
       */
      inline int GetNumWriteResults() { return numWriteResults_; }

      /**
       * @return Number of iterations until steady-state convergence
       */
      inline int GetNumIterations() {return numIterations_; }

      /** return adjoint transformation matrix for conversion of adjoint pdfs into momentum space*/
      inline const Matrix<double>& GetAdjTransformation() {return adjTransformation; }

      inline StdVector<PDFDirectionVector>& GetPDFDirectionVectors() { return microVelDirections; }
      inline StdVector<Direction>& GetinvPDFDirections() { return invPDFDirections; }

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
          void CalcVelocities(const Vector<double>& pdfs, double& ux, double& uy, double& uz);

          /**
           * Calculates macroscopic density for given element
           */
          double CalcDensity(const Vector<double>& pdfs);

          /** Calculates residual as difference between two iteration solutions.
           *  Parameter adjoint indicates whether residual should be calculated for pdfs or adjoint pdfs array.
           */
          inline double CalcResidual(int cur, int next, bool adjoint)
          {
            double R = 0.0;
            double res = -1.0;

            for (int elem = 0; elem < nNodes_; elem++) {
              //            index = k * m_sizeX * m_sizeY + j * m_sizeX + i;
              for(int  dir = 0; dir < n_q_; dir++) {
                if (adjoint)
                {
                  res = APDF(next, elem, dir) - APDF(cur, elem, dir);
                }
                else
                  res = PDF(next, elem, dir) - PDF(cur, elem, dir);
                R += res * res;
              }
            }
            return sqrt(R);
          }

          /** Calculates the two Darcy force vectors at given node in accordance to te proposed porosity model of Geng Liu et al. (2014)*/
          void CalcDarcyForce(const Vector<double>& moments, int elemId, Vector<double>& f1, Vector<double>& f2);

          /** Calculates resistance coefficient \alpha used in porosity model which is expressed by RAMP */
          inline double CalcResistanceCoeff(int elemId)
          {
//            return alpha_max_ * (1 - scales[elemId]) / (1 + scales[elemId]); //RAMP
            return 1-scales[elemId];
          }

          /** Calculates dissipation contribution of given node */
          double CalcDissipation(const Vector<double>& moments, const Vector<double>& eqMoments, double fx, double fy);

          /** Calculates adjoint collision matrix after solving of primal fluid field */
          void CalcAdjointCollMatrix(int elemId, const Vector<double>& moments);

          /** Calculates sensitivity of dissipation (objective function) with respect to LBM moments */
          void d_diss_d_moments(int elemId, const Vector<double>& moments);

          /** set enumerations for directions and boundaries*/
          void SetEnums();

          /** Set LBM weights for DxQy model*/
          void InitWeights();

          /** Set lookup table for inverse directions */
          void SetInvDirections();

          /** Set basis vectors of DmQn model */
          void SetMicroVelocities();

          /**
           * Initialize transformation matrix M for momentum space and corresponding relaxation rates matrix S
           * M ist orthogonal.
           */
          void InitTransformMatrix();

          /** get inverse direction of D2Q9 direction */
          // depends on numbering of directions
          inline int GetInvDirection(Direction dir)
          {
            assert(directions.IsValid(dir));
            return invPDFDirections[dir];
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

          inline bool LbmNodeIsObstacle(double value)
          {
            return value == LBM_NODE_TYPE_OBSTACLE;
          }

          // calculates array index for given grid coordinate in 3D
          inline int GetIndex(int x, int y, int z) const
          {
            return z * sizeX_ * sizeY_ + y * sizeX_ + x;
          }

          // calculates index in pdf array for given elem id and pdf direction
          inline unsigned int GetPdfIndex(unsigned int id, unsigned int dir) const {
            return id * n_q_ + dir;
          }

          //--------------- array of structures----------------------------------
          inline double& PDF(int cur, int x, int y, int z, int direction)
          {
            return pdfs_[cur][direction + n_q_ * GetIndex(x, y, z)];
          }

          inline const double PDF(int cur, int x, int y, int z, int direction) const
          {
            return pdfs_[cur][direction + n_q_ * GetIndex(x, y, z)];
          }

          inline double& PDF(int cur, int elem, int direction)
          {
            return pdfs_[cur][direction + n_q_ * elem];
          }

          inline const double PDF(int cur, int elem, int direction) const
          {
            return pdfs_[cur][direction + n_q_ * elem];
          }

          inline double& APDF(int cur, int x, int y, int z, int direction)
          {
            return adjMoments_[cur][direction + n_q_ * GetIndex(x, y, z)];
          }

          inline const double APDF(int cur, int x, int y, int z, int direction) const
          {
            return adjMoments_[cur][direction + n_q_ * GetIndex(x, y, z)];
          }

          inline double& APDF(int cur, int elem, int direction)
          {
            return adjMoments_[cur][direction + n_q_ * elem];
          }

          inline const double APDF(int cur, int elem, int direction) const
          {
            return adjMoments_[cur][direction + n_q_ * elem];
          }

          void CreateOutput(const char * file, int cur);

          inline bool PointsToBoundary(int x, int y, int z, int dir)
          {
            PDFDirectionVector tmpDir = microVelDirections[dir];
            int tmp_x = x + tmpDir.off_x;
            int tmp_y = y + tmpDir.off_y;
            int tmp_z = z + tmpDir.off_z;

            return tmp_x < 0 || tmp_x >= sizeX_ || tmp_y < 0 || tmp_y >= sizeY_ || tmp_z < 0 || tmp_z >= sizeZ_;
          }

           /** Calculate vector of equilibrium moments based on current pdf array
            * m_eq the return value */
          void CalcEquilMoments(const Vector<double>&  moments,  Vector<double>& m_eq) const;

          /**
           * LBM operators in 2D
           */
          void Prop_coll_step2D(int cur, int next);

          void Prop_coll_velinlet2D(int cur);

          void Prop_coll_bounce_back2D(int cur);

          void Prop_coll_densoutlet2D(int cur);

          void AdjointCollision(int cur);

          void AdjointPropagation(int cur, int next);

          /**
           * LBM operators in 3D
           */
          void Prop_coll_step3D(int cur, int next);

          void Prop_coll_velinlet3D(int cur);

          void Prop_coll_bounce_back3D(int cur);

          void Prop_coll_densoutlet3D(int cur);

          int dim_;
          int sizeX_;
          int sizeY_;
          int sizeZ_;
          int nNodes_;     // Number of total nodes (fluid + obstacle)
          double ux_;      // Inlet x velocity
          double uy_;      // Inlet y velocity
          double uz_;      // Inlet z velocity
          double omega_nu_;
          int maxIter_;
          double maxTol_;


          /** plot the residuum over lbm iterations */
          bool plot_;
          // indicates whether LBM simulation results should be written to hdf5 file every m_writeFrequency'th step or not
          // set to true, if no optimization is done (e.g. design is prescribed by density file)
          bool writeIntermediateResults_;
          int writeFrequency_;
          // counts how many intermediate steps we have already written to hdf5 file
          int numWriteResults_;
          // how many iterations until steady-state convergence
          int numIterations_;
          // indicates whether SRT or MRT model should be used
          bool srt_;
          // additional relaxation rates for MRT model
          double omega_e_, omega_eps_, omega_q_;
          double alpha_max_; // parameter used in porosity model of MRT simulation

          // number of microscopic velocities in LBM model, e.g. 9 for D2Q19 or 19 for D3Q19
          int n_q_;

          int lbmCalls_; //counts how often solver was called

          StdVector<double> scales;

          StdVector<double> weights;

          StdVector< StdVector<double> > u_in; // inflow x-velocities in case of parabolic profile

          StdVector< StdVector<double> > pdfs_;
          StdVector< StdVector<double> > adjMoments_;
          StdVector<double> tmpPdfs_;

          // store moments and equilibrium moments of steady state solution
          // need this for adjoint LBM simulation
//          StdVector<double> moments_;
//          StdVector<double> eqMoments_;

          // stores microscopic velocities (directions) of D3Q19 model: e.g. for Q_N: e_N = (0,1,0)
          StdVector<PDFDirectionVector> microVelDirections;
          // lookup table to get inverse directions of the pdfs
          StdVector<Direction> invPDFDirections;
          int cur_, adjCur_;
          int next_, adjNext_;

          // Transformation matrix M for momentum space
          Matrix<double> transformation;
          Matrix<double> invTransformation;
          Matrix<double> adjTransformation;
          // Store multiplication of backtransformation M^-1 with relaxation rates matrix S
          Matrix<double> invM_S;
          // Relaxation rates matrix S is diagonal, thus we only store the diagonal entries
          StdVector<double> relax_rates;

          StdVector<int> inlet;
          StdVector<int> outlet;
          StdVector<int> bb;
          StdVector<int> rel; // indices of the fluid m_nodes
          StdVector<int > obst; // indices of obstacle nodes
          StdVector<Matrix<double> > adjCollision; // adjoint collision matrices
          StdVector<Vector<double> > d_diss_d_m; // partial derivatives of dissipation function w.r.t moments

          // function pointers to LBM operators (propagation, collision); use these to avoid many if-statements to distinguish 2D from 3D case
          void (LatticeBoltzmann::*prop_coll_step)(int, int);
          void (LatticeBoltzmann::*prop_coll_velinlet)(int);
          void (LatticeBoltzmann::*prop_coll_bounce_back)(int);
          void (LatticeBoltzmann::*prop_coll_densoutlet)(int);

  }; // end LatticeBoltzmann


} // namespace



#endif /* LATTICEBOLTZMANN_HH_ */
