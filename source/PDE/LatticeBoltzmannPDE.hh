// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_EXTLBMPDE
#define FILE_EXTLBMPDE

#include <string>
#include <boost/numeric/ublas/matrix_sparse.hpp>

//#include "DataInOut/ParamHandling/ParamNode.hh"
//#include "General/Enum.hh"
//#include "General/defs.hh"
//#include "General/Environment.hh"
//#include "MatVec/Matrix.hh"
//#include "MatVec/Vector.hh"
#include "PDE/LatticeBoltzmannSolver/LatticeBoltzmann.hh"
#include "SinglePDE.hh"
//#include "Utils/nodestoresol.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{
class BaseResult;
class BaseBDBInt;
class Grid;
class PDECoupling;
class DesignElement;
class Function;
class LatticeBoltzmann;
class TransferFunction;
class Function;
class DesignSpace;

using boost::numeric::ublas::compressed_matrix;
using boost::numeric::ublas::mapped_matrix;
using boost::numeric::ublas::generalized_vector_of_vector;


//! Class for mechanic equation (no adaptivity)
class LatticeBoltzmannPDE: public SinglePDE , LatticeBoltzmannBase
{

public:
  typedef enum { INTERNAL, EXTERNAL } Iface;

  //!  Constructor. here we read integration parameters
  /*!
      \param aGrid pointer to grid
   */
  LatticeBoltzmannPDE( Grid* grid, PtrParamNode pn, PtrParamNode infoNode, shared_ptr<SimState> simState, Domain* domain );

  ~LatticeBoltzmannPDE();

  //    //!  Deconstructor
  //    virtual ~ExtLBMPDE() {;};

  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators();

  /** @see SinglePDE::DefineSurfaceIntegrators() */
  void DefineSurfaceIntegrators() {} ;

  /** @see SinglePDE::DefineNcIntegrators() */
  void DefineNcIntegrators() {} ;

  virtual std::map<SolutionType, shared_ptr<FeSpace> > CreateFeSpaces( const std::string&  formulation, PtrParamNode infoNode );


  //! define the SoltionStep-Driver
  virtual void DefineSolveStep();

  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);

  //! initialize time stepping: nothing to do in smoother!
//  void InitTimeStepping();

  //! set time step
  //! \param dt Current time step
//  virtual void SetTimeStep(const Double dt){};

  //! calculate coupling terms
//  virtual void CalcOutputCoupling();

  /** actually calls LBM. */
  void Solve();

  /** Sensitivity analysis (gradient) necessary for the optimization.
    Performs the missing propagation step for the sensitivity analysis. */
  void SensitivityAnalysis(TransferFunction* tf, Function* f, DesignSpace* space);
//
//  /** For the olt multi-file interface, reads the files x_pardiso2.dat and dRdp.mtx and computes the gradient */
//  void SetPrecalculatedGradient(StdVector<DesignElement*>& design, Function* f);

  /** Perform postprocessing on data. */
  // Deprecated
//  void CalcResults( shared_ptr<BaseResult> result );

  /** implementation of objective function */
  double CalcPressureDrop();

  // returns how often CalcResults() was called. We need this to write out the last simulation step
  inline int GetNumWriteResults() {return numWriteResults_;}

  // returns how many iterations were required by LBM solver to reach steady-state in last LBM call
  inline int GetNumIterations() {return numIterations_;}

  // return respective id of element with elemId in LBM world
  inline int GetLbmId(unsigned int elemId) { return elem_to_idx[elemId]; }

  virtual void ReadSpecialResults();

  std::string ToString(const StdVector<double>& elements, bool x_fast, bool as_int) const;

  /** exports a boost compressed matrix in Matrix-Market format */
  // static void ToFile(const std::string& file, const compressed_matrix<double>& M);
  static void ToFile(const std::string& file, const mapped_matrix<double>& M);

  void create_output(const char * file);

  Iface GetIface() const { return iface_; }

  //////////////////////////////////////////////////////// functions calculating results from PDFs //////////////////////////////////////////////
  /** Calculate the LBM Density of an element idx */
  inline double CalcLBMDensity(unsigned int idx) const
  {
    double density = 0.0;
    for(unsigned int h = 0; h < n_q_; h++)
      density += GetPdf(idx,h);

    return density;
  }

  /** Calculate pressure for given element idx */
  double CalcPressure(unsigned int idx) const;

  /** Calculate velocity components for given density and element idx */
  inline double CalcVelocityX(unsigned int idx, double density) const
  {
    if (n_q_ == 9) {
      double val = ((GetPdf(idx, Q_E) - GetPdf(idx, Q_W)) + (GetPdf(idx, Q_NE) - GetPdf(idx, Q_SW)) + (GetPdf(idx, Q_SE) - GetPdf(idx, Q_NW))) / density;
      double term1 = GetPdf(idx, Q_E) - GetPdf(idx, Q_W);
      double term2 = GetPdf(idx, Q_NE) - GetPdf(idx, Q_SW);
      double term3 = GetPdf(idx, Q_SE) - GetPdf(idx, Q_NW);
      double test = term1+term2+term3;
      return val + test - test;
    }
    else
      return ((GetPdf(idx, Q_E) - GetPdf(idx, Q_W)) +( GetPdf(idx, Q_NE) - GetPdf(idx, Q_SW)) + (GetPdf(idx, Q_SE) - GetPdf(idx, Q_NW)) + (GetPdf(idx, Q_TE) - GetPdf(idx, Q_BW)) + (GetPdf(idx, Q_BE) - GetPdf(idx, Q_TW))) / density;
  }

  inline double CalcVelocityY(unsigned int idx, double density) const
  {
    if (n_q_ == 9)
      return ((GetPdf(idx, Q_N) - GetPdf(idx, Q_S)) + (GetPdf(idx, Q_NE) - GetPdf(idx, Q_SW)) + (GetPdf(idx, Q_NW) - GetPdf(idx, Q_SE))) / density;
    else
      return ((GetPdf(idx, Q_N) - GetPdf(idx, Q_S)) + (GetPdf(idx, Q_NW) - GetPdf(idx, Q_SE)) + (GetPdf(idx, Q_NE) - GetPdf(idx, Q_SW)) + (GetPdf(idx, Q_TN) - GetPdf(idx, Q_BS)) + (GetPdf(idx, Q_BN) - GetPdf(idx, Q_TS))) / density;
  }

  inline double CalcVelocityZ(unsigned int idx, double density) const
  {
    if (n_q_ == 9)
      return 0;
    else
      return ((GetPdf(idx, Q_T) - GetPdf(idx, Q_B)) + (GetPdf(idx, Q_TW) - GetPdf(idx, Q_BE)) + (GetPdf(idx, Q_TE) - GetPdf(idx, Q_BW)) + (GetPdf(idx, Q_TN) - GetPdf(idx, Q_BS)) + (GetPdf(idx, Q_TS) - GetPdf(idx, Q_BN))) / density;
  }

  inline void ExtractIntermediateSolution() {
    pdfs_ = lbm->GetPdfs();
  }

  //! Calculate macroscopic velocities
  Vector<Double> CalcVelocities(unsigned int idx);

  // extract probability distributions for output
  Vector<Double> ExtractDistribution(unsigned int idx);

  Matrix<Double> couplingNodes_;

private:

  /** part of the constructor */
  void InitRegions(PtrParamNode pn, Grid* grid);

  /** setup elemen mapping
   * @see elem_to_idx
   * @see ids_to_elem */
  void SetupElementMapping(Grid* grid);

  //! Define available result types
  void DefinePostProcResults();

  //! Define available primary result types
  void DefinePrimaryResults();

  inline double GetPdf(unsigned int idx, int dir) const  {
    return pdfs_[idx * n_q_ + dir];
  };

  inline double& GetPdf(unsigned int idx, int dir) {
    return pdfs_.GetPointer()[idx * n_q_ + dir];
  };

  inline unsigned int GetIndex(unsigned int x, unsigned int y, unsigned int z) const {
    return z * n_x_ * n_y_ + y * n_x_ + x;
  }

  /** Perform position of matrix element in linearized matrix*/
  inline unsigned int GetMatrixElemId(unsigned int row, unsigned int col, unsigned int ncols) {
    return row * ncols + col;
  }

  inline void GetCoordinates(unsigned int index, unsigned int& x, unsigned int& y, unsigned int&z) const{
    x =  index % n_x_;
    index /= n_x_;
    y = index % n_y_;
    index /= n_y_;
    if (dim_ == 3)
      z = index;
  }

  inline unsigned int GetPdfIndex(unsigned int x, unsigned int y, unsigned int z, unsigned int dir) const {
    return GetIndex(x,y,z) * n_q_ + dir;
  }

  inline unsigned int GetPdfIndex(unsigned int index, unsigned int dir) const {
    return index * n_q_ + dir;
  }

  inline bool PointsToBoundary(unsigned int x, unsigned int y, unsigned int z, unsigned int dir)
  {
    LatticeBoltzmann::PDFDirectionVector tmp = microVelDirections_[dir];
    int tmp_x = x + tmp.off_x;
    int tmp_y = y + tmp.off_y;
    int tmp_z = z + tmp.off_z;
    return (tmp_x < 0 || tmp_x >= (int)n_x_ || tmp_y < 0 || tmp_y >= (int)n_y_ || tmp_z < 0 || tmp_z >= (int)n_z_) ;
  }

  // testing PointsToBoundary()
  void TestPointsToBoundary();

  // reads discrete velocities from extern LBM simulation
  void ReadProbabilityDistribution(const std::string& file);

  /** Setup structure for calling solver */
  void SetupElements();

  void SetupParabolicInflow();

  /** write the single CFS2LBM.dat dat (new interface) */
  void ExportCFS2LBM(const StdVector<double>& elements);

  /** write the files for the old interface: data.dat, obst.dat, por.dat, non_sing_dat. These are e.g. generated by thomas guess's matlab optimization */
  void ExportMultipleFiles(const StdVector<double>& elements);

  /** writes given matrix in file with given name*/
  void WriteMatrix(const std::string& file, const compressed_matrix<double> & M);

  /**  This method computes the indices of the adjoint system which avoid a singular Jacobian. Based on Georg Pingen and Thomas Guess, see non_singularities_new.m */
  void SetNonSingualrityIndices();

  // void DeleteSingularities(const compressed_matrix<double> & M,compressed_matrix<double> & output);
  void DeleteSingularities(const mapped_matrix<double> & M, compressed_matrix<double> & output);

  /** Sets up local data for SensitivityAnalysis() */
  void SetupSensitivityAnalysis(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz, StdVector<double>& dcol, StdVector<double>& weights);

  void d_collision_step_d_f(unsigned int index, Matrix<double>& block, StdVector<double>& dfeqdux, StdVector<double>& dfeqduy, StdVector<double>& dfeqduz, const StdVector<double>& ux, const StdVector<double>& uy, const StdVector<double>& uz, const StdVector<double>& dcol, const StdVector<double>& weight);

  void d_bounceback_d_f(int index, Matrix<double>& block);
  void d_inflow_d_f(int index, Matrix<double>& block, StdVector<double>& weight);
  void d_outflow_d_f(int index, Matrix<double>& block, StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz, StdVector<double>& dloc, StdVector<double>& weight);

  // void d_propagate_d_rho(compressed_matrix<double>& Jprop, const compressed_matrix<double> & J);
  void d_propagate_d_f(mapped_matrix<double>& Jprop, const mapped_matrix<double> & J);

  Vector<double> d_pressuredrop_d_f(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz);

  void matrix_sparse_to_crs(compressed_matrix<double>& M, double* a, unsigned int* ia, unsigned int* ja);
//  void matrix_sparse_to_crs(mapped_matrix<double>& M, double* a, unsigned int* ia, unsigned int* ja);

  //! Method of smoothing
  std::string method_;

  //! Flag indicating if PDE is assembled for first time
//  bool firstTurn_;

  //! Vector storing factors for adapted pseudo mechanic bulk modulus
  Vector<Double> factor_;

  /** this stores the boundary region id */
  RegionIdType boundary_reg_;
  /** this stores the obstacle region id */
  RegionIdType obstacle_reg_;

  /** this stores the design region */
  StdVector<RegionIdType> design_reg_;

  /** LBM requires the data do be ordered in a defined way (lexicographic). The LBM data covers all regions (design and boundary).
   * To allow access by region (CalcResults()) this does the mapping from CFS elemNum to lbm idx and lbm idx to elemNum. lbm idx is 0-based
   * elem_to_idx has the size of the largest elem_num and contains the idx or -1 where there is no mapping.
   * @see SetupElemMapping() */
  StdVector<int> elem_to_idx;

  /** size n_elems and contains the CFS elemNr
   * @see elem_to_idx */
  StdVector<int> idx_to_elem;

  /** number of discrete velocities: 9 for D2Q9 or n_q_ for D3Q19 */
  unsigned int n_q_;

  /** extents of computational grid */
  unsigned int n_x_, n_y_, n_z_;

  /** total number of elements */
  unsigned int n_elems;

  // stores how often CalcResults() was called. We need this for StoreResults() for last time step
  unsigned int numWriteResults_;

  // stores how many iterations were necessary to solve one flow problem with LBM
  unsigned int numIterations_;

  /** the complete structure, special elements (negative) and porosities */
  StdVector<double> elements;

  /** This are the indices of the later Jacobian matrix entries which are non-singular.
   * @see SetNonSingualrityIndices() */
  StdVector<unsigned int> non_sing;

  /** storage for particle distribution. This is the simulation result for the function evaluation. */
  StdVector<double> pdfs_;

  /** storage for adjoint particle distribution. This is the result of an adjoint SRT simulation. */
  StdVector<double> adjPdfs;

  /** these are the indices of the inlet elements */
  StdVector<unsigned int> inlet;

  /** these are the indices of the outlet elements */
  StdVector<unsigned int> outlet;

  StdVector<int> bb;  // indices of bounce back nodes
  StdVector<int> rel; // indices of the fluid m_nodes

  StdVector<unsigned int> obstacles;

  double omega_; /** molecular collision frequency */
  double Re_; /** Reynold's number of flow problem */
  double maxWallTime_;
  unsigned int maxIter_;
  double convergence_; /** value for convergence criterion */
  unsigned int writeFrequency_;

  // inlet velocities for constant inflow at all inlet nodes
  double u_max_x_;
  double u_max_y_;
  double u_max_z_;

  unsigned int dim_; // number of spatial dimensions of domain

  bool parabolicInflow_; // indicates if inflow velocity profile is parabolic

  // inlet velocities for parabolic inflow profile
  StdVector< StdVector<double> > u_in_;

  StdVector<LatticeBoltzmann::PDFDirectionVector> microVelDirections_;
  StdVector<LatticeBoltzmannBase::Direction> invPDFDirections_;

  StdVector<Matrix<double> > adjSRTCollision; // adjoint SRT collision matrices
  StdVector<StdVector<double> > adjCollisions;
  StdVector<StdVector<double> > d_pdrop_d_f;

  /** external lbm */
  std::string executable;

  /** internal lbm solver */
  LatticeBoltzmann* lbm;

  /** type of interface */
  Iface iface_;
  Enum<Iface> iface;

  /** Use Iface enums for indicating external adjoint solving or internal adjoint SRT LBM simulation */
  Iface adjSRT_;

  /** total time of state problems */
  Timer state_;

  /** total time of adjoint solution */
  shared_ptr<Timer> forwardSim_;
  shared_ptr<Timer> backwardSim_;
  shared_ptr<Timer> setupAdjoint_;
  shared_ptr<Timer> setupPardiso_;
  shared_ptr<Timer> solveLSE_;
  shared_ptr<Timer> dProp_;
  shared_ptr<Timer> dPressDrop_;
  shared_ptr<Timer> dColl_;
  shared_ptr<Timer> rhsSetup_;
  shared_ptr<Timer> delSing_;
};

} // end of namespace
#endif
