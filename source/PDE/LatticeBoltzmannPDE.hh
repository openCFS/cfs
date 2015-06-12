// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_EXTLBMPDE
#define FILE_EXTLBMPDE

#include <string>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "PDE/NonFEM/LatticeBoltzmann.hh"
#include "SinglePDE.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{
class BaseResult;
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
  typedef enum { INTERNAL, EXT_MATLAB, EXT_CFSxLBM } Iface;

  //!  Constructor. here we read integration parameters
  /*!
      \param aGrid pointer to grid
   */
  LatticeBoltzmannPDE(Grid* grid, PtrParamNode pn);

  ~LatticeBoltzmannPDE();

  //    //!  Deconstructor
  //    virtual ~ExtLBMPDE() {;};

  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators();

  //! define the SoltionStep-Driver
  virtual void DefineSolveStep();

  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);

  //! initialize time stepping: nothing to do in smoother!
  void InitTimeStepping();

  //! set time step
  //! \param dt Current time step
  virtual void SetTimeStep(const Double dt){};

  //! calculate coupling terms
  virtual void CalcOutputCoupling();

  /** actually calls LBM. */
  void Solve();

  /** Sensitivity analysis (gradient) necessary for the optimization.
    Performs the missing propagation step for the sensitivity analysis. */
  void SensitivityAnalysis(TransferFunction* tf, Function* f, DesignSpace* space);

  /** For the olt multi-file interface, reads the files x_pardiso2.dat and dRdp.mtx and computes the gradient */
  void SetPrecalculatedGradient(StdVector<DesignElement*>& design, Function* f);

  /** Perform postprocessing on data. */
  void CalcResults( shared_ptr<BaseResult> result );

  /** implementation of objective function */
  double CalcPressureDrop();

  /** Another objective function: Mass flow rate*/
  double CalcFlowRate();

  //! returns if PDE can compute the quantity
//  virtual bool HasOutput(SolutionType output);

  // returns how often CalcResults() was called. We need this to write out the last simulation step
  int GetNumWriteResults();

  virtual void ReadSpecialResults();

  std::string ToString(const StdVector<double>& elements, bool x_fast, bool as_int) const;

  /** exports a boost compressed matrix in Matrix-Market format */
  // static void ToFile(const std::string& file, const compressed_matrix<double>& M);
  static void ToFile(const std::string& file, const mapped_matrix<double>& M);

  void create_output(const char * file);

  Iface GetIface() const { return iface_; }

  //! Contains LBM velocity
  NodeStoreSol<Double> solDeriv1_;

  Matrix<Double> couplingNodes_;

private:

  /** part of the constructor */
  void InitRegions(PtrParamNode pn, Grid* grid);

  /** setup elemen mapping
   * @see elem_to_idx
   * @see ids_to_elem */
  void SetupElementMapping(Grid* grid);

  //! Define available result types
  void DefineAvailResults();

  /** Calculate the LBM Density of an element idx */
  inline double CalcLBMDensity(unsigned int idx) const;

  inline double CalcVelocityX(unsigned int idx, double density) const;

  inline double CalcVelocityY(unsigned int idx, double density) const;

  inline double CalcVelocityZ(unsigned int idx, double density) const;

  inline double CalcPressure(unsigned int idx) const;

  /** Calculate volumetric flow rate (u*A) through surface of element with index idx. pdfDir is the pdf direction in which the outflow surface lies*/
  inline double CalcVolFlowRate(unsigned int idx, unsigned int pdfDir) const;

  inline double GetPdf(unsigned int idx, int dir) const  {
    return pdfs[idx * n_q_ + dir];
  };

  inline double& GetPdf(unsigned int idx, int dir) {
    return pdfs.GetPointer()[idx * n_q_ + dir];
  };

  inline unsigned int GetIndex(unsigned int x, unsigned int y, unsigned int z ) const {
    return z * n_x_ * n_y_ + y * n_x_ + x;
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
    LatticeBoltzmann::PDFDirectionVector tmp = (*microVelDirections_)[dir];
    int tmp_x = x + tmp.off_x;
    int tmp_y = y + tmp.off_y;
    int tmp_z = z + tmp.off_z;
    return (tmp_x < 0 || tmp_x >= (int)n_x_ || tmp_y < 0 || tmp_y >= (int)n_y_ || tmp_z < 0 || tmp_z >= (int)n_z_) ;
  }

  // testing PointsToBoundary()
  void TestPointsToBoundary();

  /** Calculate if given pdf direction points to  an element's surface */
  inline bool PointsToSurface(unsigned int pdfDir)
  {
    // A pdf points to a surface if its direction vector is 1 --> unit vector
    return fabs((*microVelDirections_)[pdfDir].off_x + (*microVelDirections_)[pdfDir].off_y + (*microVelDirections_)[pdfDir].off_z) == 1;
  }

  //! Calculate macroscopic velocities
  void CalcVelocities(shared_ptr<BaseResult> res);

  //! Calculate densities
  void CalcDensities(shared_ptr<BaseResult> res);

  void CalcPressures(shared_ptr<BaseResult> res);

  // reads discrete velocities from extern LBM simulation
  void ReadProbabilityDistribution(const std::string& file);

  // extract probability distributions for output
  void ExtractDistribution(shared_ptr<BaseResult> base_result);

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

  /** Compute and store PDF direction of an outlet element pointing to a outlet surface. We need this information to calculate the flow rate. */
  void SetOutlfowPDFs();

  // void DeleteSingularities(const compressed_matrix<double> & M,compressed_matrix<double> & output);
  void DeleteSingularities(const mapped_matrix<double> & M, compressed_matrix<double> & output);

  /** Sets up local data for SensitivityAnalysis() */
  void SetupSensitivityAnalysis(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz, StdVector<double>& dcol, StdVector<double>& weights);

  void d_collision_step_d_f(unsigned int index, Matrix<double>& block, StdVector<double>& dfeqdux, StdVector<double>& dfeqduy, StdVector<double>& dfeqduz, const StdVector<double>& ux, const StdVector<double>& uy, const StdVector<double>& uz, const StdVector<double>& dcol, const StdVector<double>& weight);

  void d_bounceback_d_f(Matrix<double>& block);
  void d_inflow_d_f(int index, Matrix<double>& block, StdVector<double>& weight);
  void d_outflow_d_f(int index, Matrix<double>& block, StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz, StdVector<double>& dloc, StdVector<double>& weight);

  // void d_propagate_d_rho(compressed_matrix<double>& Jprop, const compressed_matrix<double> & J);
  void d_propagate_d_f(mapped_matrix<double>& Jprop, const mapped_matrix<double> & J);

//  void d_propagate_d_f_inDir(mapped_matrix<double>& Jprop, const mapped_matrix<double> & J,int dir);

  Vector<double> d_pressuredrop_d_f(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz);

  //! derivative of flow rate objective function
  Vector<double> d_flowrate_d_f(StdVector<double>& ux, StdVector<double>& uy, StdVector<double>& uz);

  void matrix_sparse_to_crs(compressed_matrix<double>& M, double* a, unsigned int* ia, unsigned int* ja);
  void matrix_sparse_to_crs(mapped_matrix<double>& M, double* a, unsigned int* ia, unsigned int* ja);

  //! Method of smoothing
  std::string method_;

  //! Flag indicating if PDE is assembled for first time
  bool firstTurn_;

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

  /** number of discrete velocities: 9 for D2Q9 or 19 for D3Q19 */
  unsigned int n_q_;

  /** extents of computational grid */
  unsigned int n_x_, n_y_, n_z_;

  /** total number of elements */
  unsigned int n_elems;

  // stores how often CalcResults() was called. We need this for StoreResults() for last time step
  unsigned int numWriteResults_;

  /** the complete structure, special elements (negative) and densities */
  StdVector<double> elements;

  /** This are the indices of the later Jacobian matrix entries which are non-singular.
   * @see SetNonSingualrityIndices() */
  StdVector<unsigned int> non_sing;

  /** storage for for particle distribution. This is the simulation result for the function evaluation. */
  StdVector<double> pdfs;

  /** these are the indices of the inlet elements */
  StdVector<unsigned int> inlet;

  /** these are the indices of the outlet elements */
  StdVector<unsigned int> outlet;

  StdVector<unsigned int> obstacles;

//  /** these are the indices of the elements og the region 'inclusion'*/
//  StdVector<unsigned int> obstacle;

  double omega_; /** molecular collision frequency */
  double Re_; /** Reynold's number of flow problem */
  double maxWallTime_;
  unsigned int maxIter_;
  double convergence_; /** value for convergence criterion */
  unsigned int writeFrequency_;

  int lbmCalls_; // counts how often LBM solver was called

  // inlet velocities for constant inflow at all inlet nodes
  double u_max_x_;
  double u_max_y_;
  double u_max_z_;

  unsigned int dim_; // number of spatial dimensions of domain

  bool parabolicInflow; // indicates if inflow velocity profile is parabolic

  // inlet velocities for parabolic inflow profile
  StdVector< StdVector<double> > u_in_;

  StdVector<LatticeBoltzmann::PDFDirectionVector>* microVelDirections_;
  StdVector<LatticeBoltzmannBase::Direction>* invPDFDirections_;
  /** Storage for PDF directions of outlet elements pointing to a boundary surface. We need this information to calculate the flow rate. */
  StdVector<LatticeBoltzmannBase::Direction>* PDFToBoundarySurfaces;

  /** external lbm */
  std::string executable;

  /** internal lbm solver */
  LatticeBoltzmann* lbm;

  /** type of interface */
  Iface iface_;
  Enum<Iface> iface;

  /** total time of state problems */
  Timer state_;

  /** total time of adjoint solution */
  Timer adjoint_;

};

} // end of namespace
#endif
