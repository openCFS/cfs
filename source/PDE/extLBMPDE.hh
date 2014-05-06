// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_EXTLBMPDE
#define FILE_EXTLBMPDE

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "SinglePDE.hh"
#include "Utils/nodestoresol.hh"

// D2Q9
#define _QN_  9
// Address PDFs via Coordinates and Direction.
#define PDF(_x, _y, _direction)  pdfs[ ((n_x) * (n_y)) * (_direction) +  (_y) * n_x + (_x) ]
// Address PDFs via Index and Direction only.
#define PDF_IDX(_index, _direction)  pdfs[(_index) * _QN_ + (_direction)]

namespace CoupledField {
class BaseResult;
class Grid;
class PDECoupling;
}  // namespace CoupledField
 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  class ExtLBMPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    ExtLBMPDE(Grid* grid, PtrParamNode pn);

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

    //! perform postprocessing on data
    void CalcResults( shared_ptr<BaseResult> result );
  
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output);

    virtual void ReadSpecialResults();

    // reads discrete velocities from extern LBM simulation
    void readData(const char * file);

    // writes data file for extern LBM solver
    virtual void writeData(const char * file);

    //! Contains LBM velocity
    NodeStoreSol<Double> solDeriv1_;
    
    Matrix<Double> couplingNodes_;

  private:

    //! Define available result types
    void DefineAvailResults();

    //! Calculate macroscopic velocities
    void CalcVelocities(shared_ptr<BaseResult> res);

    //! Calculate densities
    void CalcDensities(shared_ptr<BaseResult> res);

    //! Method of smoothing
    std::string method_;

    //! Flag indicating if PDE is assembled for first time
    bool firstTurn_;

    //! Vector storing factors for adapted pseudo mechanic bulk modulus
    Vector<Double> factor_;

    /** this stores the boundary region id */
    RegionIdType boundary_;

    //extents of computational grid
    UInt n_x,n_y,n_z;
    //total number of elements
    UInt n_elems;
    //Storage for PDFs
    Double * pdfs;

    std::string executable;
    std::string sim_type;

  };
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class ExtLBMPDE
  //! 
  //! \purpose 
  //! This class implements a pseudo PDE for integrating an extern LBM simulation
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace
#endif
