// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SMOOTHPDE
#define FILE_SMOOTHPDE

#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "SinglePDE.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField {
class BaseResult;
class Grid;
class PDECoupling;
}  // namespace CoupledField
 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  class SmoothPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    SmoothPDE( Grid* Grid, PtrParamNode paramNode );

    //!  Deconstructor
    virtual ~SmoothPDE() {;};

    void InitStabParams();

    //! Initialize NonLinearities
    virtual void InitNonLin();


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
  
    virtual void ReadSpecialResults();
    
    //! calculate coupling terms
    virtual void CalcOutputCoupling();

    //! perform postprocessing on data
    void CalcResults( shared_ptr<BaseResult> result );
  
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output);

    //! Contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;

    
    Matrix<Double> elastFactors_;
    Matrix<Double> couplingNodes_;

    
    UInt stressDim_;
    
  private:

    //! Define available result types
    void DefineAvailResults();

    //! Method of smoothing
    std::string method_;

    std::string elastWeight_;
    
    //! Flag indicatingm if PDE is assembled for first time
    bool firstTurn_;

    //! Vector storing factors for adapted pseudo mechanic bulk modulus
    Vector<Double> factor_;

    Double characteristicLength_;
    Double exponent_;	
  };
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SmoothPDE
  //! 
  //! \purpose 
  //! This class implements a simple, quasi-mechanic grid smoother.
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
