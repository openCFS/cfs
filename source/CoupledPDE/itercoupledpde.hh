// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ITERCOUPLEDPDE_2003
#define FILE_ITERCOUPLEDPDE_2003

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/exception.hh"
#include "PDE/basePDE.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class Assemble;
class BaseSolveStep;
}  // namespace CoupledField

namespace CoupledField
{

  // forward class declaration
  class IterSolveStep;
  class PDECoupling;
  class SinglePDE;
  class StdPDE;

  //! This class iteratively solves a list of given SinglePDEs 
  class IterCoupledPDE : public BasePDE
  {

    // friend declaration
    friend class IterSolveStep;

  public:

    //! Constructor
    IterCoupledPDE(StdVector<StdPDE*> & PDEs,
                   StdVector<SinglePDE*> & sinlgePDEs,
                   StdVector<PDECoupling*> & Couplings,
                   PtrParamNode paramNode); 

    //! Destructor
    ~IterCoupledPDE();

    //! calculates coupling interfaces
    void InitCoupling();
  
    //! write general defines (BCs, loads, etc.) to info-file
    void WriteGeneralPDEdefines();

    Assemble * getPDE_assemble() {
      EXCEPTION( "Get Assemble-Object makes no sense for itercoupledPDE");
    };

    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();

    //! write a restart file "simname_pdename.restart"
    void WriteRestart( );


    //! read a restart file "simname_pdename.restart"
    void ReadRestart(UInt &startStep );

    //! write results in file
    void WriteResultsInFile(const UInt kstep = 0,
                            const Double asteptime = 0.0 );

    //! set time step
    //! \param dt Current time step
    void SetTimeStep(const Double dt);

  protected:

    /** Write coupling info. TODO -> check for overloading! */ 
    void ToInfo(PtrParamNode in);

    //! calculates the norm of a vector
    //Double CalcNorm(NormType normtype, SingleVector & val, SingleVector & oldval);

    UInt miniter_;                        //!< minimum number of iterations per time step
    UInt maxiter_;                        //!< maximum number of iterations per time step
    StdVector<Double> norms_;              //!< norm of coupling values

    //! pointer to SolveStep classes
    IterSolveStep * solveStep_;

    //! Parameter node for "iterative" coupling
    PtrParamNode myParam_;

    //! Flag for nonlinear logging
    bool nonLinLogging_;
  
    // general PDE parameters
    StdVector<StdPDE *> PDEs_;         //!< list of belonging PDEs
    StdVector<SinglePDE*> singlePDEs_;
    StdVector<PDECoupling*> Couplings_; //!< vector of coupling objects
    UInt NumPDEs_;                   //!< number of PDEs 
  

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class IterCoupledPDE
  //! 
  //! \purpose 
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
