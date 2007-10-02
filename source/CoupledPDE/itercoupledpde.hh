// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ITERCOUPLEDPDE_2003
#define FILE_ITERCOUPLEDPDE_2003

#include "PDE/basePDE.hh"
#include "General/exception.hh"

namespace CoupledField
{

  // forward class declaration
  class IterSolveStep;
  class StdPDE;
  class SinglePDE;
  class PDECoupling;
  class ParamNode;

  //! This class iteratively solve a list of given SinglePDEs 
  class IterCoupledPDE : public BasePDE
  {

    // friend declaration
    friend class IterSolveStep;

  public:

    //! Constructor
    IterCoupledPDE(StdVector<StdPDE*> & PDEs,
                   StdVector<SinglePDE*> & sinlgePDEs,
                   StdVector<PDECoupling*> & Couplings,
                   ParamNode* paramNode); 

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

    //! Solve static step
    //void SolveStepStatic(const UInt kstep, const Double asteptime, 
    //           const bool updatesysmat);
  
    //! solve transient step
    //void SolveStepTrans(const UInt kstep, const Double asteptime,
    //                  const bool updatesysmat);
  
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

    //! Write information about coupling interfaces
    //! and coupling setup into ostream
    void WriteCouplingInfo(std::ostream &out);

    //! calculates the norm of a vector
    //Double CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval);

    UInt maxiter_;                        //!< maximum number of iterations per time step
    StdVector<Double> norms_;              //!< norm of coupling values

    //! pointer to SolveStep classes
    IterSolveStep * solveStep_;

    //! Parameter node for "iterative" coupling
    ParamNode * myParam_;

    //! Flag for nonlinear logging
    bool nonLinLogging_;
  
    // general PDE parameters
    AnalysisType analysistype_;         //!< type of analysis
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
