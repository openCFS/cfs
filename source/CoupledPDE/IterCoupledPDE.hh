// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ITERCOUPLEDPDE_2003
#define FILE_ITERCOUPLEDPDE_2003

#include "PDE/BasePDE.hh"
#include "General/Exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

  // forward class declaration
  class IterSolveStep;
  class StdPDE;
  class SinglePDE;
  class ParamNode;
  class EntityList;

  //! This class iteratively solve a list of given SinglePDEs 
  class IterCoupledPDE : public BasePDE
  {

    // friend declaration
    friend class IterSolveStep;

  public:

//    
//    //class for storing the convergence criterions
//    class CplQuantity {
//    public:
//      
//      //! Constructor with accumulator class
//      CplQuantity(SolutionType solType);
//      virtual ~CplQuantity();
//      
//      void AddCoefFct( PtrCoefFct coef, shared_ptr<EntityList> list );
//      
//      
//      
//    protected:
//
//      SolutionType solType_;
//      
//      StdVector<PtrCoefFct> 
//    };
    
    
    //! Constructor
    IterCoupledPDE(StdVector<StdPDE*> & PDEs,
                   PtrParamNode paramNode); 

    //! Destructor
    ~IterCoupledPDE();

    //! calculates coupling interfaces
    void InitCoupling();
  
    //! write general defines (BCs, loads, etc.) to info-file
    void WriteGeneralPDEdefines();

    Assemble * GetAssemble() {
      EXCEPTION( "Get Assemble-Object makes no sense for itercoupledPDE");
    };

    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();

    //! Obtain coupling quantity
    
    //! This method returns a given coefficient function
    //! from a contained SinglePDE. Internally, it creates an additional
    //! surrounding struct to calculate some norm for determining an
    //! stopping criterion when evaluating the CoefFunction.
    PtrCoefFct GetCouplingCoefFct( SolutionType type,
                                   shared_ptr<EntityList>  list,
                                   const std::string& pdeName,
                                   bool& updatedGeo );
    
    //! Update PDE due to updated step in multistep solution strategy
    virtual void UpdateToSolStrategy();
    
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

    //! Dump information to InfoNode 
    void ToInfo(PtrParamNode in);

    //! Maximum number of iterations per time step
    UInt maxiter_;                        

    //! pointer to SolveStep classes
    IterSolveStep * solveStep_;

    //! Parameter node for "iterative" coupling
    PtrParamNode myParam_;

    //! Flag for nonlinear logging
    bool nonLinLogging_;
  
    //! Type of analysis
    BasePDE::AnalysisType analysistype_;
    
    //! Pointer to pdes
    StdVector<StdPDE *> PDEs_;
    
    //! Pointer to SinglePDEs
    StdVector<SinglePDE*> singlePDEs_;
    
    //! Number of PDEs
    UInt numPDEs_; 
    
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
