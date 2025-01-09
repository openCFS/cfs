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
  
  class EntityList;
  class CoefFunctionAccumulator;
  class DirectCoupledPDE;
  
  //! This class iteratively solve a list of given SinglePDEs 
  class IterCoupledPDE : public BasePDE
  {
    // friend declaration
    friend class IterSolveStep;

  public:


    
    //! Constructor
    IterCoupledPDE( StdVector<SinglePDE*>& singlePDEs,
                    StdVector<DirectCoupledPDE*>& cplPde,
                    PtrParamNode paramNode,
                    PtrParamNode infoNode,
                    shared_ptr<SimState> simState, Domain* domain ); 

    //! Destructor
    ~IterCoupledPDE();
    
    //! Retrieve info pointer
    PtrParamNode GetInfoNode() {
      return infoNode_;
    }

    //! Write general defines (BCs, loads, etc.) to info-file
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

    //! Triggers the finalization within the iterSolveStep
    void TriggerFinalize();

    //! Function similar to GetCouplingCoefFct but only returns updatedGeo
    void GetUpdateGeoForPDE( SolutionType type,
                             shared_ptr<EntityList>  list,
                             const std::string& pdeName,
                             bool& updateGeo );
    
    //! Update PDE due to updated step in multistep solution strategy
    virtual void UpdateToSolStrategy();
    
    //! write results in file
    void WriteResultsInFile(const UInt kstep = 0,
                            const Double asteptime = 0.0 );

    //! set time step
    //! \param dt Current time step
    void SetTimeStep(const Double dt);

  protected:

    //! Dump information to InfoNode 
    void ToInfo(PtrParamNode in);

    //! pointer to SolveStep classes
    IterSolveStep * solveStep_;

    //! Parameter node for "iterative" coupling
    PtrParamNode myParam_;

    //! Type of analysis
    BasePDE::AnalysisType analysistype_;
    
    //! Pointer to pdes
    StdVector<StdPDE *> PDEs_;
    
    //! Pointer to SinglePDEs
    StdVector<SinglePDE*> singlePDEs_;
    
    //! Pointer to IterativeCoupledPDEs
    StdVector<DirectCoupledPDE*> coupledPDEs_;

    //! Map storing the coeffunction accumulators
    std::map<SolutionType, shared_ptr<CoefFunctionAccumulator> > accu_;
    
    //! Number of PDEs
    UInt numPDEs_; 
    
    //! Info node for logging 
    PtrParamNode infoNode_; 

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
