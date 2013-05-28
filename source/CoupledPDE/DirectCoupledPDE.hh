// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DIRECT_COUPLED_PDE
#define FILE_DIRECT_COUPLED_PDE

#include "PDE/StdPDE.hh"

namespace CoupledField {

  // forward class declaration
  class SinglePDE;
  class BasePairCoupling;
  class Assemble;
  class AlgebraicSys;
  class StdPDE;
  
  //! This class implements the direct coupling of StdPDEs.
  class DirectCoupledPDE : public StdPDE
  {
  public:

    typedef StdVector<shared_ptr<BaseFeFunction> > FeFunctionList;
    typedef StdVector<shared_ptr<FeSpace> > FeSpaceList;
  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
    */
    DirectCoupledPDE( Grid *aptgrid, PtrParamNode paramNode,
                      PtrParamNode infoNode,
                      shared_ptr<SimState> simState, Domain* domain ); 
    
    //! Destructor
    virtual ~DirectCoupledPDE();
  
    //! Set SinglePDEs
    void SetSinglePDEs( const StdVector<SinglePDE*> &pdes);

    //! returns all single PDEs
    StdVector<SinglePDE*>& GetSinglePDEs()
    { return singlePDEs_;};

    //! Set Coupling objects
    void SetCouplings( const StdVector<BasePairCoupling*> &couplings);

    //! Get couplings object
    StdVector<BasePairCoupling*>* GetCouplingsObject() 
    { return &couplings_;};

    //! Initialization routine
    void Init(UInt sequenceStep );
            
    //! write general defines (BCs, loads, etc.) to info-file
    void WriteGeneralPDEdefines();
  

    //! set boundary condition
    void SetBCs();

    //!
    virtual void InitTimeStepping();
    

    //! Update PDE due to updated step in multistep solution strategy
    virtual void UpdateToSolStrategy();

 
    // ======================================================
    // POSTPROC SECTION
    // ======================================================
  
    //@{
    //! \name Methods performing post-processing
  
    //! perform cleanup and do last computations
    void Finalize();

    //! write results in file
    void WriteResultsInFile(const UInt kstep = 0,
                            const Double asteptime = 0.0 );

    //@}

  private:

    //! define the SolutionStep-Driver
    virtual void DefineSolveStep();

    //! Reads ncInterfaces defined in the XML file
    virtual void ReadNcInterfaces()
    { EXCEPTION("Not implemented"); };

    //! References to SinglePDEs
    StdVector<SinglePDE*> singlePDEs_;

    //! References to pairwise coupling objects
    StdVector<BasePairCoupling*> couplings_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class DirectCoupledPDE
  //! 
  //! \purpose This class is implements the concept of directly coupling two
  //! or more SinglePDE's together by the use of BasePairCoupling 
  //! objects. It has the same interface as a StdPDE, so the solveStep-classes
  //! can work transparently with it.
  //! 
  //! \collab This class is instantiated in the class Domain. It gets one or 
  //! more instances of SinglePDE (which get assembled on a main diagonal 
  //! position in the resulting algebraic system) and the according instances
  //! of BasePairCoupling, which implement the pairwise coupling between each 
  //! two SinglePDEs.
  //! 
  //! \implement
  //! 
  //! \status In use
  //! 
  //! \unused
  //! 
  //! \improve 
  //! - At the moment only a Newmark timestepping scheme is applied. In some
  //! cases this will not be sufficient.
  //! - Harmonic coupling has not been tested yet, as well as nonlinear 
  //! analysis
  //!

#endif

} // end of namespace

#endif
