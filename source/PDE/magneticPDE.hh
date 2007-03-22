// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MAGNETICPDE
#define FILE_MAGNETICPDE

#include "SinglePDE.hh" 
#include "Forms/magforceop.hh"

namespace CoupledField
{

  // Forward declaration of classes
  class Coil;

  //! Class for defining the magnetic field in 2D
  class MagPDE : public SinglePDE
  {
  public:

    //! Constructor
    MagPDE( Grid * aptgrid, ParamNode* paramNode );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagPDE();


    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs();

    //! Read special store results
    void ReadSpecialResults();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();


    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! do PostProcessing step
    void CalcResults( shared_ptr<BaseResult> result );

    // ======================================================
    // COUPLING SECTION
    // ======================================================


    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);
  
    //! calculate coupling terms
    void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);

  protected:

    // =======================================================================
    //  Initialization
    // =======================================================================

    //! Define availabe result types
    void DefineAvailResults();

    //! Query parameter object for information on coils
    void ReadCoils();
    
    //! Query parameter object for information on permanent magnets
    void ReadMagnets();

    //! Init the time stepping
    void InitTimeStepping();


    // =======================================================================
    //  POSTPROCESSING
    // =======================================================================
    void WriteL2File(Vector<Double>& uiSD);

    //! computes the magnetic energy
    template<class TYPE>
    void CalcEnergy( shared_ptr<BaseResult> result );
    
    //! computes the magnetic flux density
    template<class TYPE>
    void CalcFluxDensity( shared_ptr<BaseResult> result );

    //! computes the eddy current denstiy
    template<class TYPE>
    void CalcEddyCurrent( shared_ptr<BaseResult> result );

    //! Calculate the total flux/flux derivative
    template<class TYPE>
    void CalcFlux( shared_ptr<BaseResult> result, 
                   bool timeDeriv );

    //! computation of force (virtual work principle)
    void CalcForceVWP( shared_ptr<BaseResult> result );

    //! computation of Lorentz force
    void CalcNodeForceLorentz(Vector<Double> & force, 
                              StdVector<StdVector<UInt> > & 
                              elemNodeToCouplingNode,
                              UInt actCoupling, 
                              UInt numCouplingNodes);
 
 
    // ---- Magnetic Force variables ---
 

    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; 

    //! force operator (for coupling as well as postprocessing)
    MagForceOp* ForceOpVWP_;

    // =======================================================================
    //   COILS
    // =======================================================================
    
    //@{ \name Attributes related to coils
    
    //! Names of coils resp. their subdomains
    StdVector<RegionIdType> coilRegionId_;  
    
    //! Parameters of the individual coils;
    StdVector<shared_ptr<Coil> > coilDef_;
    
    //! Write the induced voltage to file
    void WriteUI2File(Vector<Double>& uiSD);

    //! file for informational output of coils
    std::ofstream * UIfile_; 
    
    //! name of file for saving current/voltage values
    std::string UIfilename_;
    
    //@}

    // =======================================================================
    //   PERMANENT MAGNETS
    // =======================================================================
    
    //@{ \name Attributes related to permanent magnets
    
    //! Subdomains containing permanent magnets
    StdVector <RegionIdType> magnetsDomain_;
    
    //! x-component of direction of magnetisation for each magnet
    
    //! x-component of direction of magnetisation for each magnet
    //! \todo As suggested by Fred Hofer, the direction of magnetisation of a
    //! permanent magnet must now be specified in the XML parameter file and
    //! no longer in the material data file. While magneticPDE already reads
    //! these data, they are not yet used in the simulation.
    StdVector<Double> magnetsOriX_;
    
    //! y-component of direction of magnetisation for each magnet
    StdVector<Double> magnetsOriY_;

    //! z-component of direction of magnetisation for each magnet
    StdVector<Double> magnetsOriZ_;

    //@}

    // =======================================================================
    //   COILS
    // =======================================================================
    std::string nonLinMethod_;
    
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagPDE
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
  //! \improve Encapsulate the defintiion of magnets within an own struct
  //! 

#endif

} // end of namespace
#endif

