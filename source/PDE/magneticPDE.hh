// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MAGNETICPDE
#define FILE_MAGNETICPDE

#include <map>
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
    MagPDE( Grid * aptgrid, PtrParamNode paramNode );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagPDE();


    //! Initialize PDE
    virtual void Init( UInt sequenceStep, PtrParamNode base = PtrParamNode() );

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

    //! do jobs, before analysis starts
    void PreparePDE4Computation();

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

    //! computes the H field according to a hysteresis operator
    void CalcHfield( shared_ptr<BaseResult> result );

    //! computes the eddy current denstiy
    template<class TYPE>
    void CalcEddyCurrent( shared_ptr<BaseResult> result );

    //! computes the eddy current power
    template<class TYPE>
    void CalcEddyPower( shared_ptr<BaseResult> result );

    //! Calculate the total flux/flux derivative
    template<class TYPE>
    void CalcFlux( shared_ptr<Coil>, 
                   TYPE& flux, 
                   bool timeDeriv );

    //! computation of force (virtual work principle)
    void CalcForceVWP( shared_ptr<BaseResult> result );

    //! computation of Lorentz force
    template<class TYPE>
    void CalcForceLorentz( shared_ptr<BaseResult> result );
    
    //! computation of nodal Lorentz force
    template<class TYPE>
    void CalcNodeForceLorentz( Vector<TYPE> & force,
                               const StdVector<RegionIdType> regionIds,
                               const std::map<UInt, UInt>& cplNodeNumPos );
    
    // ---- Magnetic Force variables ---
 
    //! map coupling node number to its position
    StdVector<std::map<UInt, UInt> > cplNodeNumPos_;

    //! assigns each coupling element node the according Coupling Node number
    
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; 

    //! force operator (for coupling as well as postprocessing)
    MagForceOp* ForceOpVWP_;
    
    //! regions on which the Lorentz force should be calculated
    std::set<RegionIdType> regionsForceL_;
    
    //! vector containing regionIds of non-conforming interfaces
    StdVector<RegionIdType> ncIFaces_;
    
    // =======================================================================
    //   HELPER METHODS FOR CALCULATING AUXILIARY QUANTITIES 
    // =======================================================================
    
    //! Calc FluxDensity in integration point

    //! Calculates the flux density B = rot(A) at the given integration point.
    //! If ip is 0, the midpoint of the element is evaluated.
    template<class TYPE>
    void CalcFluxDensityAtIP( EntityIterator it,
                              UInt ip,
                              Vector<TYPE>& field );

    //! Calculates the magnetic field H at the given integration point.
    //! If ip is 0, the midpoint of the element is evaluated.
    void CalcHfieldAtIP( EntityIterator it, UInt ip,
                         Vector<Double>& field );

    //! Calc EddyCurrent in integration point

    //! Calculates the eddy current (density) at the given integration point.
    //! If ip is 0, the midpoint of the element is evaluated.
    template<class TYPE>
    void CalcEddyCurrentAtIP( EntityIterator it,
                              UInt ip,
                              Vector<TYPE>& field );
    

    // =======================================================================
    //   COILS
    // =======================================================================
    
    //@{ \name Attributes related to coils
    
    //! Names of coils resp. their subdomains
    StdVector<RegionIdType> coilRegionId_;  
    
    //! Parameters of the individual coils;
    StdVector<shared_ptr<Coil> > coilDef_;
    
    //! Write the induced voltage to file
    template<class TYPE>
    void WriteUI2File( );
    
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
    
    
  private:
    //!true, if 3d computation
    bool is3d_;
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

