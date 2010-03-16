// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MAGNETICEDGEPDE
#define FILE_MAGNETICEDGEPDE

#include <map>
#include "SinglePDE.hh" 
#include "Forms/magforceop.hh"


namespace CoupledField
{

  // Forward declaration of classes
  class Coil;

  //! Class for 3D magnetics using edge elements
  class MagEdgePDE : public SinglePDE
  {
  public:

    //! Constructor
    MagEdgePDE( Grid * aptgrid, PtrParamNode paramNode );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgePDE();


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
    //! computes the magnetic flux density
    template<class TYPE>
    void CalcFluxDensity( shared_ptr<BaseResult> result );

    //! computes the eddy current denstiy
    template<class TYPE>
    void CalcEddyCurrent( shared_ptr<BaseResult> result );

    //! computes the magnetic energy
    template<class TYPE>
    void CalcEnergy( shared_ptr<BaseResult> result );
    
    //! computes the divergence of the magnetic potential
    template<class TYPE>
    void CalcMagPotentialDiv( shared_ptr<BaseResult> result );
    
    // ---- Magnetic Force variables ---
 
    //! map coupling node number to its position
    StdVector<std::map<UInt, UInt> > cplNodeNumPos_;

    //! assigns each coupling element node the according Coupling Node number
    
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; 

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

    //! Calc EddyCurrent in integration point

    //! Calculates the eddy current (density) at the given integration point.
    //! If ip is 0, the midpoint of the element is evaluated.
    template<class TYPE>
    void CalcEddyCurrentAtIP( EntityIterator it,
                              UInt ip,
                              Vector<TYPE>& field );
                              
    //! computation of Lorentz force
    void CalcNodeForceLorentz( Vector<Double> & force, 
                               std::map<UInt, UInt>& cplNodeNumPos,
                               UInt actCoupling, 
                               UInt numCouplingNodes );
    

    // =======================================================================
    //   COILS
    // =======================================================================
    
    //@{ \name Attributes related to coils
    
    //! Names of coils resp. their subdomains
    StdVector<RegionIdType> coilRegionId_;  
    
    //! Parameters of the individual coils;
    StdVector<shared_ptr<Coil> > coilDef_;
    
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
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgePDE
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

