// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_MAGNETIC_PDE_HH
#define FILE_CFS_MAGNETIC_PDE_HH

#include <map>

#include "SinglePDE.hh"

namespace CoupledField
{


  //! Class for magnetic field calculation (eddy current case)
  
  //! This class is utilized to solve the magnetic problem (eddy current case)
  //! in 2D and 3D using nodal finite elements. The primary unknown is the 
  //! magnetic vector potential \[\vec{A}\]. 
  //! In case of 3D time dependent / harmonic simulations an additional scalar
  //! potential \[\Phi\] is utilized to ensure that the solution is embedded
  //! in the space HCurl.
  class MagneticPDE: public SinglePDE {

  public:

    //! Constructor
    MagneticPDE( Grid *aGrid, PtrParamNode paramNode );

    //!  Destructor
    virtual ~MagneticPDE();

    //! Initialize NonLinearities
    void InitNonLin();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs();
    
    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){};

    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();
    
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Read special results definition
    void ReadSpecialResults();

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! calculate coupling terms
    void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> results );
    
    
    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );
  protected:

    //! Flag for mixed formulation
    
    //! In case of a transient / harmonic 3D simulation, we need 
    //! an additional scalar potential to ensure HCurl compatibility.
    bool isMixed_;

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
    
    
    //! Query parameter object for information on coils
    void ReadCoils();

    //! Query parameter object for information on permanent magnets
    void ReadMagnets();
    
    //! Initialize time stepping method
    void InitTimeStepping();
    
    //! read in softening types
    void ReadSoftening();

    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

  };

} // end of namespace
#endif

