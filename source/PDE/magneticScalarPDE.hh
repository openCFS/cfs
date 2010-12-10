// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MAGNETIC_SCALARPDE_HH
#define FILE_MAGNETIC_SCALARPDE_HH

#include "SinglePDE.hh" 

namespace CoupledField
{

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  
  //! Magnetic Field PDE using the total magnetic scalar potential
  
  //! This class implements the magnetostatic field problem by utilizing
  //! the total scalar magnetic potential. The class does not consider
  //! any time derivatives, i.e. eddy currents. The excitation of the magnetic
  //! field can either be incorporated using Dirichlet boundary conditions
  //! of by defining e.g. piecewise current sticks, using the law of
  //! Biot-Savart.
  class MagScalarPDE : public SinglePDE {

  public:
    //! Constructor
    /*!
      \param 
      \param aGrid pointer to grid
      \param aGrid pointer to class Grid
    */
    MagScalarPDE( Grid* aptgrid, PtrParamNode paramNode );

    //! Destructor
    virtual ~MagScalarPDE(){};
        
    //! Initialize NonLinearities
    void InitNonLin();    

    //! Define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: will be implemented in the next step
    void InitTimeStepping() {;};

    //! Nothing to do
    void SetTimeStep(const Double dt) {;};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> result );
    
    // ======================================================
    // COUPLING SECTION
    // ======================================================
    
    //! Turn the magnetostriction coupling on

    //! Triggers the correct assembly of the magnetic block in a 
    //! magnetostrictive-coupled simulation, because the coupled magnetic block
    //! is negative compared to the normal one
    void SetMagStrictCoupling();


    //! Initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! Calculate coupling terms
    void CalcOutputCoupling();

    //! Returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);
  
  protected:

    // *****************
    //  POSTPROCESSING
    // *****************

    //! Define availabe result types
    void DefineAvailResults();

    //! Calculate magnetic  field intensity (H) 
    template <class TYPE>
    void CalcMagneticField( shared_ptr<BaseResult> vals );
  
//    //! Calculate H_s 
//    template <class TYPE>
//    void CalcHs( shared_ptr<BaseResult> vals );
//    
//    //! Calculate H_m 
//    template <class TYPE>
//    void CalcHm( shared_ptr<BaseResult> vals );
        
    //! Calculate magnetic  flux density (B)
    template <class TYPE>
    void CalcMagneticFluxDensity( shared_ptr<BaseResult> vals );
       
//    //! Computes the magnetic energy for each subdomain
//    template <class TYPE>
//    void CalcEnergy( shared_ptr<BaseResult> vals );

    
  private:

    //! Initialize the defined and available results types
    void InitStoreResults();
    
    //!true, if 3d computation
    bool is3d_;
    
    //! flag for Magnetostriction-coupling
    bool isMagnetostrictiveCoupled_; 

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class magneticscalarPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. 
  
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


