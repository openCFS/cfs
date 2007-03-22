// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOUSTICXYZPDE_2005
#define FILE_ACOUSTICXYZPDE_2005

#include "SinglePDE.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODESolve/ODESolver_RKF45.hh"
 
namespace CoupledField {

  //! Class for acoustic equation (no adaptivity)


  class AcousticXYZPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    AcousticXYZPDE(Grid *aGrid, ParamNode* paramNode);

    //! Destructor
    virtual ~AcousticXYZPDE(){};

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! perform postprocessing on data
    void CalcResults( shared_ptr<BaseResult> result );


    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling)
    {Error("InitCoupling not implemented",__FILE__,__LINE__);};
  
    //! calculate coupling terms
    void CalcOutputCoupling()
    {Error("CalcOutputCoupling not implemented",__FILE__,__LINE__);};

    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output) { return false; }


  protected:

    //! Init the time stepping
    void InitTimeStepping();

    //! Define available result types
    void DefineAvailResults();

    // ========================
    // time stepping
    // ========================
    // solving of nonlinear acoustics
    NodeStoreSol<Double> sol_der1Array_, sol_der2Array_;

    // ========================
    // Postprocessing results
    // ========================

    BaseNodeStoreSol * press_;    //!< pressure

    NodeStoreSol<Double> solDeriv1_; //!< contains 1st derivative of solution
    NodeStoreSol<Double> solDeriv2_; //!< contains 2nd derivative of solution

  private:

    // reads in the PML data
    void ReadDataPML(std::string& typePML, Matrix<Double>& inner, 
		     Double& dampPML, ParamNode * actNode);
    
    //! computes the PML layer dimensions
    void GetPMLLayerData(Matrix<Double>& inner, Matrix<Double>& outer,
			 RegionIdType regionId );

    //!
    bool savePress_;

    //!
    bool savePressHist_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class AcousticXYZPDE
  //! 
  //! \purpose 
  //! This class is derived from class SinglePDE.
  //! It is used for solving acoustic equation on one time step.  
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
