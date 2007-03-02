#ifndef FILE_BUBBLE_PDE_HH
#define FILE_BUBBLE_PDE_HH

#include "SinglePDE.hh" 

namespace CoupledField
{
  
  // forward class declaration
  class BubbleODE;
  class BaseODESolver;
  

  //! Class for wrapping ODE for bubble radius computation
  class BubblePDE : public SinglePDE {

  public:

    //! Constructor. here we read integration parameters
    /*!
      \param 
      \param aGrid pointer to grid
      \param aGrid pointer to class Grid
    */
    BubblePDE (Grid * aptgrid, ParamNode* paramNode );

    //! Destructor
    ~BubblePDE();

    //! Define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping() {;};

    //! Nothing to do
    void SetTimeStep(const Double dt) {;};

    // ======================================================
    // SOLVING SECTION
    // ======================================================

    //! Solve one time step of the ODE 
    void Solve();
    

    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! Do PostProcessing step
    void CalcResults( shared_ptr<BaseResult> result );

     // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! Initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! Fill in input coupling terms
    void CalcInputCoupling();

    //! Calculate coupling terms
    void CalcOutputCoupling();

    //! Returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);
  
  private:

    //! Define available results
    void DefineAvailResults();

    //! Calculate RHS values for acoustic PDE
    void CalcAcouRHS( StdVector<std::string>& regions, StdVector<UInt>& nodes,
                      Vector<Double>& values );

    //! Allocate the radius and velocity solutions to the right coupling elements
    void CalcBubbleRadius( StdVector<Elem*> & elems,
				      Vector<Double>& couplVals,
				      SolutionType solType );

    //! Output stream for writing results of one element
    std::ofstream outValues_;

    //! Pointer to class of bubbledynamical method for each local element
    StdVector<BubbleODE*> ptBubble_;

    //! Pointer to the solver class for ode's
    BaseODESolver *ptODESolver_;

    //! Vector contains on input the initial data of radius und velocity
    //! for the bubbledynamic, on the output the solution 
    StdVector<Double> bubbleValues_;

    //! Vector contains radius of bubbles for each element
    StdVector<Double> radius_;
    //! Vector contains radius of bubbles at iteration end of last time step
    StdVector<Double> radiusOldStep_;

    //! Vector contains velocity of bubble walls for each element
    StdVector<Double> velocity_;
    //! Vector contains velocity of bubble walls at iteration end of last time step
    StdVector<Double> velocityOldStep_;

    // In actual cavitation model: bubbledensity is constant in time, 
    // but may vary in space; in complexer models it can vary in time too.
    // Attention: for now assumption that density is constant in space. 
    //! Density of number of bubbles per unit volume
    Double bubbleDensity_;

    //! Attribute describing model for bubble dynamics
    BubbleDynType bubbleDynType_;


    //! Parameters for the bubble dynamics
    //! \param initRadius_ Initial radius of bubble
    //! \param initVel_ Initial velocity of bubble
    //! \param density_ Density of fluid
    //! \param sonicVel_ Sonic velocity of fluid
    //! \param pStatic_ Static pressure of fluid
    //! \param pVapour_ Vapour pressure of fluid
    //! \param surfaceTension_ Surface tension of fluid
    //! \param polytrop_ Polytropic index of fluid
    //! \param viscosity_ Viscosity of fluid
    Double initRadius_;
    Double initVel_;
    Double density_;
    Double sonicVel_;
    Double pStatic_;
    Double pVapour_;
    Double surfaceTension_;
    Double polytrop_;
    Double viscosity_;


    //! Parameter for the dimensionless case
    StdVector<Double> y_;
    Double tDim_;
    Double t_;
    Double pressure_;
    Double pressureDeriv_;  

    //! Suggested step size for ODESolverfor each element
    StdVector<Double> hTry_;

    // ====== Excitation variables ======
    
    //! \param pressureAmpl_ Amplitude of pressure, needed only
    //! if pressure is directly computed
    //! \param frequency_ Frequency of pressure, needed only
    //! if pressure is directly computed
    Double pressureAmpl_;
    Double frequency_;


    Double dt_;
    Double dtStep_;
    Double steptime_;

    // buffer for pressure input
    Vector<Double> * pressBuf_;

    // buffer for first derivative of pressure input
    Vector<Double> * pressDerivBuf_;

    //! vector with coupling elements
    StdVector<Elem*> * couplElems_;

    // ======================================================
    // STORERESULTS SECTION
    // ======================================================    

    //! Storesolution for vectorial result ( radius, velocity, volumeFraction )
    ElemStoreSol<Double> elemResult_;

    //! Storesolution for vecotrial result( RHS, 2nd derivative of radius )
    ElemStoreSol<Double> addElemResult_;

    StdVector<std::string> saveElemBubbleHist_;//!< name of elements to be saved
    StdVector<std::string> saveElemBubbleHistRegion_;//!< name of elements to be saved
    StdVector<std::string> saveElemRHSHist_;//!< name of elements to be saved
    StdVector<std::string> saveElemRHSHistRegion_;//!< name of elements to be saved

    //! Flags indicating writing of results
    bool writeValues_, writeRHS_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class BubblePDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! electrostatic equation in 3D. 
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
