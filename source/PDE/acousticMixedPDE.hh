#ifndef FILE_NEWBASEACOUSTICMIXEDPDE
#define FILE_NEWBASEACOUSTICMIXEDPDE

#include "SinglePDE.hh"
 
namespace CoupledField
{

  //! Forward class declarations
  class BaseForm;


  //! Class for acoustic PDE using velocity and pressure
  class AcousticMixedPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
//     Grid *aGrid; //pointer to grid
//     TimeFunc *aTimeFunc; //pointer to class WriteResults. output data.
//     WriteResults *aOutFile; //pointer to class TimeFunc


    AcousticMixedPDE(Grid* aGrid, ParamNode* paramNode );

    //!  Deconstructor
    virtual ~AcousticMixedPDE();

   //! Define availabe result types
    void DefineAvailResults();

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( );

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //! Read special boundary conditions
    void ReadSpecialBCs();
    
    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> results );

   //! do PostProcessing step
    virtual void PostProcess( );

    //! initalize PDE coupling
    virtual void InitCoupling(PDECoupling * Coupling) {;};

    //! Fill in input coupling terms
    virtual void CalcInputCoupling() {;};
  
    //! calculate coupling terms
    virtual void CalcOutputCoupling() {;};

  protected:
 
    //! Init the time stepping
    void InitTimeStepping();

  private:

    //! true, if solution should be written to result file
    bool saveSolVel_;
    bool saveSolPres_;

  };

} // end of namespace
#endif

