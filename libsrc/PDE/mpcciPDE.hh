#ifndef FILE_MPCCIPDE_NEW
#define FILE_MPCCIPDE_NEW

#include "scalarnodeEQN.hh"
#include "SinglePDE.hh" 

#ifdef MpCCI
#include <MpCCIcpl/MpCCIexch.hh>
#endif

namespace CoupledField
{

  //! Class for coupling a pde via MpCCI
  /*! 
    This class is derived from class BasePDE. It is used for coupling CFS++ via MpCCI
  */

  class MpcciPDE : public SinglePDE {

  public:

    //! Constructor. here we read integration parameters
    /*!
      \param 
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    MpcciPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut);

    //! Deconstructor
    virtual ~MpcciPDE(){};

    virtual void Init(Integer sequenceStep = 0,
		      std::string  bcSequenceTag = "anyTag");


    virtual void PreparePDE4Computation();

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators();

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //! return size of solution
    virtual Integer getSize() const 
    { return numPDENodes_*dofspernode_;}

    //! initialize time stepping: 
    //! nothing to do in mpcci!
    virtual void InitTimeStepping(){ ENTER_FCN( "MpcciPDE::InitTimeStepping"); };

    //! set current time step
    virtual void SetTimeStep(const Double dt){ ENTER_FCN( "MpcciPDE::SetTimeStep");};

   
    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! do PostProcessing step
    virtual void PostProcess();

    virtual void WriteGeneralPDEdefines()
    { ENTER_FCN( "MpcciPDE::WriteGeneralPDEdefines");};

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  

    virtual void WriteResultsInFile(const Integer kstep = 0,
				    const Double asteptime = 0.0,
				    Integer stepOffset = 0,
				    Double timeOffset = 0.0);
    

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! initalize PDE coupling
    virtual void InitCoupling(PDECoupling * Coupling);

    virtual void CalcInputCoupling();

    //! calculate coupling terms
    virtual void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output);
  

  protected:

    StdVector<StdVector<Elem*> > F_Interface_; //!<vector of vectors conaining Elements with acting force
    StdVector<StdVector<StdVector<ShortInt> > > isBoundaryNode_; //!< vector containing flag array for element boundary nodes
    StdVector<StdVector<StdVector<Integer> > > elemNodeToCouplingNode_; //!< assigns each coupling element node the according Coupling Node number
    
    // *****************
    //  POSTPROCESSING
    // *****************

    // for check: own solver
    Boolean SolverCFS_; //<! parameter indicator: TRUE, if you want to use Solver CFS. reading from config-file
    Matrix<Double> sysmat_;
    Vector<Double> vecrhs_;

  private:
    //!MpCCI
   StdVector<Integer> mapSD_;
#ifdef MpCCI
    MpCCIexch * ptMpCCIexch_;
#endif
    Boolean flagFirstTimeStep_; // flag for first time stewp
    Integer MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
    void ReadStoreResults();

  };

} // end of namespace

#endif
