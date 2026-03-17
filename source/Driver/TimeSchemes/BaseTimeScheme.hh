// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TimeSchemeGLM.hh
 *       \brief    Basic Timestepping interface classes
 *
 *       \date     Jan 3, 2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef BASETIMESCHEME_HH_
#define BASETIMESCHEME_HH_

#include "General/Environment.hh"
#include "MatVec/SingleVector.hh"

namespace CoupledField{
class Domain;

class BaseTimeScheme{

  public:

    bool isAdaptive;

    BaseTimeScheme(){
      //This is the default which corresponds to effectiveStiffness formulation
      solOrder_ = 0;
    }


    virtual ~BaseTimeScheme(){

    }

    /*! Initializes all required vectors and variables
     *  \param[in] solVec pointer to the feFunction coefficient vector
     *  \param[in] dt current timestep size
     */
    virtual void Init(SingleVector* fncCoef,Double dt)=0;

    /*! Modifies the initialization (for initialState usage)
     * Here we can modify the behaviour before the actual initialization based on user configuration
     */
    virtual void ModifyInit(bool extrapolateStatic)=0;

    /*! Begin a new step
     * 
     * \param[in] updatePredictor  Flag, if predictor values get re-calculated. In case of
     *                             nonlinear solution, this has to be set to false
     * \param[in] storeInitialIterGlmVector Flag, if the initial iteration GLM vector should be stored.
     *                                      This is used for the first iteration of a coupled problem or a nonlinear problem
     */
    virtual void BeginStep(bool updatePredictor = true, bool storeInitialIterGlmVector=false)=0;

    /*!
     *   Setter to Allow PDE`s to set the Domain, so there is MathParser access in TimeSchemeGLM. -> Only set by Smooth PDE , Needed for controll of adaptive Timestepping.
     */
    void SetDomain(Domain* d) { domain_ = d; }

    /*!
     *   Computes the effective RHS based on the GLM vector and preceeding stage solutions
     *   \param[in] actStage The current stage number (used to determine row in GLM tableau)
     *   \param[in] derivId The current derivative i.e. matrix stiffness->derivid = 0, damping->derivId=1 etc.
     *   \param[out]rhsVec OutputVector
     *   \param[in] subIdx not used right now, preparing for local time stepping
    */
    virtual void ComputeStageRHS(UInt actStage, Integer derivId, SingleVector* rhsVec, Integer subIdx=-1,bool forceIncremental=false)=0;

    virtual void UpdateStageRHSWithVector(UInt actStage, Integer derivId, SingleVector* rhsVec,
                                 SingleVector* UpdateVector, Double factor, bool forceReset = false)=0;
    
    /// Update function called at the end of the solvestep
    virtual void FinishStep()=0;


    // Update function that processes the glmVector in the case of a GLM-scheme
    virtual void ProcessGlmVec(bool converged=false)=0;

    // Triggers the update of the glmVector
    virtual void ResetGlmVector()=0;

    /*! Change the formulation of the scheme e.g. from effective mass to effective stiffness
     *  \param[in] order The order of time derivative of the systems solution
     *  \param[in] timeStepSize the current time step size
    */
    virtual void SetSolutionTimeDerivOrder(UInt order,Double timeStepSize)=0;

    /*!
     *  obtain coefficients for the construction of the effective system matrix
     *  \param[in] stage the current stage
     *  \param[in] matMap mapping of FeMatrixTypes to a certain time derivative order
     *  \param[out] matFactors The resulting map to which the coefficients of the scheme get added
    */
    virtual void AddMatFactors(UInt stage, const std::map<FEMatrixType,Integer> & matMap
                               , std::map<FEMatrixType,Double> & matFactors)=0;

    /// Obtain number of stages of the scheme
    virtual UInt GetNumStages()=0;

    /// Obtain reference to current stage vector to avoid copying of elements
    virtual SingleVector * GetStageVector(UInt stage)=0;

    //! Obtain reference to current GLM vector to avoid copying of elements
    virtual SingleVector* GetGLMVector(UInt numSol)=0;

    //! Obtain reference to current GLM vector to avoid copying of elements
    virtual SingleVector* GetInitialIterGLMVector(UInt numSol)=0;

    //! Obtain the size of the GLM vector
    virtual UInt GetSizeGLMVector()=0;

    /*! Transforms a given boundary condition value according to the given scheme
     *   \param[out] transVal the output value
     *   \param[in] initValue input value of BC
     *   \param[in] initDerivOrder Does the BC specify the unknown itselfs or the initDerivOrder-th time derivative
     *   \param[in] eqnNumber to obtain the correct values from GLM vector
    */
    virtual void AdaptBC(Double& transVal, Double initValue,
                         UInt initDerivOrder, Integer eqnNumber) = 0;
    virtual void AdaptBC(Complex& transVal, Complex initValue,
                         UInt initDerivOrder, Integer eqnNumber) {
      // nothing to do for complex valued peoblems
    }

    /*! Obtain the time derivative of the solution
     *   \param[in] order requested order of time derivative
     */
    virtual SingleVector* GetTimeDerivative(UInt order)=0;

    /*! Set time derivative vector for a certain order
     *   \param[in] order of time derivative
     *   \param[in] coefVector coefficient vector
     */
    virtual void SetTimeDerivVector(UInt order,SingleVector * coefVector)=0;

    /// Obtain the time derivative order of the systems solution
    UInt GetSolutionTimeDerivOrder(){
      return solOrder_;
    }

    /// Give the timestep the possibility to initialize
    virtual void InitStage(UInt aStage,Double aTime,Domain* domain)=0;

    virtual bool isIncremental()=0;

    virtual void forceIncremental()=0;

    virtual void ExportGLM(string pdeName, int feFctId, int curStep, int coupleIter)=0;
    
  protected:

    // Current time derivative order of the solution
    UInt solOrder_;

    //Domain needed for accces to Mathparser in TimeScheme, (Used in Adaptive timestepping)
    Domain* domain_ = nullptr;


};

/*! \class BaseTimeScheme
 *     \brief Basic interface for all TimeStepping schemes in CFS++
 *     @author A. Hueppe
 *     @date 01/2012
 *
 *    The idea is that every linear timestep can be subdivided into a stage loop and an
 *    update step. For each stage we do the following:
 *    \li Construct a system matrix
 *    \li Create a RHS given by the PDE or the user
 *    \li Add entrys to the RHS according to the time scheme (create effective matrix)
 *    \li Transform and Impose boundary conditions
 *    \li Solve The system
 *
 *    Finally we update everything and the step is finished.
 *
 *    How the single steps are performed determines the timescheme. Right now we try to put as
 *    much as possible in a general linear method (GLM) which is an abstract notation applicable
 *    to many different time schemes. Anyhow it is still possible to overlaod this bas class by a
 *    specialized one for a specific scheme. The big advantage ist that the main loop in
 *    StdSolveStep does not change
 */

}

#endif /* BASETIMESCHEME_HH_ */
