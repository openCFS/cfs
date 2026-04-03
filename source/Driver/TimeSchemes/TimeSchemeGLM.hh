// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TimeSchemeGLM.hh
 *       \brief    Provides class for a GLM based time stepping scheme
 *
 *       \date     03/01/2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef TIMESCHEMEGLM_HH_
#define TIMESCHEMEGLM_HH_

#include <fstream>
#include <set>

#include "BaseTimeScheme.hh"
#include "GLMSchemeLib.hh"
#include "Utils/mathParser/mathParser.hh"


namespace CoupledField{


//see below for description
class TimeSchemeGLM : public BaseTimeScheme{
  public:

    typedef enum{
      NONE,
      INCREMENTAL,
      TOTAL
    } NonLinType;



    /*!
     *  Constructor of the GLM scheme
     *  \param[in] type The TimeScheme to be used. Newmark, Trapezoidal, etc.
     *  \param[in] solDerivOrder The time derivative order of the solution to the effective system
     */
    TimeSchemeGLM(GLMScheme::SchemeType type, UInt solDerivOrder=0, TimeSchemeGLM::NonLinType nlType=NONE);
    
    
    /*!
     * Alternative constructor to directly pass a pre-constructed GLM-scheme
     * \param[in] scheme Externally created time scheme. Ownership gets handed to this class.
     * \param[in] solDerivOrder The time derivative order of the solution to the effective system
     */
    TimeSchemeGLM(GLMScheme* scheme, UInt solDerivOrder=0, TimeSchemeGLM::NonLinType nlType=NONE);
    
    //! Copy constructor
    TimeSchemeGLM(const TimeSchemeGLM& ts);

    virtual ~TimeSchemeGLM();

    //! \copydoc BaseTimeScheme::Init(SingleVector*,Double)
    virtual void Init(SingleVector* solVec,Double dt);

    //! \copydoc BaseTimeScheme::ModifyInit(bool)
    virtual void ModifyInit(bool extrapolateStatic);

    //! \copydoc BaseTimeScheme::BeginStep(bool)
    virtual void BeginStep(bool updatePredictor=true, bool storeInitialIterGlmVector=false);

    //This function is pretty messy right now and we need to reconsider
    // mainly because of the many if clauses to realize a optional predictor scheme...
    //! \copydoc BaseTimeScheme::ComputeStageRHS(UInt,Integer,SingleVector*,Integer)
    virtual void ComputeStageRHS(UInt actStage, Integer derivId, SingleVector* rhsVec, Integer subIdx=-1, bool skipIncremental=false);

    virtual void UpdateStageRHSWithVector(UInt actStage, Integer derivId, SingleVector* rhsVec,
                                            SingleVector* UpdateVector, Double factor, bool forceReset = false);

    /// Update the GLM Vectors according to new solution
    virtual void FinishStep( );

    // In the case that we did not converge, we have to reset the glmVec (and the last solution stored to the feFunction in order to start with the correct values in the next sub-iteration
    // If the iteration did converge, we simply free some memory by deleting the old glmVec
    virtual void ProcessGlmVec(bool converged=false);

    // Triggers the update of the glmVector
    void ResetGlmVector() {
      resetGlmVector_ = true;
    }

    //! \copydoc BaseTimeScheme::SetSolutionTimeDerivOrder(UInt,Double)
    virtual void SetSolutionTimeDerivOrder(UInt order,Double timeStepSize){
      solOrder_ = order;
      curScheme_->ComputeCoefficients(order,timeStepSize);
    }

    //! Obtain number of stages of the scheme
    UInt GetNumStages(){
      return curScheme_->numStages_;
    }

    //! Obtain reference to current stage vector to avoid copying of elements
    SingleVector * GetStageVector(UInt stage){
      assert(stage < curScheme_->numStages_);
      return stageVector_[stage];
    }

    //! Obtain reference to current GLM vector to avoid copying of elements
    SingleVector* GetGLMVector(UInt numSol){
      return glmVector_[numSol];
    }

    //! Obtain reference to initial GLM vector to avoid copying of elements
    SingleVector*  GetInitialIterGLMVector(UInt numSol){
      return initialIterGlmVector_[numSol];
    }

    //! Obtain the size of the GLM vector
    UInt GetSizeGLMVector(){
      return curScheme_->sizeGLMVec_;
    }

    void reset_dt();

    //! \copydoc BaseTimeScheme::AddMatFactors(UInt,const std::map<FEMatrixType,Integer> &,std::map<FEMatrixType,Double> &)
    virtual void AddMatFactors(UInt stage, const std::map<FEMatrixType,Integer> & matMap,
                                  std::map<FEMatrixType,Double> & matFactors);

    //! \copydoc BaseTimeScheme::AdaptBC(Double&,Double,UInt,Integer)
    virtual void AdaptBC(Double& transVal, Double initValue,UInt initDerivOrder, Integer eqnNumber){
      transVal = curScheme_->TransformBC(glmVector_,initValue,initDerivOrder, eqnNumber);
    }

    virtual SingleVector* GetTimeDerivative(UInt order);

    virtual void SetTimeDerivVector(UInt order,SingleVector * coefVector);

    /// Give the timestep the possibility to initialize
    virtual void InitStage(UInt aStage,Double aTime,Domain* domain){
      curScheme_->PrepareStage(aStage,aTime, domain);
    };
    
    bool isIncremental(){
      if(nLinType_ == INCREMENTAL){
        return true;
      } else {
        return false;
      }
    }
    
    void forceIncremental(){
      nLinType_ = INCREMENTAL;
    }

    double standartStepsize(bool * accepted);
    double smoothStepsize(bool * accepted);

  protected:
    // Pointer to Mathpaerser Instance
    MathParser * mathparser_ = nullptr;

    void InitGLMs();

    ///map of all available time schemes
    std::map<GLMScheme::SchemeType, GLMScheme*> availSchemes;

    /// pointer to time scheme
    GLMScheme* curScheme_;

    /// type of timescheme
    GLMScheme::SchemeType curType_;

    //!This is the input vector of the GLM how its components get a meaning
    //!in combination with the GLM scheme used namely the parameters numOldSol_ numsolDerivs_ etc.
    StdVector< SingleVector* > glmVector_;

    // This is just a copy of the initial glmVector from the sub-iteration which is needed for the transition between iterations
    // More detailed description: The finishStep function overwrites the feFunction with the current (usually not converged) solution in order to calculate the norms.
    // In the next time step, this solution is used for the time scheme, although e.g. a standard linear iteration should use the same vectors/matrices as before and only update the RHS from a linear form.
    // In order to avoid this, we store the old glmVec seperately and "undo" the last part of finishStep by resetting the glmVec
    StdVector< SingleVector* > initialIterGlmVector_;

    // Bool to change the behaviour when reading in an initialState from a previous sequence step
    // If it is set to true, we assume that all primary values for all steps needed for the time scheme are equal to the static one.
    bool extrapolateStatic_ = false;

    // Bool if we have to reset the glmVector
    bool resetGlmVector_ = false;


    ///Stores for each stage, for each time derivative the stage values
    StdVector< SingleVector* > stageVector_;

    ///For the basic schemes we can store some references to the right hand sides which are in fact
    ///the predictor values
    StdVector< SingleVector* > predictors_;

    ///Stores the index in the GLM vector which is connected to the stage vector
    Integer avoidUpdateIdx_;

    ///Stores for each step if the predictors are calculated
    StdVector<bool> predictorCalculated_;

    std::set<UInt> avoidFreeingIdx_;

    ///Store the type of nonlinearity to be considered in the scheme
    NonLinType nLinType_;

    // Stores previus y_[-1], used when adaptive Timestepping
    SingleVector* prevPrevSol_;

    int adaptiveStepCount_;

  private:

    void LTELocalErrorEstimation();
    bool ComputeAdaptiveStepSize();



    ///just export the scheme to a file
    void ExportGLM(string pdeName, int feFctId, int curStep, int coupleIter){
      std::string fname = "glmExport_" + pdeName + "_feFctId" + std::to_string(feFctId) +  "_step" + std::to_string(curStep) + "_coupleIter" + std::to_string(coupleIter) + ".txt";
      std::fstream myfile(fname,  std::ios::out);
      myfile << "This is the GLM Vector" << std::endl;
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        myfile << "Index " << i << std::endl;
        myfile << glmVector_[i]->ToString(TS_NONZEROS,"\n") << std::endl;
        myfile << "Finish GLM Vector" << std::endl;
      }
      myfile << std::endl;

      myfile << "This is the initial GLM Vector of this time step" << std::endl;
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        myfile << "Index " << i << std::endl;
        myfile << initialIterGlmVector_[i]->ToString(TS_NONZEROS,"\n") << std::endl;
        myfile << "Finish initial GLM Vector" << std::endl;
      }
      myfile << std::endl;
      myfile.close();
    }

};

/*! \class TimeSchemeGLM
 *     \brief This class represents methods to realize a GLM.
 *     @author A. Hueppe
 *     @date 01/2012
 *
 *    The idea begind a general linear method is that many time stepping schemes can be formulated as a
 *    tableau of the structure
 *    \f[
 *     \begin{array}{c|c}
 *     \mathbf{a} & \mathbf{U} \\
 *     \hline
 *     \mathbf{b} & \mathbf{V} \\
 *     \end{array}
 *    \f]
 *    The entries in \f$\mathbf{a}\f$ correspond to the stages of the scheme and the column entries in
 *    \f$\mathbf{U}\f$ correspond to the entries of the GLM vector. The concrete Schemes are defined in
 *    the GLMSchemeLib.hh by giving their tableaus. The implementation makes use of two optimization flags
 *    which minimize the memory requirements and are provided by the GLMScheme
 *
 *    \li \c lastStageIsSolution_ If this flag is set, the stage vector of the last stage points directly to the corresponding
 *           entry in the GLMVector.
 *    \li \c usePredictors_ If this flag is set, we store the result of the RHS calculation in an internal variable
 *            In the Finish stage method we reuse these to update the entrys in the GLM vector thus saving
 *            some recalculation
 *
 *    The known schemes: Trapezoidal,Newmark and HHT use both flags no matter which formulation used
 *    (effMass,effStiff) etc.
 */
}

#endif /* TIMESCHEMEGLM_HH_ */
