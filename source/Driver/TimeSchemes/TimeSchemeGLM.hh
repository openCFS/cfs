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

#include "BaseTimeScheme.hh"
#include "GLMSchemeLib.hh"

namespace CoupledField{

//see below for desciption
class TimeSchemeGLM : public BaseTimeScheme{
  public:

    /// Enumeration for each GLM scheme available
    typedef enum{
      TRAPEZOIDAL = 1,
      NEWMARK = 2
     // RK4 = 3,
     // BDF2 = 4
    } SchemeType;


    /*!
     *  Constructor of the GLM scheme
     *  \param[in] type The TimeScheme to be used. Newmark, Trapezoidal, etc.
     *  \param[in] solDerivOrder The time derivative order of the solution to the effective system
     */
    TimeSchemeGLM(SchemeType type, UInt solDerivOrder=0);

    virtual ~TimeSchemeGLM();

    //! \copydoc BaseTimeScheme::Init(SingleVector*,Double)
    virtual void Init(SingleVector* solVec,Double dt);

    //This function is pretty messy right now and we need to reconsider
    // mainly because of the many if clauses to realize a optional predictor scheme...
    //! \copydoc BaseTimeScheme::ComputeStageRHS(UInt,Integer,SingleVector*,Integer)
    virtual void ComputeStageRHS(UInt actStage, Integer derivId, SingleVector* rhsVec, Integer subIdx=-1);

    /// Update the GLM Vectors according to new solution
    virtual void FinishStep();

    //! \copydoc BaseTimeScheme::SetSolutionTimeDerivOrder(UInt,Double)
    virtual void SetSolutionTimeDerivOrder(UInt order,Double timeStepSize){
      solOrder_ = order;
      curScheme_->ComputeCoefficients(order,timeStepSize);
    }

    /// Obtain number of stages of the scheme
    UInt GetNumStages(){
      return curScheme_->numStages_;
    }

    /// Obtain reference to current stage vector to avoid copying of elements
    SingleVector * GetStageVector(UInt stage){
      assert(stage < curScheme_->numStages_);
      return stageVector_[stage];
    }

    //! \copydoc BaseTimeScheme::AddMatFactors(UInt,const std::map<FEMatrixType,Integer> &,std::map<FEMatrixType,Double> &)
    virtual void AddMatFactors(UInt stage, const std::map<FEMatrixType,Integer> & matMap
                               , std::map<FEMatrixType,Double> & matFactors);

    //! \copydoc BaseTimeScheme::AdaptBC(Double&,Double,UInt,Integer)
    virtual void AdaptBC(Double& transVal, Double initValue,UInt initDerivOrder, Integer eqnNumber){
      transVal = curScheme_->TransformBC(glmVector_,initValue,initDerivOrder, eqnNumber);
    }

  protected:

    void InitGLMs();

    ///map of all available time schemes
    std::map<SchemeType, GLMScheme*> availSchemes;

    /// pointer to time scheme
    GLMScheme* curScheme_;

    /// type of timescheme
    SchemeType curType_;

    //!This is the input vector of the GLM how its components get a meaning
    //!in combination with the GLM scheme used namely the parameters numOldSol_ numsolDerivs_ etc.
    StdVector< SingleVector* > glmVector_;

    ///Stores for each stage, for each time derivative the stage values
    StdVector< SingleVector* > stageVector_;

    ///For the basic schemes we can store some references to the right hand sides which are in fact
    ///the predictor values
    StdVector< SingleVector* > predictors_;

    ///Stores the index in the GLM vector which is connected to the stage vector
    Integer avoidUpdateIdx_;

    ///Stores for each step if the predictors are calculated
    StdVector<bool> predictorCalculated_;

  private:

    ///just export the scheme to a file
    void ExportGLM(){
      std::fstream myfile("glmExport.txt",  std::ios::out);
      myfile << "This is the GLM Vector" << std::endl;
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        myfile << "Index " << i << std::endl;
        myfile << glmVector_[i]->ToString(1,'\n') << std::endl;
        myfile << "Finish GLM Vector" << std::endl;
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
