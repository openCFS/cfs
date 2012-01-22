/*
 * TimeSchemeGLM.hh
 *
 *  Created on: Jan 3, 2012
 *      Author: ahueppe
 */

#ifndef TIMESCHEMEGLM_HH_
#define TIMESCHEMEGLM_HH_

#include <fstream>

#include "BaseTimeScheme.hh"
#include "GLMSchemeLib.hh"

namespace CoupledField{

class TimeSchemeGLM : public BaseTimeScheme{
  public:

    typedef enum{
      TRAPEZOIDAL = 1,
      NEWMARK = 2,
      RK4 = 3,
      BDF2 = 4
    } SchemeType;


    //Solderiv order is default set to zero which is the effective stiffness formulation
    TimeSchemeGLM(SchemeType type, UInt solDerivOrder=0);

    virtual ~TimeSchemeGLM();

    virtual void Init(SingleVector* solVec,Double dt);

    virtual void ComputeStageRHS(UInt actStage, Integer derivId, SingleVector* rhsVec, Integer subIdx=-1);

    virtual void FinishStep();

    virtual void SetSolutionTimeDerivOrder(UInt order,Double timeStepSize){
      solOrder_ = order;
      curScheme_->ComputeCoefficients(order,timeStepSize);
    }

    UInt GetNumStages(){
      return curScheme_->numStages_;
    }

    SingleVector * GetStageVector(UInt stage){
      //assert(stage < curScheme_->numStages_);

      return stageVector_[stage];
    }

    virtual void AddMatFactors(UInt stage, const std::map<FEMatrixType,Integer> & matMap
                               , std::map<FEMatrixType,Double> & matFactors);

    virtual void AdaptBC(Double& transVal, Double initValue,UInt initDerivOrder, Integer eqnNumber){
      transVal = curScheme_->TransformBC(glmVector_,initValue,initDerivOrder, eqnNumber);
    }

  protected:

    void InitGLMs();

    std::map<SchemeType, GLMScheme*> availSchemes;

    GLMScheme* curScheme_;

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

  private:

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
}

#endif /* TIMESCHEMEGLM_HH_ */
