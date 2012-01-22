/*
 * BaseTimeScheme.hh
 *
 *  Created on: Jan 3, 2012
 *      Author: ahueppe
 */

#ifndef BASETIMESCHEME_HH_
#define BASETIMESCHEME_HH_

#include "General/Environment.hh"
#include "MatVec/SingleVector.hh"

namespace CoupledField{


class BaseTimeScheme{

  public:

    BaseTimeScheme(){
      //This is the default which corresponds to effectiveStiffness formulation
      solOrder_ = 0;
    }

    virtual ~BaseTimeScheme(){

    }

    virtual void Init(SingleVector* fncCoef,Double dt)=0;

    virtual void ComputeStageRHS(UInt actStage, Integer derivId, SingleVector* rhsVec, Integer subIdx=-1)=0;

    virtual void FinishStep()=0;

    virtual void SetSolutionTimeDerivOrder(UInt order,Double timeStepSize)=0;

    virtual void AddMatFactors(UInt stage, const std::map<FEMatrixType,Integer> & matMap
                               , std::map<FEMatrixType,Double> & matFactors)=0;

    //virtual SingleVector* GetSolution(UInt actStage)=0;

    virtual UInt GetNumStages()=0;

    virtual SingleVector * GetStageVector(UInt stage)=0;

    virtual void AdaptBC(Double& transVal, Double initValue,UInt initDerivOrder, Integer eqnNumber) = 0;

    UInt GetSolutionTimeDerivOrder(){
      return solOrder_;
    }

  protected:

    UInt solOrder_;


};

}

#endif /* BASETIMESCHEME_HH_ */
