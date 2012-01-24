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


class BaseTimeScheme{

  public:

    BaseTimeScheme(){
      //This is the default which corresponds to effectiveStiffness formulation
      solOrder_ = 0;
    }

    virtual ~BaseTimeScheme(){

    }

    /*! Initializes all required vectors and variables
     *  \param[in] solVec pointer to the feFunction ceofficient vector
     *  \param[in] dt current timestep size
     */
    virtual void Init(SingleVector* fncCoef,Double dt)=0;

    /*!
     *   Computes the effective RHS based on the GLM vector and preceeding stage solutions
     *   \param[in] actStage The current stage number (used to determine row in GLM tableau)
     *   \param[in] derivId The current derivative i.e. matrix stiffness->derivid = 0, damping->derivId=1 etc.
     *   \param[out]rhsVec OutputVector
     *   \param[in] subIdx not used right now, preparing for local time stepping
    */
    virtual void ComputeStageRHS(UInt actStage, Integer derivId, SingleVector* rhsVec, Integer subIdx=-1)=0;

    /// Update function called at the end of the solvestep
    virtual void FinishStep()=0;

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

    /*! Transforms a given boundary condition value according to the given scheme
     *   \param[out] transVal the output value
     *   \param[in] initValue input value of BC
     *   \param[in] initDerivOrder Does the BC specify the unknown itselfs or the initDerivOrder-th time derivative
     *   \param[in] eqnNumber to obtain the correct values from GLM vector
    */
    virtual void AdaptBC(Double& transVal, Double initValue,UInt initDerivOrder, Integer eqnNumber) = 0;

    /// Obtain the time derivative order of the systems solution
    UInt GetSolutionTimeDerivOrder(){
      return solOrder_;
    }

  protected:

    /// Current time derivative order of the solution
    UInt solOrder_;


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
