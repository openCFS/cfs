#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>

#include "SolveStepHyst.hh"
#include "Driver/Assemble.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Results/BaseResults.hh"
#include "Utils/EvalIntegrals/BiotSavart.hh"
#include "Utils/Timer.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
//#include "MatVec/SingleVector.hh"

namespace CoupledField {
  
  // declare logging stream
  DECLARE_LOG(solvestephyst)
  DEFINE_LOG(solvestephyst, "solvestephyst")
  
  SolveStepHyst::SolveStepHyst(StdPDE & apde)
  :StdSolveStep(apde)
  {
  }
	
  //! Destructor
  SolveStepHyst::~SolveStepHyst() {
    //logFile_.close();
  }
  
//  void SolveStepHyst::SetLastItOrLastTSSBMVectors(bool lastTS){
//    /*
//     *  Function needed in NonLinHysteresis case
//     *
//     *  it stores the current values of rhs, residual, solution, ...
//     *  for the next iteration or next ts by copying it to the corresponding
//     *  SBMVectors
//     */
//    
//    if(lastTS){
//      // oldTS > values after last iteration of previous TS
//      // to be stored after iteration suceeded
//      oldTSLinRhsVec_ = currentLinRhsVec_;
//      oldTSNonLinRhsVec_ = currentNonLinRhsVec_;
//      oldTSSolVec_ = solVec_;
//    } else {
//      // oldIt > values of the current TS but from previous It
//      // during first iteration of a new TS, these vectors contain the values
//      // after the last iteration of the previous TS (similar as oldTS...)
//      oldItNonLinRhsVec_ = currentNonLinRhsVec_;
//      oldItResVec_ = currentResVec_;
//    }
//  }
//  
  
  void SolveStepHyst::CalcResidualAndConfigSystemForHysteresis(SBM_Vector& oldSolution,SBM_Vector& solIncrement,SBM_Vector& stageSol, Double usedEta,
          UInt stage, UInt callingCnt, UInt evalVersion, bool trans, bool forceReevaluation,
          bool skipReassembly, bool debugOutput, bool reset){
    /*
     * new parateter added 21.7: skipReassembly
     * > if true, residual will be calculated but new system will not get assembled
     * > useful during line search as we do not use the newly assembled system in that case (we simply overwrite it during
     * next eta)
     * > set this flag to false after the final eta has been found
     *
     *
     * reset will have the same effect as callingCnt = 0
     */
    
    //UInt callingCntBackup = callingCnt;
    
    
    if(debugOutput){
      std::cout << "<<<< CalcResidualAndConfig >>>>" << std::endl;
      std::cout << "CallingCnt: " << callingCnt << std::endl;
      std::cout << "evalVersion: " << evalVersion << std::endl;
    }
    
    // IMPORTANT: reallow reevaluation of hyst operator (as input may have changed)
    PDE_.SetFlagInCoefFncHyst("resetReeval",true);
    
    
    if(reset){
      callingCnt = 0;
    }
    
    /*
     * flags to be set in hysteresis operator to enable different
     * evalVersions
     */
    /*
     * tensorToReturn---
     * 1: assemble matrix with epsInitial / nuInitial
     * 2: assemble matrix with eps0 / nu0
     * 3: assemble deltaMat
     */
    UInt tensorToReturnLHS;
    UInt tensorToReturnRHS;
    
    /*
     * tensorToAdd---
     * (only needed if tensorToReturn = 3)
     * 1: add epsInitial to deltaP/deltaE / add nuInitial to deltaM/deltaH
     * 2: add eps0 to deltaP/deltaE / add nu0 to deltaM/deltaH
     */
    UInt tensorToAddLHS;
    UInt tensorToAddRHS;
    
    /*
     * vectorToReturn
     * 1: put P(E_k) / M(H_k) on rhs
     * 2: put P=0 / M = 0 on rhs
     * 3: put P(E_k) - P(E_lastTS) on rhs / M(H_k) - M(H_lastTS)
     * 4: put only P(E_lastTS) / M(H_lastTS) on rhs
     */
    UInt vectorToReturn;
    
    /*
     * useLastTS
     * (only needed if tensorToReturn = 3)
     * true: deltaMat_k = (P(E_k) - P(E_lastTS))/(E_k - E_lastTS)
     *        or for magnetics with M and H
     *
     * false: deltaMat_k = (P(E_k) - P(E_k-1))/(E_k - E_k-1)
     *         or for magnetics with M and H
     */
    bool useLastTS;
    
    /*
     * evaluationPurpose
     * 1: assemble
     * 2,3: store
     * 4: output
     *
     * > here of course only option 1 is used
     * (for more details see coefFncHyst)
     */
    UInt evaluationPurpose = 1;
    PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
    
    /*
     * flag indicating if we need to reassemble system
     * (only needed if matrixes on rhs and lhs are different)
     */
    bool reassemble;
    
    /*
     * flag for easier switching between different rhs versions
     * 1: full (re-)evaluation
     *      electrostatics:
     *        res_k+1 = f_lin + f_nonlin( E_k+1 ) - div( eps0*E_k+1 ) - div( P_k+1 )
     *
     *      electromagnetics
     *        res_k+1 = J_lin + J_nonlin( B_k+1 ) - rot( nu0*B_k+1 ) + rot( M_k+1 )
     *
     * 2: incremental evaluation
     *      electrostatics:
     *        res_k+1 = res_k - f_nonlin( E_k ) + f_nonlin( E_k+1 ) - div( [eps0 + deltaP_k+1/deltaE_k+1]*deltaE_k+1 )
     *
     *      electromagnetics:
     *        res_k+1 = res_k - J_nonlin( B_k ) + J_nonlin( B_k+1 ) - rot( [nu0 + deltaM_k+1/deltaB_k+1]*deltaB_k+1 )
     */
    UInt REScase;
    
    /*
     * 0: RHS_k = RES_k
     * 1: further cases not needed at the moment, but have to be considered, when performing a full stepping
     *    i.e. we compute the increment towards the solution from the previous timestep instead of the value from the
     *    last iteration
     */
    UInt RHScase;
    
    
    /*
     * EvalVersion = 0
     *  > hystoperator does not affect solution process
     *  > solve linear system with eps0 / nu0
     *  > activate evaluation of hystoperator during output evaluation
     *
     *  for all iterations solve (rhs may have nonlinear, nonhyst parts):
     *    electrostatics:
     *
     *      div( eps0*E^n+1_k ) = f^n+1_k-1
     *
     *    electromagnetics
     *
     *      rot( nu0*B^n+1_k ) = J^n+1_k-1
     */
    if(evalVersion == 0){
      // DEBUG version; hysteresis not considered during solution process
      tensorToReturnLHS = 2; //always solve with eps0/nu0
      tensorToReturnRHS = 2; //eps0/nu0
      
      tensorToAddLHS = 1; //does not matter here
      tensorToAddRHS = tensorToAddLHS; //does not matter here but should be same as for LHS for reassemble flag to become false
      useLastTS = false; //does not matter here
      
      vectorToReturn = 2; //hystOperator will return 0 instead of P/M;
      REScase = 1; //incremental computation of res makes no sense as we have no deltaMat at all
      RHScase = 0; //RHS = RES
      
    }
    /*
     * EvalVersion = 1
     *  > during first iteration solve:
     *    electrostatics:
     *
     *      div( eps0*deltaE_1 ) = f^n+1_0 - div( eps0*E^n+1_0 ) - div( P^n+1_0 )
     *
     *        deltaE_1 = E^n+1_1 - E^n+1_0
     *
     *    electromagnetics
     *
     *      rot( nu0*deltaB_1 ) = J^n+1_0 - rot( nu0*B^n+1_0 ) + rot( M^n+1_0 )
     *
     *        deltaB_1 = B^n+1_1 - B^n+1_0
     *
     *  > during later iterations solve:
     *    electrostatics:
     *
     *      div( [eps0 + deltaP_k-1/deltaE_k-1]*deltaE_k ) = f^n+1_k-1 - div( eps0*E^n+1_k-1 ) - div( P^n+1_k-1 )
     *
     *        deltaE_k = E^n+1_k - E^n+1_k-1
     *
     *    electromagnetics
     *
     *      rot( [nu0 + deltaM_k-1/deltaB_k-1]*deltaB^k ) = J^n+1_k-1 - rot( nu0*B^n+1_k-1 ) + rot( M^n+1_k-1 )
     *
     *      deltaB_k = B^n+1_k - B^n+1_k-1
     *
     */
    else if(evalVersion == 1){
      // DELTA_MAT; rhs and lhs are affected by output of hysteresis operator
      //  during first iteration use no deltaMatrix
      if(callingCnt == 0){
        tensorToReturnLHS = 5; // 2 for eps0/nu0 // 1 to start with epsInitial instead // 4 start with estimated eps
        tensorToAddLHS = 2; // not needed here
      } else {
        tensorToReturnLHS = 3; // deltaMatrix
        tensorToAddLHS = 2; //add eps0 / nu0
      }
      
      useLastTS = false;
      if(forceReevaluation || (callingCnt == 0)){
        // full evaluation of residual; matrix has to be assemble with eps0/nu0; P/M has to
        // be put on rhs
        tensorToReturnRHS = 2; //eps0
        vectorToReturn = 1; //return current state of hystOperator (P/M)
        tensorToAddRHS = 1; //does not matter here
        REScase = 1;
      } else {
        tensorToReturnRHS = 3; //assemble deltaMat
        vectorToReturn = 2; //P is already included in deltaMat > return 0 instead
        tensorToAddRHS = 2; //assemble matrix with eps0 added to deltaP/deltaE / nu0 added deltaM/deltaB
        REScase = 2; //incremental residual
      }
      
      RHScase = 0; // RHS = RES
    }
    /*
     * EvalVersion = 2
     *  > during all iterations solve:
     *    electrostatics:
     *
     *      div( [eps0 + deltaP_k-1/deltaE_k-1]*deltaE_k ) = f^n+1_k-1 - div( eps0*E^n+1_k-1 ) - div( P^n+1_k-1 )
     *
     *        deltaE_k = E^n+1_k - E^n+1_k-1
     *        deltaE_0 = E^n+1_0 - E^n  > useLastTS = true
     *
     *    electromagnetics
     *
     *      rot( [nu0 + deltaM_k-1/deltaB_k-1]*deltaB^k ) = J^n+1_k-1 - rot( nu0*B^n+1_k-1 ) + rot( M^n+1_k-1 )
     *
     *      deltaB_k = B^n+1_k - B^n+1_k-1
     *      deltaB_0 = B^n+1_0 - B^n    > useLastTS = true
     *
     */
    else if(evalVersion == 2){
      // DELTA_MAT; rhs and lhs are affected by output of hysteresis operator
      tensorToReturnLHS = 3; // deltaMatrix
      tensorToAddLHS = 2; //add eps0 / nu0
      
      if(callingCnt == 0){
        useLastTS = true;
      } else {
        useLastTS = false;
      }
      
      if(forceReevaluation || (callingCnt == 0)){
        // full evaluation of residual; matrix has to be assemble with eps0/nu0; P/M has to
        // be put on rhs
        tensorToReturnRHS = 2; //eps0
        vectorToReturn = 1; //return current state of hystOperator (P/M)
        tensorToAddRHS = 1; //does not matter here
        REScase = 1;
      } else {
        tensorToReturnRHS = 3; //assemble deltaMat
        vectorToReturn = 2; //P is already included in deltaMat > return 0 instead
        tensorToAddRHS = 2; //assemble matrix with eps0 added to deltaP/deltaE / nu0 added deltaM/deltaB
        REScase = 2; //incremental residual
      }
      
      RHScase = 0; // RHS = RES
    }
    else{
      EXCEPTION("Evalversion > 3 not implemented yet");
    }
    
    if((tensorToReturnRHS == 3)&&(tensorToReturnLHS == 3)){
      /*
       * in case of delta mat on lhs and rhs, check also if the addTensor variable is the same
       */
      reassemble = (tensorToAddLHS != tensorToAddRHS);
    } else {
      reassemble = (tensorToReturnLHS != tensorToReturnRHS);
    }
    
    /*
     * #### PREPARATIONS ####
     */
    /*
     * set current solution by overwriting solVec_
     */
		
    //		std::cout << "oldSolution: " << oldSolution.ToString() << std::endl;
    //		std::cout << "solVec_.ToString: "  << solVec_.ToString() << std::endl;
		
    if(usedEta != 0){
      //      std::cout << "solVec before add " << solVec_.ToString() << std::endl;
      //      std::cout << "oldSolution before add " << oldSolution.ToString() << std::endl;
      //      std::cout << "solIncrement before add " << solIncrement.ToString() << std::endl;
      //      std::cout << "usedEta " << usedEta << std::endl;
      
      solVec_.Add( 1.0, oldSolution, usedEta, solIncrement);
      //      std::cout << "solVec after add " << solVec_.ToString() << std::endl;
      //      std::cout << "oldSolution after add " << oldSolution.ToString() << std::endl;
      //      std::cout << "solIncrement after add " << solIncrement.ToString() << std::endl;
      //      std::cout << "usedEta " << usedEta << std::endl;
      
    } else {
      solVec_ = oldSolution;
    }
    
    /*
     * we have to set the stageSol to the actual solutionVector;
     * this is needed for the handling of the time stepping scheme
     * (see further below)
     * by doing so, we can tell the ComputeStageRHS function to include
     * the current solution and thus consider an incremental stepping
     * 
     * Attention: this is only valid if we really want to subtract the whole
     * solution from the rhs! incremental computation of residual is not possible
     * thus!
     */
    stageSol = solVec_;
		
    SBM_Vector outputRHS(BaseMatrix::DOUBLE);
    if(debugOutput){
      algsys_->GetRHSVal( outputRHS );
      //      std::cout << "RHS before setting up residual: " << outputRHS.ToString() << std::endl;
      //
      //      std::cout << "old solution (oldSol): " << oldSolution.ToString() << std::endl;
      //      std::cout << "solInc (solInc): " << solIncrement.ToString() << std::endl;
      //      std::cout << "usedEta: " << usedEta << std::endl;
      //      std::cout << "current solution (oldSol + eta*solInc): " << solVec_.ToString() << std::endl;
      //
      //      std::cout << "Input parameter:" << std::endl;
      //      std::cout << "evalVersion: " << evalVersion << std::endl;
      //      std::cout << "forceReevaluation: " << forceReevaluation << std::endl;
      //      std::cout << "callingCnt: " << callingCnt << std::endl;
      //
      //      std::cout << "Derived flags:" << std::endl;
      //      std::cout << "REScase: " << REScase << std::endl;
      //
      //      std::cout << "tensorToReturnLHS: " << tensorToReturnLHS << std::endl;
      //      std::cout << "tensorToReturnRHS: " << tensorToReturnRHS << std::endl;
      //
      //      std::cout << "tensorToAddLHS: " << tensorToAddLHS << std::endl;
      //      std::cout << "tensorToAddRHS: " << tensorToAddRHS << std::endl;
      //
      //      std::cout << "vectorToReturn: " << vectorToReturn << std::endl;
      //      std::cout << "useLastTS: " << useLastTS << std::endl;
      //
      //      std::cout << "Reassembling needed: " << reassemble << std::endl;
      
    }
    
    /*
     * evaluate and store current linear RHS
     * > has to be done only once for each TS > during first call
     */
    if(callingCnt == 0){
      if(debugOutput){
        std::cout << "Assemble lin rhs" << std::endl;
      }
      currentLinRhsVec_.Init();
      algsys_->InitRHS(currentLinRhsVec_); //load storage into algsys
      assemble_->AssembleLinRHS(); //write current linear rhs to algsys
      algsys_->GetRHSVal( currentLinRhsVec_ ); //retrieve value
    }
    
    //				std::cout << "currentLinRhsVec_ " << currentLinRhsVec_.ToString() << std::endl;
    //		
		
    /*
     * evaluate and store current NON-linear RHS (without hysteresis!)
     * > has to be done during each iteration as it depends on the current solution
     */
    if(debugOutput){
      std::cout << "Assemble nonlin rhs (without hysteresis)" << std::endl;
    }
    currentNonLinRhsVec_.Init();
    algsys_->InitRHS(currentNonLinRhsVec_); //load storage into algsys
    PDE_.SetFlagInCoefFncHyst("vectorToReturn",2); // deactivate hysteresis!
    assemble_->AssembleNonLinRHS(); //write current linear rhs to algsys
    algsys_->GetRHSVal( currentNonLinRhsVec_ ); //retrieve value
    
    //		std::cout << "currentNonLinRhsVec_ " << currentNonLinRhsVec_.ToString() << std::endl;
    //		
		
    /*
     * #### EVALUATE RESIDUAL ACCORDING TO PREVIOUSLY SET FLAGS ####
     */
    /*
     * Assemble matrices
     */
    // setup flags for matrix computation
    PDE_.SetFlagInCoefFncHyst("tensorToReturn",tensorToReturnRHS);
    PDE_.SetFlagInCoefFncHyst("tensorToAdd",tensorToAddRHS);
    PDE_.SetFlagInCoefFncHyst("useLastTS",useLastTS);
    assemble_->AssembleMatrices(false);
    
    /*
     * Setup residual
     */
    currentResVec_.Init();
    
    /* REScase:
     * 1: full evaluation
     *        res_k+1 = f_lin + f_nonlin( E_k+1 ) - div( eps0*E_k+1 ) - div( P(E_k+1) )
     *
     * 2: incremental computation
     *        res_k+1 = res_k + f_nonlin( E_k+1 ) - f_nonlin( E_k ) - div( deltaMat_k*deltaE_k )
     *
     */
    //debugOutput = true;
    if(REScase == 1){
      if(debugOutput){
        std::cout << "RES - full reevaluation" << std::endl;
      }
      
      /* REScase:
       * 1: full evaluation
       *        res_k+1 = f_lin + f_nonlin( E_k+1 ) - div( P(E_k+1) ) - div( eps0*E_k+1 )
       */
      // add f_lin
      currentResVec_.Add(1.0, currentLinRhsVec_);
      
      // put res to rhs
      algsys_->InitRHS(currentResVec_);
      
      if(debugOutput){
        algsys_->GetRHSVal( outputRHS );
        std::cout << "RES - linear part: " << outputRHS.ToString() << std::endl;
      }
      
      // TODO clean up!
      // temporary fix
      // for full reevaluation, we have to evaluate the rhs loading with
      // the current state of the hyst operator
      PDE_.SetFlagInCoefFncHyst("forceCurrentTS",1); 

      
      // + f_nonlin( E_k+1 ) - div( P(E_k+1) )
      PDE_.SetFlagInCoefFncHyst("vectorToReturn",vectorToReturn); // activate hysteresis!
      assemble_->AssembleNonLinRHS(); //add current nonlinear rhs to algsys
      
      // reset
      PDE_.SetFlagInCoefFncHyst("forceCurrentTS",0); 
      
      if(debugOutput){
        algsys_->GetRHSVal( outputRHS );
        std::cout << "RES - linear+nonlinear+hyst part: " << outputRHS.ToString() << std::endl;
      }
			
      //			std::cout << "solVec_ "  << solVec_.ToString() << std::endl;
      //			std::cout << "solVec_[0].GetSize() " << solVec_(0).GetSize() << std::endl;
      //      algsys_->GetRHSVal( outputRHS );
      //			std::cout << "outputRHS "  << outputRHS.ToString() << std::endl;
      //			std::cout << "outputRHS[0].GetSize()" << outputRHS(0).GetSize() << std::endl;
      
      // - div( eps0*E_k+1 )
      solVec_.ScalarMult(-1.0);
			
			algsys_->UpdateRHS(STIFFNESS,solVec_,true);
			
			if(!trans){
				algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
			}
      //			
      //			if(!trans){
      //				algsys_->UpdateRHS(STIFFNESS,solVec_,true);
      //			} else {
      //				// for systems without mass/damping etc it was enough to just update
      //				// stiffness and stiffness update
      //				// howerver; for piezo, magnetic and magstrict, we have to update with
      //				// mass, damping etc, too
      //				// otherwise, the contributions from those matrices won't get to 0
      //			 // if(trans){
      //			 //   algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
      //			 // }
      //				// i.e. we have to update rhs with Keff*oldSolution or Keff*currentSolution
      //				/*
      //				 * However, the following lines do not work as intended, as the do not
      //				 *					compute Keff but Kmod = K*D*M instead
      //				 *				 (i.e. the time integration factors are missing in front of
      //				 *					D and M)
      //				 * However 2, this proceducre might not be needed at all as the
      //				 *				 there is the possibility to add the current solution
      //				 *				 directly to the computation of the stage RHS
      //				 *				 see TimeSchemeGLM.cc   void TimeSchemeGLM::ComputeStageRHS(UInt actStage, Integer derivId,
      //                                      SingleVector* rhsVec, Integer subIdx)
      //				 *				According to this function, the stage RHS will be considered
      //				 *				automatically if the nonlin stepping type is INCREMENTAL
      //				 * Open Question:
      //				 *				Does it only update with damping and mass or also with stiffness?
      //				 *
      //				 */
      //				
      //				std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
      //				std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
      //				std::map<FEMatrixType,Integer>::iterator matIt;
      //	//      UInt pos = 0;
      //
      //				// update RHS with all matrices!
      //				for(matIt = matrices.begin();matIt != matrices.end();matIt++){
      //					if(matIt->second < 0)
      //						continue;
      //
      //					algsys_->UpdateRHS(matIt->first,solVec_,true);
      //				}
      //			}
      
      solVec_.ScalarMult(-1.0);
      
      if(debugOutput){
        algsys_->GetRHSVal( outputRHS );
        std::cout << "RES - linear part+nonlinear+hyst part - solution: " << outputRHS.ToString() << std::endl;
        std::cout << "solution: " << solVec_.ToString() << std::endl;
        std::cout << "solIncrement: " << solIncrement.ToString() << std::endl;
      }
    } else if (REScase == 2){
      EXCEPTION("Incremental computation not fixed yet");
//			
//			if(debugOutput){
//        std::cout << "RES - incremental evaluation" << std::endl;
//        std::cout << "oldItResVec: " << oldItResVec_.ToString() << std::endl;
//      }
//      
//      /* REScase:
//       * 4: incremental computation
//       *        res_k+1 = res_k + f_nonlin( E_k+1 ) - f_nonlin( E_k ) - div( deltaMat_k*deltaE_k )
//       */
//      // add res_k
//      currentResVec_.Add(1.0, oldItResVec_);
//      
//      // + f_nonlin( E_k+1 ) - f_nonlin( E_k )
//      currentResVec_.Add(1.0, currentNonLinRhsVec_);
//      currentResVec_.Add(-1.0, oldItNonLinRhsVec_);
//      
//      // put res to rhs of system
//      algsys_->InitRHS(currentResVec_);
//      
//      //- div( deltaMat_k*deltaE_k )
//      solIncrement.ScalarMult(-1.0*usedEta);
//			
//			// test: update rhs with keff*solincrement
//			// however this might not be possible as we do not get the 
//			// correct scalings for damping and mass matrix 
//			// i.e. by the following, we apply
//			// (K+M*D)*solIncrement instead of
//			// (K+1/beta*dt^2 M + gamma/dt D)*solIncrement
//			// > here we might need an additional function that
//			// sclaes the solution vector by the correct entries
//			
//			
//			
//			if(!trans){
//				algsys_->UpdateRHS(STIFFNESS,solIncrement,true);
//			} else {
//				// for systems without mass/damping etc it was enough to just update
//				// stiffness and stiffness update
//				// howerver; for piezo, magnetic and magstrict, we have to update with
//				// mass, damping etc, too
//				// otherwise, the contributions from those matrices won't get to 0
//        // if(trans){
//        //   algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
//        // }
//				// i.e. we have to update rhs with Keff*oldSolution or Keff*currentSolution
//				std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
//				std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
//				std::map<FEMatrixType,Integer>::iterator matIt;
//        //      UInt pos = 0;
//        
//				// update RHS with all matrices!
//				for(matIt = matrices.begin();matIt != matrices.end();matIt++){
//					if(matIt->second < 0)
//						continue;
//          
//					algsys_->UpdateRHS(matIt->first,solIncrement,true);
//				}
//			}
//			
//      //			// old:
//      //      algsys_->UpdateRHS(STIFFNESS,solIncrement,true);
//      //
//      //      if(trans){
//      //        algsys_->UpdateRHS(STIFFNESS_UPDATE,solIncrement,true);
//      //      }
//      
//      solIncrement.ScalarMult(-1.0/usedEta);
//      
    } else {
      EXCEPTION("REScase unknown");
    }
    
    //retrieve value
    algsys_->GetRHSVal( currentResVec_ );
    
    if(debugOutput){
      algsys_->GetRHSVal( outputRHS );
      std::cout << "RES after setting up residual; before considering timestepping: " << outputRHS.ToString() << std::endl;
    }
    
    /*
     * common steps and all methods and callingCnt
     */
    /*
     * update with stage rhs / updating accoring to time stepping
     */
    if(trans){
      
      std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
      std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
      std::map<FEMatrixType,Integer>::iterator matIt;
      //      UInt pos = 0;
      
      // update RHS according to time stepping
			// in order to correctly subtract the current solution from the rhs
			// we have to make sure that the nonLinType is INCREMENTAL
			// If this is the case, ComputeStageRHS will also add the correctly
			// scaled stage solution to the stageRHS and thus we get
			// -(Keff-K)*stageSol on the rhs (note: K is usually not included in ComputeStageRHS
			// as the corresponding factor = 0 (for Newmark); so we have to subtract
			// K*stageSol separately (see above)
			bool forceIncremental = true;
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
        if(matIt->second < 0)
          continue;
        
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(fctId),-1,forceIncremental);
        }
        algsys_->UpdateRHS(matIt->first,stageRHS_,true);
      }
			
			// same issue maybe as with the initialisation > for coupled pde, pos does not give
			// the right subindex
      //      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
      //
      //        if(matIt->second < 0)
      //          continue;
      //        for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
      //          fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(pos));
      //        }
      //        algsys_->UpdateRHS(matIt->first,stageRHS_,true);
      //      }
    }
    
    if(debugOutput){
      algsys_->GetRHSVal( outputRHS );
      std::cout << "RES after setting up residual; after considering timestepping: " << outputRHS.ToString() << std::endl;
    }
    debugOutput = false;
    /*
     * update current residual with information of time stepping
     * TODO: not actually sure about this step; should the time stepping update count to residual or not?
     *       pro:
     *        via this step we take mass matrix * second time derivative into account
     *        including these values will actually lead to decreasing residual which is not the case without the addition
     *
     */
    algsys_->GetRHSVal( currentResVec_ );
    
      //  std::cout << "reassemble: " << reassemble << std::endl;
      //  std::cout << "skipReassembly: " << skipReassembly << std::endl;
    /*
     * #### REASSEMBLE SYSTEM ####
     */
    if((reassemble)&&(skipReassembly == false)){
      
      if(debugOutput){
        std::cout << "Re-Assembling system for next iteration" << std::endl;
      }
      
      //        std::cout << "tensorToReturnLHS: " << tensorToReturnLHS << std::endl;
      PDE_.SetFlagInCoefFncHyst("tensorToReturn",tensorToReturnLHS);
      PDE_.SetFlagInCoefFncHyst("tensorToAdd",tensorToAddLHS);
      PDE_.SetFlagInCoefFncHyst("useLastTS",useLastTS);
      assemble_->AssembleMatrices(false);
    }
    					
				
    /*
     * #### REPLACE RHS (full stepping only) ####
     */
    if(RHScase != 0){
      EXCEPTION("Currently only incremental stepping implemented");
    }
    
  }
    
  Double SolveStepHyst::LineSearchHyst(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double& etaLineSearch, UInt evalVersion,
          UInt callingCnt, bool trans, bool performLineSearch,
          bool forceReevaluation, bool debugOutput, bool reset, UInt allowedSteps)  {
    
    static Timer* lineSearchTimer = new Timer();
    static UInt numCalls = 0;
    if(debugOutput){
      std::cout << "LineSearchHyst" << std::endl;
      lineSearchTimer->Start();
      numCalls++;
    }
    
    /*
     * backup nonlinear rhs, residual and solution vector
     */
    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = solVec_;
    
    SBM_Vector resOld(BaseMatrix::DOUBLE);
    resOld = currentResVec_;
    
    
    const UInt nrEtas = 13;
    //Double eta[nrEtas] = {1.0,0.5,0.2,0.1};
    Double eta[nrEtas] = {1.0,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.0025,0.001,0.0005,0.00025,0.0001};
    //Double eta[nrEtas] = {1.0,0.316,0.1,0.0316,0.01,0.00316,0.001,0.000316,0.0001,0.0000316,0.00001,0.0000316,0.00001};
    //Double eta[nrEtas] = {1.0,0.075,0.05,0.025,0.1};
    //Double eta[nrEtas] = {1.0,0.316,0.1,0.0316,0.01};
    //Double eta[nrEtas] = {1.0,0.75,0.5,0.25,0.1};
    bool assumeQuadratic = false;
    
    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;
    
    if(performLineSearch == false){
      /*
       * only take first eta = 1
       * -> we do not have to lock hysteresis nor do we have to backup
       *    the solution or the residual as we directly take the first value
       */
      etaOpt = 1.0;
    } else {
      /*
       * lock hysteresis memory
       * NOTE: in the current evaluation, the locking is done automatically
       * in coefFnc hyst if the system is called during assembly;
       * only if SetPreviousHystValues is used with setLastTS = true, the
       * hysteresis memory is unlocked
       */
      
      /*
       * test different eta
       */
      bool skipReassembly = true;
      UInt nrEtasToTest = std::min(allowedSteps,nrEtas);
      for( UInt i=0; i<nrEtasToTest; i++) {
        if(debugOutput){
          std::cout << "testing eta: " << eta[i] << std::endl;
        }
        
        CalcResidualAndConfigSystemForHysteresis(solOld,solIncrement,stageSol,eta[i], 0, callingCnt, evalVersion,
                trans, forceReevaluation, skipReassembly, debugOutput, reset);
        
        /*
         * calculate error norm as criterion
         */
        Double residualL2Norm = currentResVec_.NormL2();
        
        if(debugOutput){
          std::cout << "L2-norm of computed residual for eta = " << eta[i] << ": " << residualL2Norm << std::endl;
        }
        
        /*
         * restore old solution vector, old residual and nonlinear rhs
         * -> needed to perform a correct next call to CalcResidual...
         */
        solVec_ = solOld;
        currentResVec_ = resOld;
        
        if (residualL2Norm < residualL2NormOpt) {
          residualL2NormOpt = residualL2Norm;
          etaOpt = eta[i];
        } else {
          if(assumeQuadratic){
            /*
             * here we assume that the residual follows a quadratic curve
             * i.e. we only have one minimum which is the global minimum;
             * in that case, if the residual increases between two iterations,
             * we can assume that we passed the minimum already and stop
             */
            break;
          }
        }
      }
      
      /*
       * unlock hysteresis memory
       * -> no longer done here but first when we store the lastTS values
       */
      
    }
    
    /*
     * now we have found etaOpt, so we can call CalcResidualAndConfigSystemForHysteresis
     * one last time to setup the correct system for next computation
     */
    etaLineSearch = etaOpt;
    
    if(debugOutput){
      lineSearchTimer->Stop();
      std::cout << "Total time needed for linesearch: " << lineSearchTimer->GetCPUTime() << std::endl;
      std::cout << "Average time needed for linesearch: " << lineSearchTimer->GetCPUTime()/numCalls << std::endl;
      std::cout << "(optimal) eta: " << etaOpt << std::endl;
    }
    
    // calculate res and assemble system once more using the optimal eta
    bool skipReassembly = false;
    CalcResidualAndConfigSystemForHysteresis(solOld,solIncrement,stageSol,etaLineSearch, 0, callingCnt, evalVersion,
            trans, forceReevaluation, skipReassembly, debugOutput, reset);
    
    if(debugOutput){
      std::cout << "minimal residual: " << residualL2NormOpt << std::endl;
    }
    
    return currentResVec_.NormL2();
  }
  
  void SolveStepHyst::StepTransNonLinHysteresis() {
    /*
     * New approach / new hope?
     *
     * Idea:
     *  1. Iteration: full step starting at 0
     *    i.e. Reset solution vector to 0
     *         Reset old iteration vector to 0
     *         Lock memory of hysteresis operator
     *          > it will act as a nonlinear curve; no branching
     *          > due to reset of solution to 0 it will be dropped from its last state
     *            (= solution state after last timestep) to a remanent state
     *            > this drop is not stored in the memory of the operator, i.e.
     *              if we would apply the solution of the last ts to the operator, we
     *              would end up at the same state as prior to the drop to remanence
     *         Assemble system using the difference between the values from lastTS
     *         and the reset solution
     *  2-n. Iterartion: incremental stepping
     *    i.e. Still lock hyst operator
     *         Assemble system using difference from current state to the lastIT state
     * */
    static UInt timestepCnt = 0;
    timestepCnt++;
    if(debugOutput_){
      std::cout << "####                   #### " << std::endl;
      std::cout << "#### BEGIN OF TIMESTEP #### " << std::endl;
      std::cout << "####   Nr: " << timestepCnt << "   #### " << std::endl;
    }
    
    LOG_TRACE(solvestephyst) << "solvestephyst::StepTransNonLinHysteresis";
    
    SBM_Vector solInc(BaseMatrix::DOUBLE);
    SBM_Vector actSol(BaseMatrix::DOUBLE);
    
    //		std::cout << "solVec_: " << solVec_.ToString() << std::endl;
    //		
    bool deltaIDBC = true;
    actSol = solVec_;
		
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    if ( numStages > 1 ){
      /*
       * maybe it is quite easy to enable multiple stages, maybe it does not work;
       * feel free to find out
       */
      EXCEPTION("StepTransNonLinHysteresis: just one stage time-stepping allowed");
    }
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;
    
    UInt pos = 0;
    
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }
    
    //do initialization
    rhsVec_.Init();
    algsys_->InitRHS(rhsVec_);
    
    LOG_DBG(solvestephyst) << "StepTransNonLinHysteresis: Stage: " << 0 ;
    
    //we obtain a reference to the stage vectors of the scheme
    /*
     * Important: Although stageSol seems not to be necessary, the TimeScheme will rely on it
     * It HAS to store the total solution (not the increment!) -> see TimeSchemeGLM.cc
     * Important2: Incremental stepping is only working (correctly) for EffectiveStiffness formulation
     */
    SBM_Vector stageSol;
    stageSol.Resize(feFunctions_.size());
		
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      FeFctIdType fctId = fncIt->second->GetFctId();
      stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(0),fctId);
      fncIt->second->GetTimeScheme()->InitStage(0,actTime_,PDE_.GetDomain());
    }
		
		// so nicht! dadurch ist die stageSol womöglich falsch aufgebaut (z.b. mech und
		// elec unknows vertauscht!
    //    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
    //      stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(0),pos);
    //      fncIt->second->GetTimeScheme()->InitStage(0,actTime_,PDE_.GetDomain());
    //    }
    stageSol.SetOwnership(false);
		
    /*
     * During the later computations, we will add all incremental solutions to
     * stageSol; unlike the case of StepTransNonLin, we add the results to the
     * current solution instead of only gathering the increments.
     * Note that we initialize stageSol with actSol in StepTransNonLin, too, but
     * there, we will Init stageSol with 0.0 during the first iteration
     */
    stageSol = actSol;
    solVec_  = actSol;
    
		/*
     Irgendetwas stimmt hier doch nicht!
     Am Anfang hat solVec_ die richtige Größe
     Nachdem actSol an stageSol zugewiesen wurde, hat solVec plötzlich 
     subVector 1 und 2 vertauscht, stageSol und actSol sind aber noch richtig
     
     Danach nachdem solVec = actSol passen die Größen von solVec und actSol
     wieder, aber stageSol hat jetzt vertauschte einträge
     * 
     * Setzt man statt solVec = actSol, solVec = stageSol so passt gar nichts
     * mehr
     */
		
    // setup right hand side
    Double loadFactor = 1.0;
    Double RhsLinL2Norm = SetLinRHS(loadFactor);
    
    /*
     * disable output of switching and rotation state as BMP images (even if flag
     * printOut is set to 1 in material file)
     */
    PDE_.SetFlagInCoefFncHyst("allowBMP",false);
    
    /*
     * START ITERATIONS
     */
    static UInt overallIterationCounter = 0;
    bool performOneMoreStep;
    bool transient = true;
    UInt iterationCounter = 0;
    UInt stage = 0;
    Double initialEta = 0.0;
    Double residualErr;
    Double incrementalErr;
    bool forceReevaluation = forceReevaluation_;
    
    /*
     * NOTE:
     * evalVersion_
     * 0 = debug fixpoint > hysteresis will be evaluated only during output
     *
     * 1 = fixpoint / deltaMat
     *       lhs and rhs with eps0/nu0 > only during first iteration!
     *       lhs = deltaMat, rhs with eps0/nu0 during later iteration
     * 2 = deltaMat during all iterations
     */
    //    std::cout << "evalVerion: " << evalVersion_ << std::endl;
    
    /*
     * ITERATION 0
     * > Prepare initial system
     */
    PDE_.SetFlagInCoefFncHyst("estimateSlope",9);
    bool skipReassembly = false; // currently not used
    
    bool reset = false; 
    UInt numResets = 0;
    
    // keep track of residual; if it is not decreasing over several iterations > reset (see below)
    std::deque<Double> resNormHistory(5, 50000);
//    
//    std::cout << "solVec_ at start of TS: " << solVec_.ToString() << std::endl;
//    std::cout << "stageSol at start of TS: " << stageSol.ToString() << std::endl;
//    std::cout << "actSol at start of TS: " << actSol.ToString() << std::endl;
//    
    CalcResidualAndConfigSystemForHysteresis(solVec_,solInc,stageSol, initialEta, stage,
            iterationCounter, evalVersion_, transient,
            forceReevaluation, skipReassembly, debugOutput_, reset);
    // after this initial setup, we have to call SetLastItOrLastTSSBMVectors, to store the current
    // res and nonlinRHS as the values from the previous iteration so that these values are correct
    // when we call CalcResidualAndConfigSystemForHysteresis during the first actual iteration
    SetLastItOrLastTSSBMVectors(false);
    /*
     * ITERATION > 0
     */
    do {
      
      overallIterationCounter++;
      iterationCounter++;
      
      LOG_DBG(solvestephyst) << "StepTransNonLinHysteresis: Iteration: " << iterationCounter ;
      
      if(debugOutput_){
        std::cout << "####                   #### " << std::endl;
        std::cout << "#### BEGIN OF ITERATION #### " << std::endl;
        std::cout << "####   Nr: " << iterationCounter << "   #### " << std::endl;
        std::cout << "####   Total Nr: " << overallIterationCounter << "   #### " << std::endl;
      }
      
      /*
       *  Assembly has already be done during CalcResidualAnd...
       *  Now build up the algebraic system
       */
      // set system matrix to zero initially, as ConstructEffectiveMatrix only
      // sums up the contributions
      matrix_factor_.clear();
      algsys_->InitMatrix(SYSTEM);
      for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
        FeFctIdType fctId = fncIt->second->GetFctId();
        fncIt->second->GetTimeScheme()->AddMatFactors(0,matrices,matrix_factor_[fctId]);
        algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
      }
      
      PDE_.SetBCs();
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
      
      bool setIDBC = false;
      
      if (( iterationCounter == 1 && couplingIter_ == 0 )||(reset)) {
        /*
         * we only need to include the idbc values once during the first iteration
         * in the new approach, we do not need deltaIDBC as we start from 0
         */
        setIDBC = true;
      }
      
      /*
       *  SOLVE SYSTEM
       */
      algsys_->Solve(setIDBC,deltaIDBC);
      
      /*
       *  RETRIEVE SOLUTION INCREMENT
       *  1, iteration: increment towards 0, i.e. full solution
       *  2+ iteration: increment towards previous iteration
       */
      algsys_->GetSolutionVal(solInc, setIDBC, deltaIDBC );
      
      /*
       * CALCULATE RESIDUAL AND SETUP SYSTEM FOR NEXT ITERATION VIA LINESEARCH
       */
      Double residualL2Norm = 0.0;
      Double etaLineSearch  = 1.0;
      bool performLinesearch;
      
      if (( lineSearch_ == "none" || iterationCounter == 1 )||(reset)) {
        /*
         * during the first step, we need a full step with eta = 1.0
         * this has to be done to ensure correct idbc values!
         */
        performLinesearch = false;
      } else {
        performLinesearch = true;
      }
      
      if(reset){
        // reset reset flag
        reset = false;
      }
      
      // the linesearching will test different etas and take the one which
      // produces the smallest residual
      // afterwards the system will be setup (assembled) with this new solution
      // new solution (oldSolution + etaOpt*solutionIncrement) is written to solVec_
      // with each reset, try also one smaller step for linesearch
      UInt allowedSteps = 7+2*numResets;
      residualL2Norm = LineSearchHyst(solInc, stageSol, etaLineSearch, evalVersion_, iterationCounter,
              transient , performLinesearch ,forceReevaluation, debugOutput_, reset, allowedSteps);
      
      // calculation of relative residual error =======================================
      //if ( RhsLinL2Norm > 0.0 ){
      //Double oldResError = residualErr;
      if ( RhsLinL2Norm > 1.0 ){
        residualErr = residualL2Norm / RhsLinL2Norm;
      } else {
        residualErr = residualL2Norm;
      }
      
      resNormHistory.push_back(residualErr);
      resNormHistory.pop_front();
      
      // calculate incremental error ========================================
      //TODO: we should have eta*solInc.NormL2() as we do not perform full steps in general
      //Double solIncrL2Norm = solInc.NormL2();
      
      Double solIncrL2Norm = etaLineSearch*solInc.NormL2();
      Double solL2Norm  = solVec_.NormL2();
      
      if ( solL2Norm > 1.0 ){
        // if ( actSolL2Norm > 1e-7 ){
        incrementalErr = solIncrL2Norm / solL2Norm;
      } else {
        incrementalErr = solIncrL2Norm;
        // incrementalErr = solIncrL2Norm/1e-7;
      }
      
      // update actSol to currently computed value
      actSol = solVec_;
      
      performOneMoreStep =
              (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
      
      if(residualErr < 1e-18){
        /*
         * in some cases e.g. E -> 0, P != 0 (remanence)
         * the residual might become very small although the solution does not converge
         * if the residual becomes too small ( ~ 1e-28) the system reports
         * error during solution (at least if using pardiso);
         * -> a complete reevaluation of the residual might help
         *
         */
        forceReevaluation = true;
      } else {
        /*
         * restore old state
         */
        forceReevaluation = forceReevaluation_;
      }
      
      if(debugOutput_){
        std::cout << "UsedEta: " << etaLineSearch << std::endl;
        std::cout << "Solution update: " << solInc.ToString() << std::endl;
        std::cout << "current solution norm: " << solL2Norm << std::endl;
        std::cout << "incrementalErr_abs: " << solIncrL2Norm << std::endl;
        std::cout << "incrementalErr_rel: " << incrementalErr << std::endl;
        std::cout << "residualErr: " << residualErr << std::endl;
      }
      
      // compare residual at front and back; should have decreased; otherwise try reset
      if(resNormHistory.front() <= resNormHistory.back()){
        std::cout << "Residual did not decrease. Reset?" << std::endl;
        std::cout << "Iteration " << iterationCounter << " of timestep " << timestepCnt << std::endl;
        std::cout << "Number of previous resets: " << numResets << std::endl;
        reset = true;
        
        // reset deque
        for(UInt i = 0; i < resNormHistory.size(); i++){
          resNormHistory[i] = 50000;
        }
      }
      
      if(reset == true){
        numResets += 1;
        //      std::cout << "Perform reset" << std::endl;
        // set solution to 0 i.e. set new starting point for iterations
        // note: the memory of the hyst operator is not deleted but as we do not
        //       write it until the end of the timestep, it is basically still the
        //       state from the start of the timestep
        solVec_.Init();
        // estimate slope around this point
        // estimate using a stepping
        PDE_.SetFlagInCoefFncHyst("estimateSlope",9+numResets);
        // for the first iteration after the reset, we have to include the boundary conditions
        // again and they must not e the deltaBCs but the full boundary conditions (as we start from 0)
        deltaIDBC = false;
        
        // now we have to reassemble the system (we basically repeat all steps till start of loop)
        // the reset flag will have the same meaning as an iterationCounter of 0 (i.e. full evaluation, computation of rhs, ...)
        CalcResidualAndConfigSystemForHysteresis(solVec_,solInc,stageSol, initialEta, stage,
                iterationCounter, evalVersion_, transient,
                forceReevaluation, skipReassembly, debugOutput_, reset);
        
        if(numResets > 4){
          EXCEPTION("Error cannot be decreases by resetting");
        }                                                 
      }
      
      
      // output of norms and data to info.xml
      if ( nonLinLogging_ == true ) {
        // get current step
        UInt actStep = PDE_.GetSolveStep()->GetActStep();
        
        if (PDE_.IsIterCoupled()) {
          WriteNonLinIterToInfoXML(pdename_, couplingIter_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
        } else {
          //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
          WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
        }
        // write norm to file
        logFile_ <<  iterationCounter << "\t"
                << residualErr << "\t"
                << incrementalErr << "\t"
                << etaLineSearch << std::endl;
      }
      
      /*
       * STORE CURRENT VALUES FOR NEXT ITERATION
       */
      bool lastTS = false;
      // Note: by calling SetPreviousHystVals with flag lastTS = false,
      //        coefFncHyst will evaluate the hysteresis operator with
      //        the current state of solVec_; for this evaluation, the memory of
      //        the hysteresis operator will be locked, i.e. the evaluation will
      //        lead to reversible changes
      PDE_.SetPreviousHystVals(lastTS);
      SetLastItOrLastTSSBMVectors(lastTS);
      
    } while(performOneMoreStep && (iterationCounter < nonLinMaxIter_));
    
    abortOnMaxIter_ = true;
    if (performOneMoreStep && (iterationCounter >= nonLinMaxIter_) && abortOnMaxIter_) {
      EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
              << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
              << "' at iteration '" << iterationCounter
              << "'.\n ==> incremental error: " << incrementalErr
              << "\n ==> residual error: " << residualErr);
    }
    
    // check solution
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator limitFeFctIt;
    limitFeFctIt = feFunctions_.find(solutionLimit_);
    if (limitFeFctIt != feFunctions_.end() ) {
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
          //const SingleVector * subv = solVec_.GetPointer(pos);
          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(pos)));
          for (UInt j=0; j < dsubVec.GetSize(); j++) {
            if (dsubVec[j] >= maxValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
                      "' is larger than the allowed maximum limit set in the XML: "
                      << maxValidValue_);
            }
            if (dsubVec[j] <= minValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
                      "' is smaller than the allowed minimum limit set in the XML: "
                      << minValidValue_);
            }
          }
        }
      }
    }
    
    //update stage
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }
    
    if(debugOutput_){
      std::cout << "Iterative scheme was successful after " << iterationCounter << " iterations" << std::endl;
    }
    
    /*
     * STORE CURRENT VALUES FOR NEXT ITERATION
     */
    bool lastTS = true;
    // Note: by calling SetPreviousHystVals with flag lastTS = true,
    //        coefFncHyst will evaluate the hysteresis operator with
    //        the final state of solVec_; for this evaluation, the memory of
    //        the hysteresis operator will be unlocked, i.e. the evaluation will
    //        lead to irreversible changes
    PDE_.SetPreviousHystVals(lastTS);
    SetLastItOrLastTSSBMVectors(lastTS);
    
    /*
     * Store IDBC values from the current time step in form of its rhs representation
     * (i.e. the effect that it will have on the rhs)
     * these values are needed to compute the deltaIDBC values
     * (currently not used)
     */
    algsys_->SetOldDirichletValues();
    
    /*
     * allow to output switching and rotation state as BMP images (if flag
     * printOut is set to 1 in material file)
     */
    PDE_.SetFlagInCoefFncHyst("allowBMP",true);
    
    /*
     * set evaluationPurpose to 4, i.e. output
     * this will lock the hysteresis operator again and will only evaluate
     * the hystOperator at the midpoint of each element regardless of the
     * evaluation depth
     * > for more infos see coefFncHyst
     */
    PDE_.SetFlagInCoefFncHyst("evaluationPurpose",4);
    
    if(debugOutput_){
      std::cout << "####                 #### " << std::endl;
      std::cout << "#### END OF TIMESTEP #### " << std::endl;
      std::cout << "####                 #### " << std::endl;
      PDE_.SetFlagInCoefFncHyst("outputDebugInfos",1);
    }
    
    //TODO: continue here
    //TODO: add defaultValue to xmlSchema!!!!!
  }
  
 
  void SolveStepHyst::StepTransNonLinHysteresisTotal() {
    
    bool performOneMoreStep;
    bool isNewton;
    Double incrementalErr;
    
    SBM_Vector newSol(BaseMatrix::DOUBLE);
    SBM_Vector oldSol(BaseMatrix::DOUBLE);
    SBM_Vector diffSol(BaseMatrix::DOUBLE);
    
    //get solution of previous time step
    SBM_Vector  prevSol(BaseMatrix::DOUBLE);
    prevSol = solVec_;
    oldSol = solVec_;
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;
    
    UInt pos = 0;
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }
    
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    if ( numStages > 1 )
      EXCEPTION("StepTransNonLinHysteresis: just one stage time-stepping allowed");
    
    for(UInt i=0;i<numStages;i++){
      //do initialization
      rhsVec_.Init();
      
      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),pos);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);
      
      //initialize solution vector for each stage
      if ( i > 0 )
        oldSol = stageSol;
      else{
        //special case of incremental non-linearity, we set the stage vector to the solution vector
        stageSol = oldSol;
      }
      
      solVec_  = oldSol;
      
      // setup right hand side
      Double loadFactor = 1.0;
      
      incrementalErr = SetLinRHS(loadFactor);
      
      // set iteration counter
      UInt iterationCounter=0;
      
      do {
        iterationCounter++;
        
        //std::cout << "Current iteration: " << iterationCounter << std::endl;
        
        // reset rhs
        //  RhsLinVal_.Init();
        //  algsys_->InitRHS(RhsLinVal_);
        //
        if(PDE_.IsHysteresis_Fixpoint() == true){
          
          // nach dem löschen muss er eigentlich sowohl die linearen als auch die nichtlinearen teile nochmals hinzufügen
          // die unterscheidung zwischen linear und nichtlinear ist also obsolet
          
          RhsLinVal_.Init();
          algsys_->InitRHS(RhsLinVal_);
          
          //     std::cout << "SetLinRHS" << std::endl;
          incrementalErr = SetLinRHS(loadFactor,true);
          //     std::cout << "Done with setting RHS" << std::endl;
          /*
           * Fixpoint iteration used to work, however, P was just postprocessed upon E
           * as we RhsLinVal_.Init() followed by SetLinRHS( ,true) we just had nonlinear terms on the rhs
           * although P depends on E it was not marked as solDependent and such we solved div(eps0*E) = 0
           * this worked of course
           * Doing a real updating on the rhs will however not work
           */
          //incrementalErr = SetLinRHS(loadFactor,false);
        } else {
          //incrementalErr = SetLinRHS(loadFactor,false);
        }
        
        // do matrices: Newton is not working for total formulation!!
        isNewton = false;
        assemble_->AssembleMatrices(isNewton);
        
        // set RHS
        algsys_->InitRHS(RhsLinVal_);
        
        //now update RHS according to time stepping
        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
          if(matIt->second < 0)
            continue;
          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
            fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
          }
          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
        }
        
        if(PDE_.IsHysteresis_Fixpoint() == false){
          //due to incremental material formulation
          //std::cout << "RHS updated with old solution" << std::endl;
          algsys_->UpdateRHS(STIFFNESS,prevSol,true);
        }
        // set system matrix to zero initially, as ConstructEffectiveMatrix only
        // sums up the contributions
        matrix_factor_.clear();
        algsys_->InitMatrix(SYSTEM);
        
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()->AddMatFactors(i,matrices,matrix_factor_[fctId]);
          algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
        }
        
        PDE_.SetBCs();
        algsys_->BuildInDirichlet();
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        
        //always set inhomog. Dirichlet BCs
        bool setIDBC = true;
        
        algsys_->Solve(setIDBC);
        algsys_->GetSolutionVal(newSol, setIDBC );
        
        // calculate incremental error ========================================
        diffSol = newSol;
        diffSol.Add( -1.0, oldSol);
        Double solIncrL2Norm = diffSol.NormL2();
        Double solNewL2Norm = newSol.NormL2();
        
        if (solNewL2Norm > 1)
          incrementalErr = solIncrL2Norm / solNewL2Norm;
        else
          incrementalErr = solIncrL2Norm;
        
        //just dummy things
        Double etaLineSearch = 1.0;
        Double residualErr = 1.0;
        // why not just try the linesearch algorithm?
        // actual linesearch is meant for the case of delta formulation
        //   i.e. we compute deltaSol using the equation system
        //   K(oldSol) * deltaSol = f - K(oldSol)*oldSol
        //   then we find out eta such that residual
        //   f - K(oldSol + eta*deltaSol)*(oldSol + eta*deltaSol) is minimal
        //   finally newSol = oldSol + eta*deltaSol
        // in our case, we compute newSol directly via
        //   K(oldSol)*newSol = f
        //   we can however try the following here
        //     compute diffSol = newSol - oldSol
        //     use linesearch to find a better eta than 1.0
        //     set newSol = oldSol + eta*diffSol
        // Note: oldSol will be overwritten!
        stageSol = oldSol;
        
        // before performing linesearch, make hysteresis operator read only,
        // i.e. do not store the effects of temporal input vectors (as they are
        // used by the linsearch function to test which eta is the best)
        //PDE_.LockHysteresis();
        
        residualErr = LineSearch(diffSol,stageSol,etaLineSearch,true);
        //    std::cout << "Used etaLineSearch: " << etaLineSearch << std::endl;
        //    std::cout << "Old solution: " << oldSol.ToString() << std::endl;
        //    std::cout << "New solution: " << stageSol.ToString() << std::endl;
        
        //PDE_.UnlockHysteresis();
        // Why is residualErr = incrementalErr?
        //Double residualErr = incrementalErr;
        
        // calculate new incremental error ========================================
        diffSol = stageSol; // this is the solution we just calculated via linesearch
        diffSol.Add( -1.0, oldSol);
        solIncrL2Norm = diffSol.NormL2();
        solNewL2Norm = newSol.NormL2();
        
        if (solNewL2Norm > 1)
          incrementalErr = solIncrL2Norm / solNewL2Norm;
        else
          incrementalErr = solIncrL2Norm;
        
        //      std::cout << "Incremental error: " << std::endl;
        //     std::cout << "rel error: " << incrementalErr << std::endl;
        
        //    std::cout << "Residual error: " << std::endl;
        //    std::cout << "rel error: " << residualErr << std::endl;
        
        // output of norms and data
        if ( nonLinLogging_ == true ) {
          WriteNonLinIterToInfoXML(pdename_, 1,iterationCounter, residualErr, incrementalErr, etaLineSearch);
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
                  << residualErr << "\t"
                  << incrementalErr << "\t"
                  << etaLineSearch << std::endl;
        }
        
        //use relaxated form
        //stageSol = newSol;
        //Double relaxfac = 0.1;
        //stageSol.Add( -relaxfac, newSol);
        //stageSol.Add( relaxfac, oldSol);
        //stageSol = newSol;
        solVec_  = stageSol;
        
        //store new solution for next iteration
        oldSol = stageSol; //newSol;
        
        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
                (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
        
        // std::cout << "IncError: " << incrementalErr << std::endl;
        // std::cout << "ResError: " << residualErr << std::endl;
        
        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
                  << "' at iteration '" << iterationCounter
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }
        
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
      
    } //stages
    
    //update stage
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }
    
    PDE_.FinalizeAfterTimeStep();
  }
  
} // end of namespace
