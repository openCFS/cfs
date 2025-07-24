

#pragma once

#include <list>

#include <boost/utility.hpp>

#include "MatVec/Vector.hh"
#include "Domain/Domain.hh"
#include "Model.hh"
#include "SMSM.hh"


namespace CoupledField {

  class MathParser;

  class EBHysteresis : public Model {

    public:
      //! Constructor
      EBHysteresis();

      //! Destructor
      virtual ~EBHysteresis();

      void Init(std::map<std::string, double> ParameterMap, std::map<std::string, string> StringParameterMap, shared_ptr<ElemList> entityList, UInt dim);

      Double ComputeMaterialParameter(Vector<Double> E, Integer ElemNum);
      Matrix<Double> ComputeTensorialMaterialParameter(Vector<Double> E, Integer ElemNum);

      // This shall act as a unified handling for updating states from one timestep to another, triggered by the SolveStepEB.
      // Originally this was handled via flags but manually updating the states (via an explicit call in SolveStepEB) is way 
      // safer and probably also faster since we do not need to check the states for every element at every iteration.
      void UpdateStates();
      void AllowUpdates(bool allow);

      // just for the computation of the residual, we do not store anything here
      Vector<Double> GetFluxDensity(Vector<Double> E, Integer ElemNum);

      Vector<Double> GetFluxDensity(Vector<Double> E, Integer ElemNum,
                                    LocPointMapped lpm, PtrCoefFct stressCoef);

      /*
      =========================================================================================
          Methods regarding the inverse scheme
      =========================================================================================
      */
      Vector<Double> ComputeMaterialParamaterDerivative(Vector<Double> E, UInt idx);

      /*
      =========================================================================================
          Methods regarding the computation of the Jacobian matrix
      =========================================================================================
      */
      Matrix<Double> EvaluateLocalMu(StdVector<Double> E, StdVector<Double> D, UInt idx);

      Matrix<Double> EvaluateLocalMuDFP(StdVector<Double> E, StdVector<Double> D, UInt idx);

      Matrix<Double> EvaluateLocalMuGBM(StdVector<Double> E, StdVector<Double> D, UInt idx);

      Matrix<Double> EvaluateLocalMuBFGS(StdVector<Double> E, StdVector<Double> D, UInt idx);

      Matrix<Double> EvaluateLocalMuFiniteDifferences(Vector<Double> E, StdVector<Double> D, UInt idx);

      Matrix<Double> EvaluateLocalMuAnhystersisOnly(Vector<Double> E, UInt idx);

      StdVector<Double> inv3x3(StdVector<Double> A);

      /*
      =========================================================================================
          Methods regarding the evaluation of the hysteresis model
      =========================================================================================
      */
      Vector<Double> Evaluate(Vector<Double> E, UInt idx);

      Vector<Double> Eval_2D_EBM_ATAN(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);
      Vector<Double> Eval_2D_VPM_ATAN(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);
      Vector<Double> Eval_3D_VPM_ATAN(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      Vector<Double> Eval_2D_VPM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);
      Vector<Double> Eval_3D_VPM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      Vector<Double> Eval_2D_EBM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);
      Vector<Double> Eval_3D_EBM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      /*
      =========================================================================================
          Methods regarding the computation of the simplified multiscale model (SMSM)
      =========================================================================================
      */
      void Calc_derivs(Double &F_prime, Double &F_prime_prime, Matrix<Double> &dM_dHrev,
                      Double MxS_prev, Double MyS_prev,
                      Double Hex_x, Double Hex_y,
                      Double HrxS_prev, Double HryS_prev,
                      Double phi, Double chi);

      Double Energy_linesearch(Double Hx, Double Hy, Double Hprev_x, Double Hprev_y, Double Mprev_x, Double Mprev_y,
                               Double phi, Double chi, Double &F_prime_orig, Double &F_prime_prime_orig);
    private:
      //==============
      std::unique_ptr<SMSM> SMSM_model_;

      UInt dim_; 

      std::map<UInt,UInt> ElemNum2Idx_;

      UInt numElems_;

      Double Ps_;
      Double A_;
      Double mu0_;
      UInt numS_;
      Double chi_factor_;
      UInt jacobian_method_;
      UInt anhyst_type_;
      std::string approx_type_;
      std::string pinning_forces_weight_;
      StdVector<Double> kappa_file_, omega_file_;


      StdVector< StdVector<Double> > Htotal_prev_;
      StdVector< StdVector<Double> > Mprev_iter_;

      StdVector< StdVector<Double> > HxS_prev_;
      StdVector< StdVector<Double> > HyS_prev_;
      StdVector< StdVector<Double> > HzS_prev_;
      StdVector< StdVector<Double> > MxS_prev_;
      StdVector< StdVector<Double> > MyS_prev_;
      StdVector< StdVector<Double> > MzS_prev_;

      StdVector< StdVector<Double> > HxS_n_;
      StdVector< StdVector<Double> > HyS_n_;
      StdVector< StdVector<Double> > HzS_n_;
      StdVector< StdVector<Double> > MxS_n_;
      StdVector< StdVector<Double> > MyS_n_;
      StdVector< StdVector<Double> > MzS_n_;

      StdVector< StdVector<Double> > HxS_n_tmp_;
      StdVector< StdVector<Double> > HyS_n_tmp_;
      StdVector< StdVector<Double> > HzS_n_tmp_;
      StdVector< StdVector<Double> > MxS_n_tmp_;
      StdVector< StdVector<Double> > MyS_n_tmp_;
      StdVector< StdVector<Double> > MzS_n_tmp_;

      // Mu-Tensor (only updated once per iteration). Reason for that is because there can be
      // several integration points inside one element and we only need to evaluate the hysteresis
      // operator once for each element.
      StdVector< Matrix<Double> > mu_;
      StdVector<bool> alreadyHasMu_;
      
      //! Pointer to math parser instance
      MathParser* mp_;

      UInt iterTracker4Mu_;
      double isMH_;

      std::string varHandle_;

      bool saveTmpStageVecs_;

    };
} //end of namespace
