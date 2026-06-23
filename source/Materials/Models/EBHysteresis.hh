

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

      //simple diagonal approximation: mu_ii = |dB_i| / |dH_i|, off-diagonals set to zero
      Matrix<Double> EvaluateLocalMu(StdVector<Double> E, StdVector<Double> D, UInt idx);

      //DFP (Davidon-Fletcher-Powell) quasi-Newton update of dB/dH (symmetric rank-2 update).
      Matrix<Double> EvaluateLocalMuDFP(StdVector<Double> E, StdVector<Double> D, UInt idx);

      //Broyden (GBM) quasi-Newton update of dB/dH (rank-1 update)
      Matrix<Double> EvaluateLocalMuGBM(StdVector<Double> E, StdVector<Double> D, UInt idx);

      //BFGS quasi-Newton update of dB/dH (symmetric rank-2 update)
      Matrix<Double> EvaluateLocalMuBFGS(StdVector<Double> E, StdVector<Double> D, UInt idx);

      //forward finite difference approximation of dB/dH. Used as initial Jacobian in the
      //first Newton iteration when no quasi-Newton history is available yet
      Matrix<Double> EvaluateLocalMuFiniteDifferences(Vector<Double> E, StdVector<Double> D, UInt idx);

      //Computes dB/dH analytically for the purely anhysteretic case using the atan anhysteretic function
      Matrix<Double> EvaluateLocalMuAnhystersisOnly(Vector<Double> E, UInt idx);

      //Computes the inverse of a 3x3 matrix using the adjugate method: A^{-1} = adj(A) / det(A)
      StdVector<Double> inv3x3(StdVector<Double> A);

      /*
      =========================================================================================
          Methods regarding the evaluation of the hysteresis model
      =========================================================================================
      */

      //used the correct Eval_*() method based on approx_type_ and anhyst_type_
      Vector<Double> Evaluate(Vector<Double> E, UInt idx);

      //2D Energy-Based Model (EBM) with atan anhysteretic function.
      //Solves the dual variational problem (Eq. 2.35) via 1D Newton iteration over the angle phi
      Vector<Double> Eval_2D_EBM_ATAN(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      //2D Vector Play Model (VPM) with analytic anhysteretic function (atan or Pacejka)
      Vector<Double> Eval_2D_VPM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      //3D Vector Play Model (VPM) with analytic anhysteretic function (atan or Pacejka)
      Vector<Double> Eval_3D_VPM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      //2D Vector Play Model (VPM) with SMSM multiscale anhysteretic model
      Vector<Double> Eval_2D_VPM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      //3D Vector Play Model (VPM) with SMSM multiscale anhysteretic model
      Vector<Double> Eval_3D_VPM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      //2D Energy-Based Model (EBM) with SMSM multiscale anhysteretic model
      Vector<Double> Eval_2D_EBM_MSM(Vector<Double> Hn, bool saveTmpStageVecs, UInt idx, StdVector<Double> weight, StdVector<Double> chi);

      //3D Energy-Based Model (EBM) with SMSM multiscale anhysteretic model
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
      Double m_sat_pacejka_;
      Double a_pacejka_;
      Double b_pacejka_;
      Double c_pacejka_;
      Double mu0_;
      UInt numS_;
      Double chi_factor_;
      UInt jacobian_method_;
      //UInt anhyst_type_;
      std::string anhyst_type_;
      std::string anhyst_formula_;
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

      //! Cached pointer to iteration counter for fast access
      Double* iterationCounterPtr_;

      UInt iterTracker4Mu_;
      double isMH_;

      std::string varHandle_;

      bool saveTmpStageVecs_;

    };
} //end of namespace
