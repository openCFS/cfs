
// ===================================================================================================
/*!
 *       \file     invEBHysteresis.hh
 *       \brief    This material model represents an inverse hysteresis operator H(B).
 *                 It uses an energy-based approach and minimizes the magnetic energy.
 *                 Therefore the following minimization problem has to be solved
 *                 
 *                        J = min{ 1/(2*mu_0)|B-J|^2 + u(J) + chi*|J-J_p| }.
 * 
 *                 Here u(J) is the internal energy, chi is a pinning force, B is the mag. flux density
 *                 and J is the mag. polarization.
 * 
 *                 Reference: On forward and inverse energy-based magnetic vector hysteresis operators
 *                            doi: 10.1109/TMAG.2025.3544507
 *
 *       \date     March, 2025
 *       \author   Lukas Domenig, TU Graz
 */
//====================================================================================================

#pragma once

#include <list>

#include <boost/utility.hpp>

#include "MatVec/Vector.hh"
#include "Domain/Domain.hh"
#include "Model.hh"


namespace CoupledField {

  class SolveStepEB;
  class MathParser;

  class invEBHysteresis : public Model {

    public:
      //! Constructor
      invEBHysteresis();

      //! Destructor
      virtual ~invEBHysteresis();

      void Init(std::map<std::string, double> ParameterMap, shared_ptr<ElemList> entityList, UInt dim);

      // This shall act as a unified handling for updating states from one timestep to another, triggered by the SolveStepEB.
      // Originally this was handled via flags but manually updating the states (via an explicit call in SolveStepEB) is way 
      // safer and probably also faster since we do not need to check the states for every element at every iteration.
      void UpdateStates();
      void AllowUpdates(bool allow);
      
      // --------------- General hysteresis operator functions --------------- //

      // recipe to evaluate the hysteresis operator
      Matrix<Double> ComputeTensorialMaterialParameter(Vector<Double> B, Integer ElemNum);

      // Evaluate the hystersis operator for all available cases
      Vector<Double> Evaluate(Vector<Double> B, UInt idx);
      
      Vector<Double> Eval_2D_invEBM_TAN(Vector<Double> Bn, bool saveTmpStageVecs, UInt idx);
      Vector<Double> Eval_3D_Brauer(Vector<Double> Bn, bool saveTmpStageVecs, UInt idx);
      Vector<Double> Eval_3D_invEBM_TAN(Vector<Double> Bn, bool saveTmpStageVecs, UInt idx);

      // --------------- FE related functions --------------- //

      // just for the computation of the residual, we do not store anything here
      Vector<Double> GetFieldIntensity(Vector<Double> B, Integer ElemNum);

      // this method is for the FE - Formulation and give the material tensor: Broyden update
      // Matrix<Double> EvaluateLocalNuBroyden(StdVector<Double> H, StdVector<Double> B, UInt idx);

      // this methods is for the FE - Formulation and give the material tensor: BFGS update
      Matrix<Double> EvaluateLocalNuBFGS(Vector<Double> dH, Vector<Double> dB, UInt idx);

      // this methods is for the FE - Formulation and give the material tensor: Broyden update
      Matrix<Double> EvaluateLocalNuGBM(Vector<Double> dH, Vector<Double> dB, UInt idx);

      // --------------- HELPER FUNCTIONS: evaluation of the hysteresis operator (START) --------------- //

      // calculates the energy of the minimization problem, i.e., the magnetic energy grad w_m
      // only needed for the line search procedure
      Double CalcMagEnergy(Vector<Double> B, Matrix<Double> J_k, Matrix<Double> J_k_prev);
      Double CalcMagEnergy3D(Vector<Double> B, Matrix<Double> J_k, Matrix<Double> J_k_prev);

      // calculates the gradient of the minimization problem, i.e., the magnetic energy grad w_m
      Vector<Double> CalcGradientMagEnergy(Vector<Double> B, Matrix<Double> J_k, Matrix<Double> J_k_prev);
      Vector<Double> CalcGradientMagEnergy3D(Vector<Double> B, Matrix<Double> J_k, Matrix<Double> J_k_prev);

      // calculates the hessian of the minimization problem, i.e., the magnetic energy grad^2 w_m
      Matrix<Double> CalcHessianMagEnergy(Matrix<Double> J_k, Matrix<Double> J_k_prev);
      Matrix<Double> CalcHessianMagEnergy3D(Matrix<Double> J_k, Matrix<Double> J_k_prev);

      // calculates the internal energy of the magnetic energy u(J)
      Vector<Double> CalcInternalEnergy(Matrix<Double> J_k);
      Vector<Double> CalcInternalEnergy3D(Matrix<Double> J_k);

      // calculates the gradient of the internal energy of the magnetic energy grad u(J)
      Vector<Double> CalcGradientInternalEnergy(Vector<Double> J);
      Vector<Double> CalcGradientInternalEnergy3D(Vector<Double> J);

      // calculates the hessian of the internal energy of the magnetic energy grad^2 u(J)

      Matrix<Double> CalcHessianInternalEnergy(Vector<Double> J);
      Matrix<Double> CalcHessianInternalEnergy3D(Vector<Double> J);

      // calculates the hessian of the dissipation energy of the magnetic energy grad^2 d(J)

      Matrix<Double> CalcHessianDissipationEnergy(Vector<Double> delta_J, Double norm_delta_J, Double chi);
      Matrix<Double> CalcHessianDissipationEnergy3D(Vector<Double> delta_J, Double norm_delta_J, Double chi);

      // applies an Armijo line search in the minimization of the magnetic energy
      Matrix<Double> LineSearchArmijoInvEB(Double phi_der0, Vector<Double> B, Matrix<Double> J_k, Matrix<Double> delta_J_k, Matrix<Double> J_k_prev);
      Matrix<Double> LineSearchArmijoInvEB3D(Double phi_der0, Vector<Double> B, Matrix<Double> J_k, Matrix<Double> delta_J_k, Matrix<Double> J_k_prev);

      // --------------- HELPER FUNCTIONS: evaluation of the hysteresis operator (END) --------------- //

      // --------------- HELPER FUNCTIONS: Solve linear systems (START) --------------- //

      // solves a linear system using gaussian elimination
      template<int DIM>
      Vector<Double> LUSolve(Matrix<Double> A, Vector<Double> b);

      template<int DIM>
      bool LUDecompose(Matrix<Double> A, Matrix<Double>& L, Matrix<Double>& U);

      template<int DIM>
      Vector<Double> LUForwardSubstitution(const Matrix<Double>& L, const Vector<Double>& b);

      template<int DIM>
      Vector<Double> LUBackwardSubstitution(const Matrix<Double>& U, const Vector<double>& y);

      // --------------- HELPER FUNCTIONS: Solve linear systems (END) --------------- //

    private:
      //==============
      UInt dim_; 

      std::map<UInt,UInt> ElemNum2Idx_;

      UInt numElems_;

      Double Js_;
      Double A_;
      Double mu0_;
      UInt numS_;
      Double chi_factor_;
      UInt jacobian_method_;
      UInt anhyst_type_;
      UInt ndx_g_;


      StdVector< StdVector<Double> > B_prev_;
      StdVector< StdVector<Double> > H_prev_;

      StdVector< StdVector<Double> > J_x_k_prev_;
      StdVector< StdVector<Double> > J_y_k_prev_;
      StdVector< StdVector<Double> > J_z_k_prev_;

      StdVector< StdVector<Double> > J_x_k_n_;
      StdVector< StdVector<Double> > J_y_k_n_;
      StdVector< StdVector<Double> > J_z_k_n_;

      StdVector< StdVector<Double> > J_x_k_n_tmp_;
      StdVector< StdVector<Double> > J_y_k_n_tmp_;
      StdVector< StdVector<Double> > J_z_k_n_tmp_;

      // Nu-Tensor (only updated once per iteration). Reason for that is because there can be
      // several integration points inside one element and we only need to evaluate the hysteresis
      // operator once for each element.
      StdVector< Matrix<Double> > nu_;
      StdVector<bool> alreadyHasNu_;

      //! Pointer to math parser instance
      MathParser* mp_;

      UInt iterTracker4Nu_;
      UInt timeStep_;
      double isMH_;

      std::string varHandle_;

      bool saveTmpStageVecs_;

    };
} //end of namespace