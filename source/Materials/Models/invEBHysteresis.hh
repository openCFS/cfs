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

      void Init(std::map<std::string, double> ParameterMap, std::map<std::string, string> StringParameterMap, shared_ptr<ElemList> entityList, UInt dim);

      // This shall act as a unified handling for updating states from one timestep to another, triggered by the SolveStepEB.
      // Originally this was handled via flags but manually updating the states (via an explicit call in SolveStepEB) is way 
      // safer and probably also faster since we do not need to check the states for every element at every iteration.
      void UpdateStates();
      void AllowUpdates(bool allow);
      
      // --------------- General hysteresis operator functions (start) --------------- //

      // recipe to evaluate the hysteresis operator
      Matrix<Double> ComputeTensorialMaterialParameter(Vector<Double> B, Integer ElemNum);

      // Evaluate the hystersis operator for all available cases
      Vector<Double> Evaluate(Vector<Double> B, UInt idx);
      
      // Available hysteresis / nonlinear models
      Vector<Double> Eval_2D_invEBM(Vector<Double> Bn, bool saveTmpStageVecs, UInt idx); // 2D exact inversed hyyteresis operator (LAVET)
      Vector<Double> Eval_2D_VSM(Vector<Double> B_n, bool saveTmpStageVecs, UInt idx); // 2D Vector Stop Model (VSM)
      Vector<Double> Eval_3D_invEBM(Vector<Double> Bn, bool saveTmpStageVecs, UInt idx); // 3D exact inversed hyyteresis operator (LAVET)
      Vector<Double> Eval_3D_VSM(Vector<Double> B_n, bool saveTmpStageVecs, UInt idx); // 3D Vector Stop Model (VSM)

      // --------------- General hysteresis operator functions (end) --------------- //

      // --------------- FE related functions --------------- //

      // just for the computation of the residual, we do not store anything here
      Vector<Double> GetFieldIntensity(Vector<Double> B, Integer ElemNum);

      // this method is for the FE - Formulation and give the material tensor: Broyden update
      // Matrix<Double> EvaluateLocalNuBroyden(StdVector<Double> H, StdVector<Double> B, UInt idx);

      // this methods is for the FE - Formulation and give the material tensor: BFGS update
      Matrix<Double> EvaluateLocalNuBFGS(Vector<Double> dH, Vector<Double> dB, UInt idx);

      // this methods is for the FE - Formulation and give the material tensor: Broyden update
      //Matrix<Double> EvaluateLocalNuGBM(Vector<Double> dH, Vector<Double> dB, UInt idx);

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

      // --------------- HELPER FUNCTIONS: evaluation of the vector stop model (START) --------------- //

      Vector<Double> Anhyst_VSM(Vector<Double> B_n);
      Double Root_Function_VSM(Double norm_J, Double norm_B_n);


      // --------------- HELPER FUNCTIONS: evaluation of the vector stop model (END) --------------- //

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

      // --------------- HELPER FUNCTIONS: Handle lookup tables for anhysteresis (START) --------------- //
      
      Double LinInterp1D(Double norm_J, std::vector<Double> J_lut, std::vector<Double> H_lut);
      Double FiniteDifference1D(Double norm_J, std::vector<Double> J_lut, std::vector<Double> H_lut);
      Double IntTrapz1D(Double norm_J, std::vector<Double> J_lut, std::vector<Double> H_lut);



      // --------------- HELPER FUNCTIONS: Handle lookup tables for anhysteresis (START) --------------- //


    private:
      //==============

      // spatial dimension of the problem (2 or 3)
      UInt dim_; 

      std::map<UInt,UInt> ElemNum2Idx_;

      // number of elements in the mesh
      UInt numElems_;
      // sautration polarization [T]
      Double Js_;
      // saturation magnetization [A/m] for atan/Pacejka anhysteretic function
      Double Ms_;
      // path to the lookup table file for the anhysteretic function
      string lookup_table_file_;
      // parameters for Pacejka anhysteretic function
      Double A_, pa_, pb_, pc_;
      // parameters for Brauer anhysteretic function
      Double p_0_, p_1_, p_2_;
      // vauum permeability [Vs/Am]
      Double mu0_;
      // number of pseudo particles in the composite model
      UInt numS_;
      // method for quasi-Newton approximation of the Jacobian: 4 = BFGS, other methods (finite differences, DFP, GBM)
      // not implemented yet
      UInt jacobian_method_;
      // type of anhysteretic model: currently only "analytic_anhysteresis" implemented
      std::string anhyst_type_;
      // analytic anhysteretic formula: tan, atan, pacejka, brauer or lookup table
      std::string anhyst_formula_;
      // inverse hysteresis approximation type: "fullInvEB" (Newton optimization) or "approxVSM" (Vector Stop Model)
      std::string approx_type_;
      // ath to the file containing pinning forces (kappa) and weights (omega) for the composite model
      std::string pinning_forces_weight_;
      UInt ndx_g_;
      // lookup table for the anhysteretic function: J values [T] and corresponding H values [A/m]
      std::vector<Double> J_lut_, H_lut_;
      StdVector<Double> kappa_file_, omega_file_;

      // Flux density B and field H from the previous Newton iteration (for quasi-Newton update)
      StdVector< StdVector<Double> > B_prev_;
      StdVector< StdVector<Double> > H_prev_;

      // B from the previous time step for the VSM and corresponding temporaty buffers
      StdVector< StdVector<Double> > B_prev_VSM_;
      StdVector< StdVector<Double> > H_prev_VSM_;      
      StdVector< StdVector<Double> > B_prev_VSM_tmp_;
      StdVector< StdVector<Double> > H_prev_VSM_tmp_;

      // internal state vectors of the VSM from previous time step and temporary buffers
      StdVector< StdVector<Double> > w_state_x_old_;
      StdVector< StdVector<Double> > w_state_y_old_;
      StdVector< StdVector<Double> > w_state_z_old_;
      StdVector< StdVector<Double> > w_state_x_tmp_;
      StdVector< StdVector<Double> > w_state_y_tmp_;
      StdVector< StdVector<Double> > w_state_z_tmp_;

      // state vectors 
      StdVector< StdVector<Double> > J_x_k_prev_;
      StdVector< StdVector<Double> > J_y_k_prev_;
      StdVector< StdVector<Double> > J_z_k_prev_;

      // magnetic polarization J_k [T] per pseudo particle and temporary buffers
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

      //! Cached pointer to iteration counter for fast access
      Double* iterationCounterPtr_;

      UInt iterTracker4Nu_;
      UInt timeStep_;
      double isMH_;

      std::string varHandle_;

      bool saveTmpStageVecs_;

    };
} //end of namespace
