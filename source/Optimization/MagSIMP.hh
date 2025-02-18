#ifndef MAGSIMP_HH_
#define MAGSIMP_HH_

#include "Optimization/SIMP.hh"


namespace CoupledField
{
/** Extension from lin elast SIMP to the (non)-linear magnetic  case */
class MagSIMP : public SIMP
{
public:
  MagSIMP();

  virtual ~MagSIMP();


  /** the relactivity nu_r is part of the material property nu_0 * nu_r. In the linear case this is region constant,
   * in the nonlinear case this is element specific. */
  double GetRelactivity(const Elem* elem, UInt dim)
  {
    // linear or nonlinear case?
    if(nonlin_[elem->regionId])
      return CalcRelactivity(elem, dim);
    else
    {
      if(lin_nu_r_[elem->regionId] < 0)
        lin_nu_r_[elem->regionId] = CalcRelactivity(elem, dim);
      return lin_nu_r_[elem->regionId];
    }
  }

  /** material constant for convenience */
  static double nu_0;

  /** also needed in DesignSpace::ApplyPhysicalDesign*().
   * Extracts nu_r from the material coefficient. In the linear case this is a material constant, in the nonlinear case
   * this is element specific
   * @param elem if NULL we assume linear */
  static double ExtractRelactivity(CoefFunction* org_mat, const Elem* elem = NULL, UInt dim = 0);

protected:

  /** @see Optimization::PostInit() */
  void PostInit() override;

  /** Get's linear or nonlinear nu_r from the elements material property */
  double CalcRelactivity(const Elem* elem, UInt dim);



  /** overloads SIMP::CalcFunction()
   * @see ErsatzMaterial::CalcFunction */
  double CalcFunction(Excitation& excite, Function* f, bool derivative) override;

  /** [1 0; 0 0] or [0 0; 0 1] for SQR_MAG_FLUX_DENS_X/Y */
  const Matrix<double>& GetSelectionMatrix(const Function* f) const;

  /** @see ErsatzMaterial::FillRealAdjointRHS() */

  bool FillRealAdjointRHS(Excitation& excite, Function* f, Vector<double>& rhs) override;
  
  /** @see ErsatzMaterial::FillComplexAdjointRHS() */

  bool FillComplexAdjointRHS(Excitation& excite, Function* f, Vector<Complex>& rhs) override;

  /** See ErsatzMaterial::SetElementK() */
  void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0) override
  {
    if(f->ctxt->IsComplex())
      SetElementK<std::complex<double>, double >( f, de, tf, app, out, derivative, calcMode, ev);
    else
      SetElementK<double, double>( f, de, tf, app, out, derivative, calcMode, ev);
  }

  /** See ErsatzMaterial::SubstractCalcU1KU2RHS() */
  void SubstractCalcU1KU2RHS(Function* f, TransferFunction* tf, DesignElement* de, DesignDependentRHS* rhs, SingleVector* mat_vec) override;

private:

  template <class T1, class T2>
  void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0);

  /** calc magnetic flux density as 1/N * sum<J,B*A> where B is the BOp from the BDBInt (here the curl operator) and A is the
     * element solution vector (scalar) and J is either [1,0] or [0,1] to select the horizontal or vertical part and N averages over
     * the sum of the N elements of f->region. The scalar product is evaluated over the integration points */
  double CalcMagFluxDensity(Excitation& excite, Function* f);

  /** calc magnetic flux density as 1/N * sum<J,B*A> where B is the BOp from the BDBInt (here the curl operator) and A is the
     * element solution vector (scalar) and J is either [1,0] or [0,1] to select the horizontal or vertical part and N averages over
     * the sum of the N elements of f->region. The scalar product is evaluated over the integration points. The result is scaled by
     * \sum_{e \in opt_space} \rho_e * vol(e) */
  double CalcMagFluxDensityLosses(Excitation& excite, Function* f);

  /** Calculate the magnetic flux density gradient. The weight is always 1 as the magnetic flux density needs to be per excitation. */
  void CalcMagFluxDensGradient(Excitation& excite, Function* f);

  /** Calculate the magnetic flux density losses gradient. The weight is always 1 as the magnetic flux density needs to be per excitation. */
  void CalcMagFluxDensGradientLosses(Excitation& excite, Function* f);

  /** magnetic flux density */
  void CalcMagFluxAdjRHS(Excitation& excite, Function* f, Vector<double>& out);

  /** magnetic flux density losses */
  void CalcMagFluxLossesAdjRHS(Excitation& excite, Function* f, Vector<double>& out);

  /** enriched shape functions which give an integration with the state vector in 2D.
   * Allows computation of <N,A>. N_e = sum_ip w_i jacdet_i N_i, where N_i is the FE-shape function.
   * @param form encodes the region we apply the excitation (current in coil)
   * @param N has size of unknowns of state but has only contributions for the region */
  void CalcN(LinearFormContext* form, Vector<double>& N);

  /** calc coupling as M^4/(L1*L2)^2 */
  double CalcMagCouplingReal(Excitation& excite, Function* f);

  /** calc coupling as TODO */
  double CalcMagCouplingComplex(Excitation& excite, Function* f);

  void CalcCouplingAdjRealRHS(Excitation& excite, Function* f, Vector<double>& out);

  void CalcCouplingAdjComplexRHS(Excitation& excite, Function* f, Vector<Complex>& out);

  /** helper for CalcCouplingAdjComplexRHS():
   *   // out = factor_N1 * N1 * Re(AB) + factor_N2 * N2 * Imag(AB) */
  void HadamardHelper(Vector<Complex>& out, Complex factor_N1, const Vector<double>& N1,  Complex factor_N2, const Vector<Complex>& A);

  void HadamardHelper2(Vector<Complex>& out, Complex factor_N1, const Vector<double>& N1,  Complex factor_N2, const Vector<double>& N2, const Vector<Complex>& AB);

  /** res = <N, real(A).^2> or <N, imag(A).^2> */
  double InnerHelper(const Vector<double>& N, const Vector<Complex>& A);


  /** Calculate the coupling gradient */
  void CalcCouplingGradient(Excitation& excite, Function* f,  TransferFunction* tf);

  /** array by index region-id which says true for nonlinear regions (here CurlCulrIntegrator-NL) */
  StdVector<bool> nonlin_;

  /** cache for region constant nu_r in the linear case. Indexed by region-id. -1.0 is uninitialized */
  StdVector<double> lin_nu_r_;

  Matrix<double> sel_x_;
  Matrix<double> sel_y_;

  /** The magnetic rhs */
  DesignDependentRHS magRHS;

  Matrix<double> sel_xy_;

  /** volume/area (axi/plane) of whole optimization domain; to be calculated either in first objective value (CalcMagFluxDensity())
   *  or first adjoint problem (CalcMagFluxAdjRHS()) calculation, depending on optimizer. To indicate that the volume is not yet
   *  calculated we set it to -1. */
  double opt_vol_ = -1.0;
};

} // end of namespace


#endif /*MAGSIMP_HH_*/
