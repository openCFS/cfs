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
      if(lin_nu_r_[elem->regionId] < 0);
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
  virtual void PostInit();

  /** Get's linear or nonlinear nu_r from the elements material property */
  double CalcRelactivity(const Elem* elem, UInt dim);



  /** overloads SIMP::CalcFunction()
   * @see ErsatzMaterial::CalcFunction */
  double CalcFunction(Excitation& excite, Function* f, bool derivative);

  /** [1 0; 0 0] or [0 0; 0 1] for SQR_MAG_FLUX_DENS_X/Y */
  const Matrix<double>& GetSelectionMatrix(const Function* f) const;



  /** @see ErsatzMaterial::FillRealAdjointRHS() */
  virtual bool FillRealAdjointRHS(Excitation& excite, Function* f, Vector<double>& rhs)
  {
    if(f->GetType() == Function::SQR_MAG_FLUX_DENS_X || f->GetType() == Function::SQR_MAG_FLUX_DENS_Y)
    {
      CalcMagFluxAdjRHS(excite, f, rhs);
      return true;
    }
    return false;
  }

  
  virtual void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0)
  {
    if(f->ctxt->IsComplex())
      SetElementK<std::complex<double>, double >( f, de, tf, app, out, derivative, calcMode, ev);
    else
      SetElementK<double, double>( f, de, tf, app, out, derivative, calcMode, ev);
  }

  virtual void SetElementRHS(DesignElement* de, const TransferFunction* tf, App::Type app, SingleVector* out, CalcMode calcMode, bool derivative = true)
  {
    /*if(f->ctxt->IsComplex())
      SetElementRHS<std::complex<double>, double >( de, tf, app, out, calcMode, derivative);
    else*/
      SetElementRHS<double, double>( de, tf, app, out, calcMode, derivative);
  }

private:

  template <class T1, class T2>
  void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0);

  template <class T1, class T2>
    void SetElementRHS(DesignElement* de, const TransferFunction* tf, App::Type app, SingleVector* out, CalcMode calcMode, bool derivative = true);

  /** magnetic flux density */
  void CalcMagFluxAdjRHS(Excitation& excite, Function* f, Vector<double>& out);

  /** calc magnetic flux density as 1/N * sum<J,B*A> where B is the BOp from the BDBInt (here the curl operator) and A is the
   * element solution vector (scalar) and J is either [1,0] or [0,1] to select the horizontal or vertical part and N averages over
   * the sum of the N elements of f->region. The scalar product is evaluates over the integration points */
  double CalcMagFluxDensity(Excitation& excite, Function* f);

  /** Calculate the magnetic flux density gradient. The weight is always 1 as the magnetic flux density needs to be per excitation */
  void CalcMagFluxDensGradient(Excitation& excite, Function* f,  TransferFunction* tf);

  /** array by index region-id which says true for nonlinear regions (here CurlCulrIntegrator-NL) */
  StdVector<bool> nonlin_;

  /** cache for region constant nu_r in the linear case. Indexed by region-id. -1.0 is uninitialized */
  StdVector<double> lin_nu_r_;

  Matrix<double> sel_x_;
  Matrix<double> sel_y_;

  /** The magnetic rhs */
  DesignDependentRHS magRHS;
};

} // end of namespace


#endif /*MAGSIMP_HH_*/
