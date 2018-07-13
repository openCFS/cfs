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
  void GetRelactivity(Elem* elem);

  /** material constant for convenience */
  double nu_0;

protected:

  /** @see Optimization::PostInit() */
  virtual void PostInit();

  /** overloads SIMP::CalcFunction()
   * @see ErsatzMaterial::CalcFunction */
  double CalcFunction(Excitation& excite, Function* f, bool derivative);

  /** [1 0; 0 0] or [0 0; 0 1] for SQR_MAG_FLUX_DENS_X/Y */
  const Matrix<double>& GetSelectionMatrix(const Function* f) const;

  /** array by index regiond-id which says true for nonlinear regions (here CurlCulrIntegrator-NL) */
  StdVector<bool> nonlin;

protected:
  
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

private:

  template <class T1, class T2>
  void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0);

  /** magnetic flux density */
  void CalcMagFluxAdjRHS(Excitation& excite, Function* f, Vector<double>& out);

  /** calc magnetic flux density as 1/N * sum<J,B*A> where B is the BOp from the BDBInt (here the curl operator) and A is the
   * element solution vector (scalar) and J is either [1,0] or [0,1] to select the horizontal or vertical part and N averages over
   * the sum of the N elements of f->region. The scalar product is evaluates over the integration points */
  double CalcMagFluxDensity(Excitation& excite, Function* f);

  /** Calculate the magnetic flux density gradient. The weight is always 1 as the magnetic flux density needs to be per excitation */
  void CalcMagFluxDensGradient(Excitation& excite, Function* f,  TransferFunction* tf);

  Matrix<double> sel_x_;
  Matrix<double> sel_y_;

};

} // end of namespace


#endif /*MAGSIMP_HH_*/
