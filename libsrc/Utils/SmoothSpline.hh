#ifndef FILE_SMOOTHSPLINE
#define FILE_SMOOTHSPLINE

#include<string>
#include "ApproxData.hh"

namespace CoupledField {

  //! \todo Replace C Arrays in constructor
  //! by normal vectors and improve error robustness,
  //! if position gets out of range
class SmoothSpline : public ApproxData
{
public:
  //! constructor getting x, y(x)
  SmoothSpline(std::string nlFncName);

  //! destructor
  virtual ~SmoothSpline();

  //computes the approximation polynom
  virtual void CalcApproximation(int start=1);

  //! computes the regularization parameter
  virtual void CalcBestParameter();

  ///
  virtual void CalcMonotoneParameter();


  //! returns y(x)
  virtual double EvaluateFunc(double x);

  //! returns  y'(x)  
  virtual double EvaluatePrime(double x);

  ///
  virtual double EvaluateFuncInv(double t);

  ///
  virtual double EvaluatePrimeInv(double t);

  ///
  int GetSize() {return nummeas;};

  ///
  double EvaluateOrigB(int i) {return y[i];};

  ///
  double EvaluateOrigNu(int i) {return x[i]/y[i];};

  ///
  void Print(double *x, double *y, double lower, double upper);

private:


  ///
  void ConstructMatrix();

  ///
  void ConstructRHS(double * y);

  ///
  void CalcCoef();
  ///
  void EvaluateInv(double v, double *f, double *p);

  ///
  double HermiteFunc(double t, int i);

  ///
  double HermitePrime(double t, int i);

  ///
  int GetInterval(double t);

  ///
  double Newton(double f,double start=1);

  ///
  void CalcStart();

  ///
  int MonotoneBH();

  ///
  int MonotoneNu();

  int size;
  int node;

  int ind;

  double mu;
  double * mat;
  double * coef;
  double * rhs;
  double * h;
  double * g;


  double alpha;
  double beta;
  double theta;
  double delta;

  //============================ just for testing
  ///
  void Read();

  ///
  void MakeOutput(double * x, double * y, double left, double right);

  ///
  void MakeOutputInv(double * x, double * y, double left, double right);


};

} //end of namespace


#endif
