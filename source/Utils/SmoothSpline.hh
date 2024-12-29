#ifndef FILE_SMOOTHSPLINE
#define FILE_SMOOTHSPLINE

#include<string>
#include "MatVec/Vector.hh"
#include "ApproxData.hh"

namespace CoupledField {

  //! computes approximation of sampled data with smoothed splines 
  //! (see Reitzinger&Kaltenabcher, SFB-Report 02-30)
  class SmoothSpline : public ApproxData
  {

  public:
    //! constructor getting x, y(x)
    SmoothSpline( std::string nlFncName, MaterialType matType );

    //! Copy Constructor
    SmoothSpline(const SmoothSpline& right)
      : ApproxData(right){
      this->size_ = right.size_;
      this->node_ = right.node_;
      this->ind_ = right.ind_;
      this->mu_ = right.mu_;
      this->mat_ = right.mat_;
      this->coef_ = right.coef_;
      this->rhs_ = right.rhs_;
      this->h_ = right.h_;
      this->g_ = right.g_;
      this->xStart_ = right.xStart_;
      this->xEnd_ = right.xEnd_;
      this->yEnd_ = right.yEnd_;
      this->theta_ = right.theta_;
      this->delta_ = right.delta_;
      this->nuMax_ = right.nuMax_;
      this->yMax_ = right.yMax_;
      this->extrapolAlpha_ = right.extrapolAlpha_;
      this->extrapolBeta_ = right.extrapolBeta_;
    }

    virtual SmoothSpline* Clone(){
      return new SmoothSpline( *this );
    }

    //! destructor
    virtual ~SmoothSpline();

    //computes the approximation polynom
    void CalcApproximation( const bool start=true ) override;

    //! computes the regularization parameter
    void CalcBestParameter();

    //! set accuracy of measured data
    void SetAccuracy( double val ) {
      delta_ = val;
    };

    //! set maximal y-value
    void SetMaxY( double val ) {
      yMax_ = val;
    };

    //! returns f(x)
    double EvaluateFunc( double x ) const override;

    //! returns  f'(x)
    double EvaluateDeriv( double x ) const override;

    //! returns grad f(x)
    Vector<double> EvaluatePrime(const Vector<double>& p) const override {
      Vector<double> ret(1,0.0);
      ret[0] = EvaluateDeriv(p[0]);
      return ret;
    }

    // returns inverse of y(x) 
    double EvaluateFuncInv( double t ) const;

    //! returns derivative of inverse of y(x)
    double EvaluatePrimeInv( double t ) const;

    //! computes the magnetic reluctivity
    double EvaluateFuncNu(double t) const override {
      // if (t > yMax_ ) {
      //   t = yMax_;
      // }
      if ( t <= yEnd_ ) {
        return ( EvaluateFuncInv(t)/t );
      }
      else {
        return ( nuMax_ + extrapolBeta_*std::exp(-extrapolAlpha_*t) );
      };
    }

    //! computes the derivative of magnetic reluctivity
    double EvaluatePrimeNu(double t) const override {
      // if (t > yMax_ ) 
      //   t = yMax_;
      if ( t <= yEnd_ ) {
        return ( (EvaluatePrimeInv(t)*t - EvaluateFuncInv(t))/(t*t) ); 
      }
      else {
        return ( -extrapolAlpha_*extrapolBeta_*std::exp(-extrapolAlpha_*t) );
      }
    }

    //! returns number of sampled data
    Integer GetSize() const { return numMeas_; };

    //! get original sampled y value
    double EvaluateOrigB( Integer i ) const { return y_[i]; };

    //! evalutes original sampled reluctivity
    double EvaluateOrigNu( Integer i ) const { return x_[i]/y_[i]; };

    //! prints out original and approximated function
    void Print() const override;

  private:


    //! sets up the system matrix
    void ConstructMatrix();

    //! computes right hand side
    void ConstructRHS( const Vector<double>& y );

    //! computes the coefficients of approximation polynom
    void CalcCoef();

    //! evalutes the Hermit functions
    double HermiteFunc( double t, Integer i ) const;

    //! evalutes the derivative of Hermit functions
    double HermitePrime( double t, Integer i ) const;

    //! gets the interval, defined by two sampled points
    Integer GetInterval( double t ) const;

    //! newton iteration for evaluation of iinverse function
    double Newton( double f, double start=1 ) const;

    //! computes the starting values for Newton-iterations; stores it in g_
    void CalcStart();

    //! checks, if BH curve is monoton 
    bool MonotoneBH();

    //! checks iv nu-B curve is monotone
    bool MonotoneNu();

    Integer size_;  //!< size of system matrix
    Integer node_;  //!< number of internal sampled points:= numMeas_ - 2

    Integer ind_;   //!< number of evaluation points for approximated curve

    double mu_;             //!< discrepancy parameter
    Vector<double> mat_;    //!< system matrix
    Vector<double> coef_;   //!< array, containing the coefficients of the spline functions
    Vector<double> rhs_;    //!< right hand side of algebraic system
    Vector<double> h_;      //!< intervals: x[i+1] - x[i]
    Vector<double> g_;      //!< contains the starting values for the Newton-iteration


    double xStart_ = -1.0;    //!< first measured x-value; should be zero
    double xEnd_ = -1.0;      //!< last measured x-value
    double yEnd_ = -1.0;      //!< kast measured y-value
    double theta_ = -1.0;     //!< increment for y; (y[0] - y[end]) / ind_
    double delta_ = -1.0;     //!< measurement error

    double nuMax_ = -1.0;     //!< maximal value of reluctivity
    double yMax_ = -1.0;      //!< maximal value for y

    double extrapolAlpha_ = -1.0;  //!< paramter for extrapolation of reluctivity
    double extrapolBeta_  = -1.0;   //!< paramter for extrapolation of reluctivity

    //============================ just for testing

    //! own method, for reading in sampled data from file
    void Read();

  };

} //end of namespace


#endif
