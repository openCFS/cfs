// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SMOOTHSPLINE
#define FILE_SMOOTHSPLINE

#include<string>
#include "MatVec/vector.hh"
#include "ApproxData.hh"

namespace CoupledField {

  //! computes approximation of sampled data with smoothed splines 
  //! (see Reitzinger&Kaltenabcher, SFB-Report 02-30)
  class SmoothSpline : public ApproxData
  {

  public:
    //! constructor getting x, y(x)
    SmoothSpline( std::string nlFncName, MaterialType matType );

    //! destructor
    virtual ~SmoothSpline();

     //computes the approximation polynom
    virtual void CalcApproximation( bool start=true );

    //! computes the regularization parameter
    virtual void CalcBestParameter();

    //! set accuracy of measured data
    virtual void SetAccuracy( Double val ) {
      delta_ = val;
    };

    //! set maximal y-value
    virtual void SetMaxY( Double val ) {
      yMax_ = val;
    };

    //! returns y(x)
    virtual Double EvaluateFunc( Double x );

    //! returns  y'(x)  
    virtual Double EvaluatePrime( Double x );

    // returns inverse of y(x) 
    virtual Double EvaluateFuncInv( Double t );

    //! returns derivative of inverse of y(x)
    virtual Double EvaluatePrimeInv( Double t );

    //! computes the magnetic reluctivity
    Double EvaluateFuncNu(Double t) {
      if (t > yMax_ ) {
        t = yMax_;
      }
      if ( t < yEnd_ ) {
        return ( EvaluateFuncInv(t)/t );
      }
      else {
        return ( nuMax_ + extrapolBeta_*std::exp(-extrapolAlpha_*t) );
      };
    }

    //! computes the derivative of magnetic reluctivity
    Double EvaluatePrimeNu(Double t) {
      if (t > yMax_ ) 
        t = yMax_;
      if ( t < yEnd_ ) {
        return ( (EvaluatePrimeInv(t)*t - EvaluateFuncInv(t))/(t*t) ); 
      }
      else {
        return ( -extrapolAlpha_*extrapolBeta_*std::exp(-extrapolAlpha_*t) );
      }
    }

    //! returns number of sampled data
    Integer GetSize() { return numMeas_; };

    //! get original sampled y value
    Double EvaluateOrigB( Integer i) { return y_[i]; };

    //! evalutes original sampled reluctivity
    Double EvaluateOrigNu( Integer i ) { return x_[i]/y_[i]; };

    //! prints out original and approximated function
    void Print();

  private:


    //! sets up the system matrix
    void ConstructMatrix();

    //! computes right hand side
    void ConstructRHS( Vector<Double>& y );

    //! computes the coefficients of approximation polynom
    void CalcCoef();

    //! computes the invers of the approximated function
    void EvaluateInv( Double v, Double& f, Double& p );

    //! evalutes the Hermit functions
    Double HermiteFunc( Double t, Integer i );

    //! evalutes the derivative of Hermit functions
    Double HermitePrime( Double t, Integer i );

    //! gets the interval, defined by two sampled points
    Integer GetInterval( Double t );

    //! newton iteration for evaluation of iinverse function
    Double Newton( Double f, Double start=1 );

    //! computes the starting values for Newton-iterations; stores it in g_
    void CalcStart();

    //! checks, if BH curve is monoton 
    bool MonotoneBH();

    //! checks iv nu-B curve is monotone
    bool MonotoneNu();

    Integer size_;  //!< size of system matrix
    Integer node_;  //!< number of internal sampled points:= numMeas_ - 2

    Integer ind_;   //!< number of evaluation points for approximated curve

    Double mu_;             //!< discrepancy parameter
    Vector<Double> mat_;    //!< system matrix
    Vector<Double> coef_;   //!< array, containing the coefficients of the spline functions
    Vector<Double> rhs_;    //!< right hand side of algebraic system
    Vector<Double> h_;      //!< intervals: x[i+1] - x[i]
    Vector<Double> g_;      //!< contains the starting values for the Newton-iteration


    Double xStart_;    //!< first measured x-value; should be zero
    Double xEnd_;      //!< last measured x-value
    Double yEnd_;      //!< kast measured y-value
    Double theta_;     //!< increment for y; (y[0] - y[end]) / ind_
    Double delta_;     //!< measurement error

    Double nuMax_;     //!< maximal value of reluctivity
    Double yMax_;      //!< maximal value for y

    Double extrapolAlpha_;  //!< paramter for extrapolation of reluctivity
    Double extrapolBeta_;   //!< paramter for extrapolation of reluctivity

    //============================ just for testing

    //! own method, for reading in sampled data from file
    void Read();

    //! prints out original and approximated functions 
    void MakeOutput(Vector<Double>& x, Vector<Double>& y);

    //! prints out the inverse of the approximated fucntion
    void MakeOutputInv(Vector<Double>& x, Vector<Double>& y);

    //! prints out the reluctivity-flux curve
    void MakeOutputNu();

  };

} //end of namespace


#endif
