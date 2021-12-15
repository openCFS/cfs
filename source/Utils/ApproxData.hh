#ifndef FILE_ApproxData
#define FILE_ApproxData

#include <string>
#include "General/Environment.hh"
#include "MatVec/Vector.hh"

namespace CoupledField {

  //! Base class for approximation of sampled data
  class ApproxData
  {
  public:

    ApproxData() {};

    //! constructor
    ApproxData( std::string nlFncName, MaterialType matType , UInt numIndep = 1);

    //! destructor: nothing to do
    virtual ~ApproxData() {};

    //! Copy Constructor
    ApproxData(const ApproxData& right) {
      //here we would also need to create a new operator
      this->factor_ = right.factor_;
      this->matType_ = right.matType_;
      this->nlFileName_ = right.nlFileName_;
      this->numIndepend_ = right.numIndepend_;
      this->numMeas_ = right.numMeas_;
      this->x1_ = right.x1_;
      this->x_ = right.x_;
      this->y_ = right.y_;
      this->z_ = right.z_;
    }

    //! reads in the sampled data form file with name fncName
    void ReadNlinFunc( std::string fncName );
    void ReadNlinFuncTwoIndep( std::string fncName );
    void ReadNlinFuncThreeIndep( std::string fncName );

    //! return file name
    inline std::string GetNlFileName() {
      return nlFileName_;
    }
    
    //! perform checks
    void PerformChecksOnInputData( std::string fncName );

    //! performs the approximation
    virtual void CalcApproximation( const bool start=true ) = 0;

    //! find lower \param klo, and higher \param khi index of the boundary values of \param x in a monotonically 
    //! increasing double dataset \param axis. 
    //! Also returns difference.
    void findBracketIndices(const double &x, const Vector<double> & axis, UInt & klo, UInt & khi, double &diff) const;

    //! evaluates the functions
    virtual double EvaluateFunc(const double x) const {
      EXCEPTION("only implemented in LinInterpolate");
      return 0;
    }

    virtual double EvaluateFunc(const double x, const double y) const {
      EXCEPTION("only implemented in BiLinInterpolate and BiCubicInterpolate");
      return 0;
    }

    virtual double EvaluateFunc(const double x, const double y, const double z) const {
      EXCEPTION("only implemented in TriLinInterpolate and TriCubicInterpolate");
      return 0;
    }

    virtual double EvaluateFunc(const Vector<double>& p) const {
      EXCEPTION("only implemented in derived classes");
      return 0;
    }

    //! returns d f(x) / dx
    virtual double EvaluateDeriv(const double x) const {
      EXCEPTION("only implemented in LinInterpolate and CubicInterpolate");
      return 0;
    }

    virtual double EvaluateDeriv(const Vector<double>& p, const int dparam) const {
      EXCEPTION("only implemented in derived classes");
      return 0;
    }

    //! evaluates the gradient of the function
    virtual Vector<double> EvaluatePrime(const Vector<double>& p) const {
      EXCEPTION("only implemented in derived classes");
      Vector<double> ret(p.GetSize(),0);
      return ret;
    }

    //! computes the magnetic reluctivity
    virtual double EvaluateFuncNu(double t) const {
      EXCEPTION(" ApproxData: EvaluateFuncNu not implemented");
      return -1.0;
    }

    //! computes the derivative of magnetic reluctivity
    virtual double EvaluatePrimeNu(double t) const {
      EXCEPTION(" ApproxData: EvaluatePrimeNu not implemented");
      return -1.0;
    }

    //! return the material type
    virtual MaterialType GetMatType() {
      return matType_;
    }

    //! prints out original and approximated function
    virtual void Print() const {};

    inline void SetFactor(double f) { factor_=f;}

  protected:

    // we use type initialization by {} to suppress the constructor warnings about uninitialized members
    Vector<double> x_{};  //!< independent value
    Vector<double> x1_{};  //!< independent value
    Vector<double> y_{};  //!< function value
    StdVector< Vector <double> > z_{}; //! function value for two independent vars
    UInt numMeas_{};      //!< number of sampled points
    UInt numIndepend_{};  //!< number of independent values (default 1)
    MaterialType matType_{}; //!< material parameter to be approximated
    std::string nlFileName_{}; //!< name of file
    StdVector< std::string > slicesFiles_{}; //! holding files containing 2d slices (needed for 3d data)

    double factor_{}; //! A factor (default 1.) which can be used to scale the value of the approximated function

  };


} //end of namespace


#endif

