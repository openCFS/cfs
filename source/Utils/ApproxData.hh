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

    //! constructor
    ApproxData( std::string nlFncName, MaterialType matType , UInt numIndep = 1);

    //! destructor: nothing to do
    virtual ~ApproxData() {;};

    //! Copy Constructor
    ApproxData(const ApproxData& right){
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
    std::string GetNlFileName() {
      return nlFileName_;
    }
    
    //! perform checks
    void PerformChecksOnInputData( std::string fncName );

    //! performs the approximation
    virtual void CalcApproximation( bool start=true ) = 0;

    //! set accuracy of measured data
    virtual void SetAccuracy( Double val ) {
      EXCEPTION(" ApproxData: SetAccuracy not implemented");
    };

    //! set maximal y-value
    virtual void SetMaxY( Double val ) {
      EXCEPTION(" ApproxData: SetMaxY not implemented");
    };


    //! find lower \param klo, and higher \param khi index of the boundary values of \param x in a monotonically 
    //! increasing double dataset \param axis. 
    //! Also returns difference.
    void findBracketIndices(const Double &x, const Vector<Double> & axis, UInt & klo, UInt & khi, Double &diff);

    //! evaluates the functions
    virtual Double EvaluateFunc( Double t ) = 0;
    virtual Double EvaluateFunc(Double x, Double y) {EXCEPTION("only implemented in BiLinInterpolate"); Double a=-1.; return a;} // bilinear
    virtual Double EvaluateFunc(Double x, Double y, Double z) {EXCEPTION("only implemented in TriLinInterpolate"); Double a=-1.; return a;} //trilinear

    //! evalutates the derivazive of the function
    virtual Double EvaluatePrime( Double t ) = 0;

    //! evalutates the invers of the function
    virtual Double EvaluateFuncInv( Double t ) = 0;

    //! valuates the inverse of the derivative of the function
    virtual Double EvaluatePrimeInv( Double t ) = 0;

    //! computes the best approximation
    virtual void CalcBestParameter() = 0;

    //! computes the magnetic reluctivity
    virtual Double EvaluateFuncNu(Double t) {     
      EXCEPTION(" ApproxData: EvaluateFuncNu not implemented");
      return -1.0; 
    }

    //! computes the derivative of magnetic reluctivity
    virtual Double EvaluatePrimeNu(Double t) {
      EXCEPTION(" ApproxData: EvaluatePrimeNu not implemented");
      return -1.0; 
    }

    //! return the material type
    virtual MaterialType GetMatType() {
      return matType_;
    }

    //! prints out original and approximated function
    virtual void Print() {;};

    void SetFactor(Double f) { factor_=f;}

  protected:

    Vector<Double> x_;  //!< independent value
    Vector<Double> x1_;  //!< independent value
    Vector<Double> y_;  //!< function value
    StdVector< Vector <Double> > z_; //! function value for two independent vars
    UInt numMeas_;      //!< number of sampled points
    UInt numIndepend_;  //!< number of independent values (default 1)
    MaterialType matType_; //!< material parameter to be approximated
    std::string nlFileName_; //!< name of file
    StdVector< std::string > slicesFiles_; //! holding files containing 2d slices (needed for 3d data)

    Double factor_; //! A factor (default 1.) which can be used to scale the value of the approxmated function

  };


} //end of namespace


#endif

