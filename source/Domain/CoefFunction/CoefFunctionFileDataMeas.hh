#ifndef COEFFUNCTIONFILEDATAMEAS_HH_
#define COEFFUNCTIONFILEDATAMEAS_HH_

#include "CoefFunction.hh"
#include "Domain/Domain.hh"

namespace CoupledField
{

/** A measured signal over time from a file, the same vector for all points. Used for the Hx/Bx
 *  components of a <fileData> element, see SinglePDE::SetBCsAndLoads().
 *  The file format is the following:
 *  * ASCII and line based, one line per time step in the order of the time steps
 *  * lines starting with hashtag # are ignored
 *  * first value is the time, followed by dim double values
 *  * valid separators are space, tab, comma, semicolon and multiples of it and mixing
 *  The attribute 'missing' of <fileData> has no meaning here, it only applies to the node based
 *  CoefFunctionFileData. */
class CoefFunctionFileDataMeas : public CoefFunction
{
public:

  /** <fileData> element for the Hx/Bx components of a transient PDE */
  CoefFunctionFileDataMeas(Domain* ptDomain, const std::string& pdename, PtrParamNode pn, int dim);

  virtual ~CoefFunctionFileDataMeas() { }

  string GetName() const override { return "CoefFunctionFileDataMeas"; }

  void GetScalar(double& scal, const LocPointMapped& lpm) override {
     EXCEPTION("GetScalar not supported");
  }

  void GetScalar(Complex& scal, const LocPointMapped& lpm)  override{
    EXCEPTION("GetScalar not supported");
  }

  /** the same vector for all points, lpm is not evaluated */
  void GetVector(Vector<Double>& vec, const LocPointMapped& lpm) override;

  void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm)  override {
    EXCEPTION("complex fileData not implemented yet");
  }

  //! Return size of vector in case coefficient function is a vector
  UInt GetVecSize() const override { return dim_; }

  std::string ToString() const override { return filename_; }


private:

  /** common init for the constructor */
  void Init(int dim);

  /** parse the data file. See class description for the data structure. */
  void ReadData(std::istream& input, int dim);

  //pointer to domain
  Domain* domain_;

  //pde name
  std::string pdeName_;

  // file name containing measured data
  std::string filename_;
  
  //dimension 
  UInt dim_;

  /** contains the time values of each time step*/
  std::list<double> time_;

  /** contains the magnetic intensity vector for each time step*/
  std::list<StdVector<Double> > dataH_;

};


} // end of namespace

#endif /* COEFFUNCTIONFILEDATAMEAS_HH_ */
