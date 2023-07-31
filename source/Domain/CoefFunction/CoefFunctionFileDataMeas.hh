#ifndef COEFFUNCTIONFILEDATAMEAS_HH_
#define COEFFUNCTIONFILEDATAMEAS_HH_

#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"
#include "Domain/Domain.hh"

namespace CoupledField
{

/** This is a more technical variant of CoefFunctionScatteredData. Here we need to give the values
 *  for the cfs element/node numbers - hence the data needs to match the mesh!
 *  Currently we handle only nodal data (named nodes) for scalar and vector, all real. Extensions shall be easy!
 *  Search for "fileData" in the tests for application examples.
 *  The file format is the following:
 *  * ASCII and line based.
 *  * lines starting with hashtag # are ignored
 *  * first value is 1-based node number
 *  * following double values (one for scalar, more for vectors)
 *  * valid separators are space, tab, comma, semicolon and multiples of it and mixing. Hence almost everything
 *  * see also StringParse in testbed.cc  */
class CoefFunctionFileDataMeas : public CoefFunction, public boost::enable_shared_from_this<CoefFunctionFileDataMeas>
{
public:

  /** <fileData> Element, e.g. for mechanic force */
  CoefFunctionFileDataMeas(Domain* ptDomain, const std::string& pdename, PtrParamNode pn, int dim);

  virtual ~CoefFunctionFileDataMeas() { }

  string GetName() const override { return "CoefFunctionFileDataMeas"; }

  void GetScalar(double& scal, const LocPointMapped& lpm) override {
     EXCEPTION("GetScalar not supported");
  }

  void GetScalar(Complex& scal, const LocPointMapped& lpm)  override{
    EXCEPTION("GetScalar not implemented yet");
  }

  /** the node of interest is transported in lpm.pl.number */
  void GetVector(Vector<Double>& vec, const LocPointMapped& lpm);

  void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm)  override {
    EXCEPTION("complex fileData not implemented yet");
  }

  //! Return size of vector in case coefficient function is a vector
  UInt GetVecSize() const { return dim_; }

  std::string ToString() const override { return filename_; }


private:

  /** common init for the constructors */
  void Init(int dim);

  /** parse the data file. See class description for data structure.
   * @param input e.g. an open ifstream or something like a strstream for testing */
  void ReadData(std::istream& input, int dim);

  //pointer to domain
  Domain* domain_;

  //pde name
  std::string pdeName_;

  // file name containing measured data
  std::string filename_;
  
  //dimension 
  UInt dim_;

  /** how to handle node numbers not given in the file */
  typedef enum { THROW_EXCEPTION, WARNING_ZERO, ZERO } Missing;
  Enum<Missing> missing;
  Missing missing_ = THROW_EXCEPTION;

  /** helper for Missing::WARNING_ZERO */
  bool warning_issued_ = false;

  /** contains the time values of each time step*/
  std::list<double> time_;

  /** contains the magnetic intensity vector for each time step*/
  std::list<StdVector<Double> > dataH_;

};


} // end of namespace

#endif /* COEFFUNCTIONFILEDATAMEAS_HH_ */
