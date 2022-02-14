#ifndef COEFFUNCTIONFILEDATA_HH_
#define COEFFUNCTIONFILEDATA_HH_

#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"

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
class CoefFunctionFileData : public CoefFunction, public boost::enable_shared_from_this<CoefFunctionFileData>
{
public:

  /** <fileData> Element, e.g. for mechanic force */
  CoefFunctionFileData(PtrParamNode scatteredDataNode, int dim);

  /** for testing purpose only */
  CoefFunctionFileData(std::istream& input, int dim);

  virtual ~CoefFunctionFileData() { }

  string GetName() const override { return "CoefFunctionFileData"; }

  void GetScalar(double& scal, const LocPointMapped& lpm) override {
    assert(data_.GetNumCols() == 1);
    data_.GetEntry(GetIndex(lpm),0, scal);
  }

  void GetScalar(Complex& scal, const LocPointMapped& lpm)  override{
    EXCEPTION("complex fileData not implemented yet");
  }

  /** the node of interest is transported in lpm.pl.number */
  void GetVector(Vector<double>& vec, const LocPointMapped& lpm)  override {
    data_.GetRow(vec, GetIndex(lpm));
  }

  void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm)  override {
    EXCEPTION("complex fileData not implemented yet");
  }

  unsigned int GetVecSize() const override { return data_.GetNumCols(); }

  std::string ToString() const override { return filename_; }

  /** as rhs integrator we don't want to be normalized by number of nodes. */
  bool DoNormalize() const override { return false; }

private:

  /** common init for the constructors */
  void Init(int dim);

  /** parse the data file. See class description for data structure.
   * @param input e.g. an open ifstream or something like a strstream for testing */
  void ReadData(std::istream& input, int dim);

  /** common helper for GetScalar/GetVector */
  unsigned int GetIndex(const LocPointMapped& lpm);


  std::string filename_;

  /** map of node numbers */
  StdVector<unsigned int> node_;

  /** this keeps the guess for fast Find in node_ */
  CfsTLS<unsigned int> node_guess_;

  /** rows is number of nodes, cols is dim (1 for scalar, 2/3 for vector, tensor would be flat if one imlements it*/
  Matrix<double> data_;


};


} // end of namespace

#endif /* COEFFUNCTIONFILEDATA_HH_ */
