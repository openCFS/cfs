#ifndef FILE_BASEFORM_
#define FILE_BASEFORM_

#include <Elements/basefe.hh>
#include <DataInOut/MaterialData.hh>

namespace CoupledField
{

//! Base class for calculation of bdb element matrices
class BaseForm 
{
public:
  //! Constructor
  BaseForm(BaseFE * aptelem, MaterialData & matData);

  //! Constructor
  BaseForm(BaseFE * aptelem);

  //! Deconstructor
  virtual ~BaseForm();

  //! Virtual function, implemented in derived classes
  virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & StiffMat);

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & result)
  {Error("CalcElemVector not implemented!",__FILE__,__LINE__);};

  //! Prints the bilinear form
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;


protected:

  //! Ptr to base element
  BaseFE  * ptelem;

  /// Ptr to material
  MaterialData * ptMaterial ;
};

} //end namespace

#endif // FILE_BASEFORM
