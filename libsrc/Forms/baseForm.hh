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
  BaseForm(MaterialData & matData);

  //! Constructor
  BaseForm(BaseFE * aptelem);

  //! Constructor
  BaseForm();

  //! Deconstructor
  virtual ~BaseForm();

  //! Virtual function, implemented in derived classes
  virtual void CalcElementMatrix(Matrix<Double>& ptCoord, Matrix<Double> & StiffMat);

  /// Calculation of vector of right hand side 
  virtual void CalcElemVector(Matrix<Double>& ptCoord, std::vector<Double> & result)
  {Error("CalcElemVector not implemented!",__FILE__,__LINE__);};

  virtual void CalcElemVector4Dip(Matrix<Double>& ptCoord, const Vector<Integer> & connecth, 
				  std::vector<Double> & Result, const std::vector<Double> gradN_x_P)
  { Error(" CalcElemVector4Dip is not implemented for this class",__FILE__,__LINE__);}

  /// Calculation of vector of right hand side given from quadrupole contribution
  virtual void CalcElemVector4Quad(Matrix<Double>& ptCoord,const Vector<Integer> & connecth,
				   const Matrix<Double> & FlowData, std::vector<Double> & Result)
  { Error(" CalcElemVector4Quad is not implemented for this class",__FILE__,__LINE__);}

  /// Extraction of element velocity values from total flowdata matrix to a matrix (connecth, dim)
  virtual void GetQttiesOfElement(Matrix<Double>& elVec,const Matrix<Double>& FlowData,
				  const Vector<Integer>& connecth, Integer matrixRow)
  { Error(" GetQttiesOfElement is not implemented for this class",__FILE__,__LINE__);} 

  //! Prints the bilinear form
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;


  //! sets pointer to actual element
  void SetElemPtr(BaseFE * elemPtr){ptelem = elemPtr;};


  //! sets pointer to actual material
  void SetMaterial(MaterialData * matPtr){ptMaterial = matPtr;};


  //! sets actual element solution
  virtual void SetActElemSol(Matrix<Double>& disp)
  {Error("SetActElemSol not implemented!",__FILE__,__LINE__);};

  //! reads the values y(x) out of the file with name fncName  
  void ReadNlinFunc(std::string fncName, std::vector<double> &xval, std::vector<Double> &yval);


protected:

  //! Ptr to base element
  BaseFE  * ptelem;

  /// Ptr to material
  MaterialData * ptMaterial ;

  Boolean isaxi_;
};

} //end namespace

#endif // FILE_BASEFORM
