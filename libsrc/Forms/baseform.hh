#ifndef FILE_BASEFORM_2001
#define FILE_BASEFORM_2001

#include "baseelem.hh"

namespace CoupledField
{

typedef Double (*RHS)(Double,Double,Double);

/// Base class for calculation of elements matrixes
template<Integer dim>
class BaseForm 
{
public:
  /// Constructor
  BaseForm(BaseElem * aptelem);

  /// Deconstructor
  virtual ~BaseForm();

  /// Virtual function, implemented in derived classes
  virtual void CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> & StiffMat);

  /// Virtual function, implemented in derived classes
  virtual void CalcElemVector(Point<dim> * ptCoord, Vector<Double> & StiffMat);

   /// Virtual function, implemented in derived classes
  virtual void CalcElemVector4InterpolatedFnc(Point<dim> * ptCoord, const Integer aComponent, Vector<Double> & aValueAtNodePoints, Vector<Double> & Result);

  /// Calculate element matrix for multiplication of derivatives respect to x of shape functions
  void CalcShapeFncDxShapeFncDx(Point<dim> * ptCoord, Matrix<Double> & Result);

  /// Calculate element matrix for multiplication of derivatives respect to y of shape functions
  void CalcShapeFncDyShapeFncDy(Point<dim> * ptCoord, Matrix<Double> & Result);

  /// Calculate element matrix for multiplication of shape functions
  void CalcShapeFncShapeFnc(Point<dim> * ptCoord, Matrix<Double> & Result);

  ///
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:

  /// 
  BaseElem  * ptelem;

  /// Calculation of function in integration points
  Double FuncAtIP(const ShortInt iIP, RHS f);  


//////////////////////////////////////////////////////////////// version 2
  /// calculation of determinat for 2 version, only for 3d
  Double CalcDet(const Integer iIntPoints, Point<3> * ptCoord);  

   /// Calculation of tramsformation coordinates in integration points
  Double TransXDxi(Integer iIntPoints, Point<3> * ptCoord) ;
  /// Calculation of tramsformation coordinates in integration points
  Double TransXDnu(Integer iIntPoints, Point<3> * ptCoord);
    /// Calculation of tramsformation coordinates in integration points
  Double TransXDgam(Integer iIntPoints, Point<3> * ptCoord);

  /// Calculation of tramsformation coordinates in integration points
  Double TransYDxi(Integer iIntPoints, Point<3> * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransYDnu (Integer iIntPoints, Point<3> * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransYDgam (Integer iIntPoints, Point<3> * ptCoord);

   /// Calculation of tramsformation coordinates in integration points
  Double TransZDxi(Integer iIntPoints, Point<3> * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransZDnu (Integer iIntPoints, Point<3> * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransZDgam (Integer iIntPoints, Point<3> * ptCoord);

};

#ifdef __GNUC__
template class BaseForm<3>;
template class BaseForm<2>;
#endif

} //end namespace

#endif // FILE_BASERFORM
