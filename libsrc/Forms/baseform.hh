#ifndef FILE_BASEFORM_2001
#define FILE_BASEFORM_2001

#include "baseelem.hh"

namespace CoupledField
{

typedef Double (*RHS)(Double,Double,Double);

/// Base class for calculation of elements matrixes
template<class Dim>
class BaseForm 
{
public:
  /// Constructor
  BaseForm(BaseElem * aptelem);

  /// Deconstructor
  virtual ~BaseForm();

  /// Virtual function, implemented in derived classes
  virtual void CalcElemMatrix(Dim * ptCoord, Matrix<Double> & StiffMat);

  /// Calculate element matrix for multiplication of derivatives respect to x of shape functions
  void CalcShapeFncDxShapeFncDx(Dim * ptCoord, Matrix<Double> & Result);

  /// Calculate element matrix for multiplication of derivatives respect to y of shape functions
  void CalcShapeFncDyShapeFncDy(Dim * ptCoord, Matrix<Double> & Result);

  /// Calculate element matrix for multiplication of shape functions
  void CalcShapeFncShapeFnc(Dim * ptCoord, Matrix<Double> & Result);
  ///
  //  virtual void SetMaterial(Material *material);
  ///
  virtual void Print(std::ostream * out, const Matrix<Double> Result) const;

protected:

  /// 
  BaseElem  * ptelem;

  /// Calculation of function in integration points
  Double FuncAtIP(const ShortInt iIP, RHS f);  


//////////////////////////////////////////////////////////////// version 2
  /// calculation of determinat for 2 version, only for 3d
  Double CalcDet(const Integer iIntPoints, Point3D * ptCoord);  

   /// Calculation of tramsformation coordinates in integration points
  Double TransXDxi(Integer iIntPoints, Point3D * ptCoord) ;
  /// Calculation of tramsformation coordinates in integration points
  Double TransXDnu(Integer iIntPoints, Point3D * ptCoord);
    /// Calculation of tramsformation coordinates in integration points
  Double TransXDgam(Integer iIntPoints, Point3D * ptCoord);

  /// Calculation of tramsformation coordinates in integration points
  Double TransYDxi(Integer iIntPoints, Point3D * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransYDnu (Integer iIntPoints, Point3D * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransYDgam (Integer iIntPoints, Point3D * ptCoord);

   /// Calculation of tramsformation coordinates in integration points
  Double TransZDxi(Integer iIntPoints, Point3D * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransZDnu (Integer iIntPoints, Point3D * ptCoord);
 /// Calculation of tramsformation coordinates in integration points
  Double TransZDgam (Integer iIntPoints, Point3D * ptCoord);

};

#ifdef __GNUC__
template class BaseForm<Point3D>;
template class BaseForm<Point2D>;
#endif

//template void BaseForm<Point2D>:: FuncAtIP<RHSfunction>(const ShortInt, RHSfunction);
//template void BaseForm<Point3D>:: FuncAtIP<RHSfunction>(const ShortInt, RHSfunction);

} //end namespace

#endif // FILE_BASERFORM
