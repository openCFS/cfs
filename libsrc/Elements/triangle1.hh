#ifndef FILE_TRIANGLE1_2001
#define FILE_TRIANGLE1_2001
 
namespace CoupledField
{
  class BaseElem;

/// Description of triangle element with 1 order interpolation 

  class Triangle1 : public BaseElem
{
public:

  //! Constructor
  Triangle1(ShortInt aintegtype);
 
  //! Deconstructor
  virtual ~Triangle1();
 
  //! Define variables of this class
  virtual void Init();

  //! Define integration points for this type of intergration
  virtual void SetIntPoints();

  //! Calculate Shape function at these intergration points
  virtual void SetShapeFncAtIntPoints();

  //! Calculate derivative of shape function at these intergration points
  virtual void SetDShapeFncAtIntPoints();
 
  //! Put information about used method(intergration type, element) in out
  // virtual void Info(ostream * out);
 
  // Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
                 const Point2D * const ptCoord, const Boolean NeedJinv=TRUE);
 
    // Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
                  const Point3D * const ptCoord, const Boolean NeedJinv=TRUE);
 
  //! Return Vector<Double> of gradient of shape func with number i at integration point ip
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

  //! Return pointer to array of value shape function at integration points
  virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

protected:

private:
  ShortInt IntegType;  //!< integration rule
  ShortInt ElemType;  //!< element type

  Vector<Double> ShFncAtIP1;
  Vector<Double> ShFncAtIP2;
  Vector<Double> ShFncAtIP3;

  Vector<Double> DxShFncAtIP1;
  Vector<Double> DxShFncAtIP2;
  Vector<Double> DxShFncAtIP3;

  Vector<Double> DyShFncAtIP1;
  Vector<Double> DyShFncAtIP2;
  Vector<Double> DyShFncAtIP3;

  Double ShapeFnc3(Double x,Double y) {return y; }
  Double ShapeFnc2(Double x,Double y) {return x; }
  Double ShapeFnc1(Double x,Double y) {return 1-x-y; }

  Double ShapeFnc3dx (Double x,Double y) { return 0.0;}
  Double ShapeFnc3dy (Double x,Double y) {return 1.0;}
  Double ShapeFnc2dx (Double x,Double y) { return 1.0;}
  Double ShapeFnc2dy (Double x,Double y) { return 0.0;}
  Double ShapeFnc1dx (Double x,Double y)  { return -1.0;}
  Double ShapeFnc1dy (Double x,Double y)  { return -1.0;}
 
};
 
}
 
#endif // FILE_TRIANGLE1_2001
