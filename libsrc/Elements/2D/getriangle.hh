#ifndef FILE_GeTRIANGLE_2002
#define FILE_GeTRIANGLE_2002

#include "baseelem.hh"
#include "jacobian.hh"
 
namespace CoupledField
{

//! Class with general procedures for triangle finite elements   
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on triangles, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights
*/

  class GeTriangle : public BaseElem
{
public:

  //! Constructor
  GeTriangle(enum IntegrationType aintegtype) { IntegType = aintegtype; }
 
  //! Deconstructor
  virtual ~GeTriangle();

protected:
 
  //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate transformation function at these integration points
  virtual void SetTransformFncAtIntPoints();

  //! Calculate derivative of transformation function at these integration points
  virtual void SetDerTransformFncAtIntPoints();
 
  //! Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
                 const Point2D * const ptCoord, const Boolean NeedJinv=TRUE);
 
    //! Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
                  const Point3D * const ptCoord, const Boolean NeedJinv=TRUE);
 
  Vector<Double> TransFncAtIP1;
  Vector<Double> TransFncAtIP2;
  Vector<Double> TransFncAtIP3;
 
  Vector<Double> DxTransFncAtIP1;
  Vector<Double> DxTransFncAtIP2;
  Vector<Double> DxTransFncAtIP3;
 
  Vector<Double> DyTransFncAtIP1;
  Vector<Double> DyTransFncAtIP2;
  Vector<Double> DyTransFncAtIP3;

  Boolean IsSet;

private:

  Double TransFnc3(Double x,Double y) { return y; }
  Double TransFnc2(Double x,Double y) { return x; }
  Double TransFnc1(Double x,Double y) { return 1-x-y; }
 
  Double TransFnc3dx (Double x,Double y) { return 0.0; }
  Double TransFnc3dy (Double x,Double y) { return 1.0; }
  Double TransFnc2dx (Double x,Double y) { return 1.0; }
  Double TransFnc2dy (Double x,Double y) { return 0.0; }
  Double TransFnc1dx (Double x,Double y)  { return -1.0; }
  Double TransFnc1dy (Double x,Double y)  { return -1.0; }

};
 
}
 
#endif // FILE_TRIANGLE1_2001
