#ifndef FILE_RECTANGLE_2002
#define FILE_RECTANGLE_2002

#include "baseelem.hh"
#include "jacobian.hh"

namespace CoupledField
{
//! Class with general procedures for quadrilateral finite elements
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights 
*/
  
class Rectangle : public BaseElem
{
public:

  //! Constructor with type of integration rule
  Rectangle(); 
  
  //! Deconstructor
  virtual ~Rectangle();

  //! return FE-Type for CLA++
  virtual Integer feType() { return QUAD;}

protected:

  //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate transformation function at these intergration points
  virtual void SetTransformFncAtIntPoints();

     //! Calculate derivative of transformation function at these intergration points
  virtual void SetDerTransformFncAtIntPoints();
 
  //! Calculate transformation function at the center of element
  virtual void SetTransformFncAtCenter();

     //! Calculate derivative of transformation function at the center of element
  virtual void SetDerTransformFncAtCenter();
  
  //! Calculation of Jacobian, inverse Jacobian and detJacobian
  virtual void CalcJacobian(Jacobian<2> & J, const Integer ip,
                 Point<2> * ptCoord, const Boolean NeedJinv=TRUE);

  //! Calculation of Jacobian, inverse Jacobian and detJacobian
  virtual void CalcJacobian(Jacobian<3> & J, const Integer ip,
                  Point<3> * ptCoord, const Boolean NeedJinv=TRUE);

 //! Calculation of Jacobian for center point in 2D
  virtual void CalcJacobianAtCenter(Jacobian<2> & J, Point<2> * ptCoord, const Boolean NeedJinv=TRUE);
  
  Vector<Double> TransFncAtIP1;
  Vector<Double> TransFncAtIP2;
  Vector<Double> TransFncAtIP3;
  Vector<Double> TransFncAtIP4;
 
  Vector<Double> DxTransFncAtIP1;
  Vector<Double> DxTransFncAtIP2;
  Vector<Double> DxTransFncAtIP3;
  Vector<Double> DxTransFncAtIP4;
 
  Vector<Double> DyTransFncAtIP1;
  Vector<Double> DyTransFncAtIP2;
  Vector<Double> DyTransFncAtIP3;
  Vector<Double> DyTransFncAtIP4;

  Double TransFncAtCenter[4];
  Double DxTransFncAtCenter[4];
  Double DyTransFncAtCenter[4];

  Boolean IsSet; //!< parameter is TRUE if we calculated transformation function at integration points
  Boolean isSetAtCenter_; //!< parameter is TRUE if we calculated tranformation function at integration points

private:
  Double TransFnc1(Double x,Double y) {return 0.25*(1+x)*(1+y); }
  Double TransFnc2(Double x,Double y) {return 0.25*(1-x)*(1+y); }
  Double TransFnc3(Double x,Double y) {return 0.25*(1-x)*(1-y); }
  Double TransFnc4(Double x,Double y) {return 0.25*(1+x)*(1-y); }

  Double TransFnc1dx (Double x,Double y)  { return 0.25*(1+y);} 
  Double TransFnc1dy (Double x,Double y)  { return 0.25*(1+x);} 
  Double TransFnc2dx (Double x,Double y)  { return 0.25*(-1-y);} 
  Double TransFnc2dy (Double x,Double y)  { return 0.25*(1-x);} 
  Double TransFnc3dx (Double x,Double y)  { return 0.25*(y-1);}
  Double TransFnc3dy (Double x,Double y)  { return 0.25*(x-1);}
  Double TransFnc4dx (Double x,Double y)  { return 0.25*(1-y);} 
  Double TransFnc4dy (Double x,Double y)  { return 0.25*(-1-x);} 
};

}

#endif // FILE_RECTANGLE
