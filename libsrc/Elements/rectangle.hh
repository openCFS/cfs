#ifndef FILE_RECTANGLE_2002
#define FILE_RECTANGLE_2002

namespace CoupledField
{
//! Class with general procedures for quadrilateral finite elements
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights 
*/
  
class Rectangle : public BaseElem
{
public:

  //! Constructor with type of integration rule
  Rectangle(enum IntegrationType aintegtype){ IntegType = aintegtype; IsSet=FALSE; }
  
  //! Deconstructor
  virtual ~Rectangle();

protected:

  //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate transformation function at these intergration points
  virtual void SetTransformFncAtIntPoints();

     //! Calculate derivative of transformation function at these intergration points
  virtual void SetDerTransformFncAtIntPoints();
 
  // Calculation of Jacobian, inverse Jacobian and detJacobian
  virtual void CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
                 const Point2D * const ptCoord, const Boolean NeedJinv=TRUE);

    // Calculation of Jacobian, inverse Jacobian and detJacobian
  virtual void CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
                  const Point3D * const ptCoord, const Boolean NeedJinv=TRUE);

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

  Boolean IsSet; //!< parameter is TRUE if we calculated transformation function at integration points
private:

  Double TransFnc3(Double x,Double y) {return 0.25*(1-x)*(1-y); }
  Double TransFnc2(Double x,Double y) {return 0.25*(1-x)*(1+y); }
  Double TransFnc1(Double x,Double y) {return 0.25*(1+x)*(1+y); }
  Double TransFnc4(Double x,Double y) {return 0.25*(1+x)*(1-y); }

  Double TransFnc3dx (Double x,Double y) { return 0.25*(y-1);}
  Double TransFnc3dy (Double x,Double y) { return 0.25*(x-1);}
  Double TransFnc2dx (Double x,Double y) { return 0.25*(-1-y);} 
  Double TransFnc2dy (Double x,Double y) { return 0.25*(1-x);} 
  Double TransFnc1dx (Double x,Double y)  { return 0.25*(1+y);} 
  Double TransFnc1dy (Double x,Double y)  { return 0.25*(1+x);} 
  Double TransFnc4dx (Double x,Double y)  { return 0.25*(1-y);} 
  Double TransFnc4dy (Double x,Double y)  { return 0.25*(-1-x);} 

//  ShortInt IntegType; //!< integration rule 

};

}

#endif // FILE_RECTANGLE
