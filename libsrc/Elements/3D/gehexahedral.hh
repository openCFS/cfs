#ifndef FILE_GEHEXAHEDRAL_2003
#define FILE_GEHEXAHEDRAL_2003

#include "baseelem.hh"

namespace CoupledField
{

//! Class with general description of hexahedral element
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on hexahedral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights
*/

class GeHexahedral : public BaseElem
{
public:
   //! Constructor with type of integration rule
   GeHexahedral();

   //! Deconstructor
  virtual ~GeHexahedral(); 

  //! return FE-Type for CLA++
  virtual Integer feType() { return HEX;}
	
///////////////////////////////////////////////////////////////////////

virtual Vector<Double> *  GetDXiShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDEtShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDZeShFncAtIP(const Integer iShFnc) ;

protected:

   void test(){ std::cout << " hexahedral " << std::endl;}

   //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate transformation function at these integration points
  virtual void SetTransformFncAtIntPoints();

  //! Calculate derivative of transformation function at these integration points
  virtual void SetDerTransformFncAtIntPoints();
  
   //! Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<2> & J, const Integer ip,
                 Point<2> *  ptCoord, const Boolean NeedJinv=TRUE); 
    
   //! Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<3> & J, const Integer ip,
                  Point<3> *  ptCoord, const Boolean NeedJinv=TRUE);

   //! Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobianAtCenterOfElem(Jacobian<3> & J,
                 Point<3> *  ptCoord, const Boolean NeedJinv=TRUE);

  Vector<Double> TransFnc1AtIP;
  Vector<Double> TransFnc2AtIP;
  Vector<Double> TransFnc3AtIP;
  Vector<Double> TransFnc4AtIP;
  Vector<Double> TransFnc5AtIP;
  Vector<Double> TransFnc6AtIP;
  Vector<Double> TransFnc7AtIP;
  Vector<Double> TransFnc8AtIP;

  Vector<Double> DXiTransFnc1AtIP;
  Vector<Double> DXiTransFnc2AtIP;
  Vector<Double> DXiTransFnc3AtIP;
  Vector<Double> DXiTransFnc4AtIP;
  Vector<Double> DXiTransFnc5AtIP;
  Vector<Double> DXiTransFnc6AtIP;
  Vector<Double> DXiTransFnc7AtIP;
  Vector<Double> DXiTransFnc8AtIP;

  Vector<Double> DEtTransFnc1AtIP;
  Vector<Double> DEtTransFnc2AtIP;
  Vector<Double> DEtTransFnc3AtIP;
  Vector<Double> DEtTransFnc4AtIP;
  Vector<Double> DEtTransFnc5AtIP;
  Vector<Double> DEtTransFnc6AtIP;
  Vector<Double> DEtTransFnc7AtIP;
  Vector<Double> DEtTransFnc8AtIP;

  Vector<Double> DZeTransFnc1AtIP;
  Vector<Double> DZeTransFnc2AtIP;
  Vector<Double> DZeTransFnc3AtIP;
  Vector<Double> DZeTransFnc4AtIP;
  Vector<Double> DZeTransFnc5AtIP;
  Vector<Double> DZeTransFnc6AtIP;
  Vector<Double> DZeTransFnc7AtIP;
  Vector<Double> DZeTransFnc8AtIP;

  Boolean IsSet;

private:

  Double TransFnc1(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1-et)*(1-ze)); }
  Double TransFnc2(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1-et)*(1-ze)); }
  Double TransFnc3(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1+et)*(1-ze)); }
  Double TransFnc4(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1+et)*(1-ze)); }
  Double TransFnc5(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1-et)*(1+ze)); }
  Double TransFnc6(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1-et)*(1+ze)); }
  Double TransFnc7(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1+et)*(1+ze)); }
  Double TransFnc8(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1+et)*(1+ze)); }

  Double TransFnc1dxi (Double xi,Double et, Double ze) { return -0.125*(1-et)*(1-ze); }
  Double TransFnc1deta (Double xi,Double et, Double ze) { return -0.125*(1-xi)*(1-ze); }
  Double TransFnc1dzeta (Double xi,Double et, Double ze) { return -0.125*(1-xi)*(1-et); }

  Double TransFnc2dxi (Double xi,Double et, Double ze) { return 0.125*(1-et)*(1-ze); }
  Double TransFnc2deta (Double xi,Double et, Double ze) { return -0.125*(1+xi)*(1-ze); }
  Double TransFnc2dzeta (Double xi,Double et, Double ze) { return -0.125*(1+xi)*(1-et); }

  Double TransFnc3dxi (Double xi,Double et, Double ze)  { return 0.125*(1+et)*(1-ze); }
  Double TransFnc3deta (Double xi,Double et, Double ze)  { return 0.125*(1+xi)*(1-ze); }
  Double TransFnc3dzeta (Double xi,Double et, Double ze)  { return -0.125*(1+xi)*(1+et); }

  Double TransFnc4dxi (Double xi,Double et, Double ze)  { return -0.125*(1+et)*(1-ze); }
  Double TransFnc4deta (Double xi,Double et, Double ze)  { return 0.125*(1-xi)*(1-ze); }
  Double TransFnc4dzeta (Double xi,Double et, Double ze)  { return -0.125*(1-xi)*(1+et); }

  Double TransFnc5dxi (Double xi,Double et, Double ze) { return -0.125*(1-et)*(1+ze); }
  Double TransFnc5deta (Double xi,Double et, Double ze) { return -0.125*(1-xi)*(1+ze); }
  Double TransFnc5dzeta (Double xi,Double et, Double ze) { return 0.125*(1-xi)*(1-et); }

  Double TransFnc6dxi (Double xi,Double et, Double ze) { return 0.125*(1-et)*(1+ze); }
  Double TransFnc6deta (Double xi,Double et, Double ze) { return -0.125*(1+xi)*(1+ze); }
  Double TransFnc6dzeta (Double xi,Double et, Double ze) { return 0.125*(1+xi)*(1-et); }

  Double TransFnc7dxi (Double xi,Double et, Double ze)  { return 0.125*(1+et)*(1+ze); }
  Double TransFnc7deta (Double xi,Double et, Double ze)  { return 0.125*(1+xi)*(1+ze); }
  Double TransFnc7dzeta (Double xi,Double et, Double ze)  { return 0.125*(1+xi)*(1+et); }

  Double TransFnc8dxi (Double xi,Double et, Double ze)  { return -0.125*(1+et)*(1+ze); }
  Double TransFnc8deta (Double xi,Double et, Double ze)  { return 0.125*(1-xi)*(1+ze); }
  Double TransFnc8dzeta (Double xi,Double et, Double ze)  { return 0.125*(1-xi)*(1+et); }
};

}
#endif //
