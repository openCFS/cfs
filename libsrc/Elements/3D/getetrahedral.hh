#ifndef FILE_GETETRAHEDRAL_2002
#define FILE_GETETRAHEDRAL_2002

#include "baseelem.hh"

namespace CoupledField
{

//! Class with general description of tetrahedral element
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on tetrahedral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights
*/

class GeTetrahedral : public BaseElem
{
public:
   //! Constructor with type of integration rule
   GeTetrahedral();

   //! Deconstructor
  virtual ~GeTetrahedral(); 

  //! return FE-Type for CLA++
  virtual Integer feType() { return TET;}
  
///////////////////////////////////////////////////////////////////////

virtual Vector<Double> *  GetDxShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDyShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDzShFncAtIP(const Integer iShFnc) ;

protected:

   void test(){ std::cout << " tetrahedral " << std::endl;}

   //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate transformation function at these integration points
  virtual void SetTransformFncAtIntPoints();

  //! Calculate derivative of transformation function at these integration points
  virtual void SetDerTransformFncAtIntPoints();
  
  //! Calculate transformation function at the center of element
  virtual void SetTransformFncAtCenter();

  //! Calculate derivative of transformation function at the center of element
  virtual void SetDerTransformFncAtCenter();

   //! Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<2> & J, const Integer ip,
                 Point<2> * ptCoord, const Boolean NeedJinv=TRUE); 
    
   //! Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<3> & J, const Integer ip,
                   Point<3> * ptCoord, const Boolean NeedJinv=TRUE);

 //! Calculation of Jacobian for center point in 3D
  virtual void CalcJacobianAtCenter(Jacobian<3> & J,  Point<3> *  ptCoord, const Boolean NeedJinv=TRUE);

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

  Vector<Double> DzTransFncAtIP1;
  Vector<Double> DzTransFncAtIP2;
  Vector<Double> DzTransFncAtIP3;
  Vector<Double> DzTransFncAtIP4;

  Double TransFncAtCenter[4];
  Double DxTransFncAtCenter[4];
  Double DyTransFncAtCenter[4];
  Double DzTransFncAtCenter[4];

  Boolean IsSet;
  Boolean isSetAtCenter_;

private:

  Double TransFnc1(Double x,Double y, Double z) { return 1-x-y-z; }
  Double TransFnc2(Double x,Double y, Double z) { return x; }
  Double TransFnc3(Double x,Double y, Double z) { return y; }
  Double TransFnc4(Double x,Double y, Double z) { return z; }

  Double TransFnc1dx (Double x,Double y, Double z) { return -1.0; }
  Double TransFnc1dy (Double x,Double y, Double z) { return -1.0; }
  Double TransFnc1dz (Double x,Double y, Double z) { return -1.0; }

  Double TransFnc2dx (Double x,Double y, Double z) { return 1.0; }
  Double TransFnc2dy (Double x,Double y, Double z) { return 0.0; }
  Double TransFnc2dz (Double x,Double y, Double z) { return 0.0; }

  Double TransFnc3dx (Double x,Double y, Double z)  { return 0.0; }
  Double TransFnc3dy (Double x,Double y, Double z)  { return 1.0; }
  Double TransFnc3dz (Double x,Double y, Double z)  { return 0.0; }

  Double TransFnc4dx (Double x,Double y, Double z)  { return 0.0; }
  Double TransFnc4dy (Double x,Double y, Double z)  { return 0.0; }
  Double TransFnc4dz (Double x,Double y, Double z)  { return 1.0; }
};

}
#endif //
