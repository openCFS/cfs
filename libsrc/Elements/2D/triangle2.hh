#ifndef FILE_TRIANGLE2_2003
#define FILE_TRIANGLE2_2003

#include "getriangle.hh"
 
namespace CoupledField
{

//! Description of triangle element with 2 order interpolation 

  class Triangle2 : public GeTriangle
{
public:

  //! Constructor with type of integration rule
  Triangle2();
 
  //! Deconstructor
  virtual ~Triangle2();
 
  //! Define variables of this class
  virtual void Init();

  //! Return Vector<Double> of gradient of shape func with number i at integration point ip
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

  //! get gradient of shape fnc for node i at center of element
  virtual  void GetGradientShFncAtCenter(Vector<Double> & ,const Integer i);

  //! Return pointer to array of value shape function at integration points
  virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

  //! test type of element
  virtual enum ElementType type(){ return Triang2;}

protected:

private:

  void SetShapeFncAtIntPoints();
  void SetDerShapeFncAtIntPoints();  
  
  Double ShapeFnc1(Double x, Double y) { return 2.*(1.-x-y)*(0.5-x-y); }
  Double DxShapeFnc1(Double x, Double y) { return 4.*(x+y)-3.; }
  Double DyShapeFnc1(Double x, Double y) { return 4.*(x+y)-3.; }
  Double ShapeFnc2(Double x, Double y)   { return x*(2.*x-1.); }
  Double DxShapeFnc2(Double x, Double y) { return 4.*x-1.; }
  Double DyShapeFnc2(Double x, Double y) { return 0.; }
  Double ShapeFnc3(Double x, Double y)   { return y*(2.*y-1.); }
  Double DxShapeFnc3(Double x, Double y) { return 0.; }
  Double DyShapeFnc3(Double x, Double y) { return 4.*y-1.; }
  Double ShapeFnc4(Double x, Double y)   { return 4.*x*y; }
  Double DxShapeFnc4(Double x, Double y) { return 4.*y; }
  Double DyShapeFnc4(Double x, Double y) { return 4.*x; }
  Double ShapeFnc5(Double x, Double y)   { return 4.*y*(1.-x-y); }
  Double DxShapeFnc5(Double x, Double y) { return -4.*y; }
  Double DyShapeFnc5(Double x, Double y) { return 4.*(1-x-2.*y); }
  Double ShapeFnc6(Double x, Double y)   { return 4.*x*(1.-x-y); }
  Double DxShapeFnc6(Double x, Double y) { return 4.*(1-2.*x-y); }
  Double DyShapeFnc6(Double x, Double y) { return -4.*x; }
  
   Vector<Double> ShapeFncAtIP1;
  Vector<Double> ShapeFncAtIP2; 
  Vector<Double> ShapeFncAtIP3; 
  Vector<Double> ShapeFncAtIP4; 
  Vector<Double> ShapeFncAtIP5; 
  Vector<Double> ShapeFncAtIP6; 

  Vector<Double> DxShapeFncAtIP1;
  Vector<Double> DxShapeFncAtIP2;
  Vector<Double> DxShapeFncAtIP3;
  Vector<Double> DxShapeFncAtIP4;
  Vector<Double> DxShapeFncAtIP5;
  Vector<Double> DxShapeFncAtIP6;

  Vector<Double> DyShapeFncAtIP1;
  Vector<Double> DyShapeFncAtIP2;
  Vector<Double> DyShapeFncAtIP3;
  Vector<Double> DyShapeFncAtIP4;
  Vector<Double> DyShapeFncAtIP5;
  Vector<Double> DyShapeFncAtIP6;
 
};
 
}
 
#endif // FILE_TRIANGLE2_2003
