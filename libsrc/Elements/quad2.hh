#ifndef FILE_DQUADRILATERAL2_2002
#define FILE_DQUADRILATERAL2_2002

namespace CoupledField
{
//! Quadrilateral finite element with four nodes (linear interpolation function)

  class Quad2 : public Rectangle
{
public:

  //! Constructor with type of integration rule
  Quad2(enum IntegrationType aintegtype);
  
  //! Deconstructor
  virtual ~Quad2();

  //! Define variables of this class
  virtual void Init();

  //! Return Vector<Double> of gradient of shape func with number i at integration point ip 
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);  
 
  //! Return pointer to array of value shape function at integration points
  virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

protected:

private:
  ShortInt ElemType;  //!< element type

  void SetShapeFncAtIntPoints();
  void SetDerShapeFncAtIntPoints();  

  Double ShapeFnc1(Double x, Double y) { return 0.25*x*y*(x-1)*(y-1);}
  Double ShapeFnc2(Double x, Double y) { return -0.25*x*y*(x+1)*(y-1);}
  Double ShapeFnc3(Double x, Double y) { return 0.25*x*y*(x+1)*(y+1);}
  Double ShapeFnc4(Double x, Double y) { return -0.25*x*y*(x-1)*(y+1);}
  Double ShapeFnc5(Double x, Double y) { return 0.5*y*(1-x*x)*(y-1);}
  Double ShapeFnc6(Double x, Double y) { return 0.5*x*(x+1)*(1-y*y);}
  Double ShapeFnc7(Double x, Double y) { return 0.5*y*(1-x*x)*(y+1);}
  Double ShapeFnc8(Double x, Double y) { return 0.5*x*(x-1)*(1-y*y);}
  Double ShapeFnc9(Double x, Double y) { return (1-x*x)*(1-y*y);}

  Double DxShapeFnc1(Double x, Double y) { return 0.25*y*(2*x-1)*(y-1);}
  Double DxShapeFnc2(Double x, Double y) { return -0.25*y*(2*x+1)*(y-1);}
  Double DxShapeFnc3(Double x, Double y) { return 0.25*y*(2*x+1)*(y+1);}
  Double DxShapeFnc4(Double x, Double y) { return -0.25*y*(2*x-1)*(y+1);}
  Double DxShapeFnc5(Double x, Double y) { return 0.5*y*(-x*2)*(y-1);}
  Double DxShapeFnc6(Double x, Double y) { return 0.5*(2*x+1)*(1-y*y);}
  Double DxShapeFnc7(Double x, Double y) { return 0.5*y*(-2*x)*(y+1);}
  Double DxShapeFnc8(Double x, Double y) { return 0.5*(2*x-1)*(1-y*y);}
  Double DxShapeFnc9(Double x, Double y) { return (-x*2)*(1-y*y);}

  Double DyShapeFnc1(Double x, Double y) { return 0.25*x*(x-1)*(2*y-1);}
  Double DyShapeFnc2(Double x, Double y) { return -0.25*x*(x+1)*(2*y-1);}
  Double DyShapeFnc3(Double x, Double y) { return 0.25*x*(x+1)*(2*y+1);}
  Double DyShapeFnc4(Double x, Double y) { return -0.25*x*(x-1)*(2*y+1);}
  Double DyShapeFnc5(Double x, Double y) { return 0.5*(1-x*x)*(2*y-1);}
  Double DyShapeFnc6(Double x, Double y) { return 0.5*x*(x+1)*(-y*2);}
  Double DyShapeFnc7(Double x, Double y) { return 0.5*(1-x*x)*(2*y+1);}
  Double DyShapeFnc8(Double x, Double y) { return 0.5*x*(x-1)*(-y*2);}
  Double DyShapeFnc9(Double x, Double y) { return (1-x*x)*(-y*2);}

  Vector<Double> ShapeFncAtIP1;
  Vector<Double> ShapeFncAtIP2; 
  Vector<Double> ShapeFncAtIP3; 
  Vector<Double> ShapeFncAtIP4; 
  Vector<Double> ShapeFncAtIP5; 
  Vector<Double> ShapeFncAtIP6; 
  Vector<Double> ShapeFncAtIP7; 
  Vector<Double> ShapeFncAtIP8; 
  Vector<Double> ShapeFncAtIP9; 

  Vector<Double> DxShapeFncAtIP1;
  Vector<Double> DxShapeFncAtIP2;
  Vector<Double> DxShapeFncAtIP3;
  Vector<Double> DxShapeFncAtIP4;
  Vector<Double> DxShapeFncAtIP5;
  Vector<Double> DxShapeFncAtIP6;
  Vector<Double> DxShapeFncAtIP7;
  Vector<Double> DxShapeFncAtIP8;
  Vector<Double> DxShapeFncAtIP9;

  Vector<Double> DyShapeFncAtIP1;
  Vector<Double> DyShapeFncAtIP2;
  Vector<Double> DyShapeFncAtIP3;
  Vector<Double> DyShapeFncAtIP4;
  Vector<Double> DyShapeFncAtIP5;
  Vector<Double> DyShapeFncAtIP6;
  Vector<Double> DyShapeFncAtIP7;
  Vector<Double> DyShapeFncAtIP8;
  Vector<Double> DyShapeFncAtIP9;

};

}

#endif // FILE_DQUADRILATERAL1
