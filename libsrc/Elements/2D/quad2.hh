#ifndef FILE_DQUADRILATERAL2_2002
#define FILE_DQUADRILATERAL2_2002

#include "rectangle.hh"

namespace CoupledField
{
//! Quadrilateral finite element with nine nodes (linear interpolation function)

  class Quad2 : public Rectangle
{
public:

  //! Constructor with type of integration rule
  Quad2();
  
  //! Deconstructor
  virtual ~Quad2();

  //! Define variables of this class
  virtual void Init();

  //! Return Vector<Double> of gradient of shape func with number i at integration point ip 
  /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(ip)\ , \ N_{i,\eta}(ip))\f$
     \param i (input) shape function number
     \param ip (input) integration point
  */
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);  
 
  //! Return pointer to array of value shape function at integration points
  /*!
     \param ShFnc (input) shape function number
  */
  virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

  //! test type of element
  virtual enum ElementType type(){ return Quadrilateral2;}

protected:

private:

  void SetShapeFncAtIntPoints();
  void SetDerShapeFncAtIntPoints();  

  Double ShapeFnc1(Double x, Double y) { return 0.25*x*y*(x-1)*(y-1);}   //!< \f$ N_1 \f$
  Double ShapeFnc2(Double x, Double y) { return -0.25*x*y*(x+1)*(y-1);}  //!< \f$ N_2 \f$
  Double ShapeFnc3(Double x, Double y) { return 0.25*x*y*(x+1)*(y+1);}   //!< \f$ N_3 \f$
  Double ShapeFnc4(Double x, Double y) { return -0.25*x*y*(x-1)*(y+1);}  //!< \f$ N_4 \f$
  Double ShapeFnc5(Double x, Double y) { return 0.5*y*(1-x*x)*(y-1);}    //!< \f$ N_5 \f$
  Double ShapeFnc6(Double x, Double y) { return 0.5*x*(x+1)*(1-y*y);}    //!< \f$ N_6 \f$
  Double ShapeFnc7(Double x, Double y) { return 0.5*y*(1-x*x)*(y+1);}    //!< \f$ N_7 \f$
  Double ShapeFnc8(Double x, Double y) { return 0.5*x*(x-1)*(1-y*y);}    //!< \f$ N_8 \f$
  Double ShapeFnc9(Double x, Double y) { return (1-x*x)*(1-y*y);}        //!< \f$ N_9 \f$

  Double DxShapeFnc1(Double x, Double y) { return 0.25*y*(2*x-1)*(y-1);}  //!< \f$ N_{1,\xi} \f$
  Double DxShapeFnc2(Double x, Double y) { return -0.25*y*(2*x+1)*(y-1);} //!< \f$ N_{2,\xi} \f$
  Double DxShapeFnc3(Double x, Double y) { return 0.25*y*(2*x+1)*(y+1);}  //!< \f$ N_{3,\xi} \f$
  Double DxShapeFnc4(Double x, Double y) { return -0.25*y*(2*x-1)*(y+1);} //!< \f$ N_{4,\xi} \f$
  Double DxShapeFnc5(Double x, Double y) { return 0.5*y*(-x*2)*(y-1);}    //!< \f$ N_{5,\xi} \f$
  Double DxShapeFnc6(Double x, Double y) { return 0.5*(2*x+1)*(1-y*y);}   //!< \f$ N_{6,\xi} \f$
  Double DxShapeFnc7(Double x, Double y) { return 0.5*y*(-2*x)*(y+1);}    //!< \f$ N_{7,\xi} \f$
  Double DxShapeFnc8(Double x, Double y) { return 0.5*(2*x-1)*(1-y*y);}   //!< \f$ N_{8,\xi} \f$
  Double DxShapeFnc9(Double x, Double y) { return (-x*2)*(1-y*y);}        //!< \f$ N_{9,\xi} \f$

  Double DyShapeFnc1(Double x, Double y) { return 0.25*x*(x-1)*(2*y-1);}  //!< \f$ N_{1,\eta} \f$
  Double DyShapeFnc2(Double x, Double y) { return -0.25*x*(x+1)*(2*y-1);} //!< \f$ N_{2,\eta} \f$
  Double DyShapeFnc3(Double x, Double y) { return 0.25*x*(x+1)*(2*y+1);}  //!< \f$ N_{3,\eta} \f$
  Double DyShapeFnc4(Double x, Double y) { return -0.25*x*(x-1)*(2*y+1);} //!< \f$ N_{4,\eta} \f$
  Double DyShapeFnc5(Double x, Double y) { return 0.5*(1-x*x)*(2*y-1);}   //!< \f$ N_{5,\eta} \f$
  Double DyShapeFnc6(Double x, Double y) { return 0.5*x*(x+1)*(-y*2);}    //!< \f$ N_{6,\eta} \f$
  Double DyShapeFnc7(Double x, Double y) { return 0.5*(1-x*x)*(2*y+1);}   //!< \f$ N_{7,\eta} \f$
  Double DyShapeFnc8(Double x, Double y) { return 0.5*x*(x-1)*(-y*2);}    //!< \f$ N_{8,\eta} \f$
  Double DyShapeFnc9(Double x, Double y) { return (1-x*x)*(-y*2);}        //!< \f$ N_{9,\eta} \f$ 

  Vector<Double> ShapeFncAtIP1; //!< contains \f$ N_1 \f$ at all integration points
  Vector<Double> ShapeFncAtIP2; //!< contains \f$ N_2 \f$ at all integration points
  Vector<Double> ShapeFncAtIP3; //!< contains \f$ N_3 \f$ at all integration points
  Vector<Double> ShapeFncAtIP4; //!< contains \f$ N_4 \f$ at all integration points
  Vector<Double> ShapeFncAtIP5; //!< contains \f$ N_5 \f$ at all integration points
  Vector<Double> ShapeFncAtIP6; //!< contains \f$ N_6 \f$ at all integration points
  Vector<Double> ShapeFncAtIP7; //!< contains \f$ N_7 \f$ at all integration points
  Vector<Double> ShapeFncAtIP8; //!< contains \f$ N_8 \f$ at all integration points
  Vector<Double> ShapeFncAtIP9; //!< contains \f$ N_9 \f$ at all integration points

  Vector<Double> DxShapeFncAtIP1; //!< contains \f$ N_{1,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP2; //!< contains \f$ N_{2,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP3; //!< contains \f$ N_{3,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP4; //!< contains \f$ N_{4,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP5; //!< contains \f$ N_{5,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP6; //!< contains \f$ N_{6,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP7; //!< contains \f$ N_{7,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP8; //!< contains \f$ N_{8,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP9; //!< contains \f$ N_{9,\xi} \f$ at all integration points

  Vector<Double> DyShapeFncAtIP1; //!< contains \f$ N_{1,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP2; //!< contains \f$ N_{2,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP3; //!< contains \f$ N_{3,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP4; //!< contains \f$ N_{4,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP5; //!< contains \f$ N_{5,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP6; //!< contains \f$ N_{6,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP7; //!< contains \f$ N_{7,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP8; //!< contains \f$ N_{8,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP9; //!< contains \f$ N_{9,\eta} \f$ at all integration points

};

}

#endif // FILE_DQUADRILATERAL1
