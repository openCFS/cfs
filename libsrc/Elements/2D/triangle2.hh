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
  /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(0,0)\ , \ N_{i,\eta}(0,0))\f$
     \param i (input) shape function number
   */
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

  //! get gradient of shape fnc for node i at center of element
  /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(0,0)\ , \ N_{i,\eta}(0,0))\f$
     \param i (input) shape function number
   */
  virtual  void GetGradientShFncAtCenter(Vector<Double> & ,const Integer i);

  //! Return pointer to array of value shape function at integration points
   /*!
     \param ShFnc (input) shape function number
  */ 
 virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

  //! test type of element
  virtual enum ElementType type(){ return Triang2;}

protected:

private:

  void SetShapeFncAtIntPoints();
  void SetDerShapeFncAtIntPoints();  
  
  Double ShapeFnc1(Double x, Double y) { return 2.*(1.-x-y)*(0.5-x-y); } //!< \f$ N_1 \f$
  Double DxShapeFnc1(Double x, Double y) { return 4.*(x+y)-3.; }         //!< \f$ N_{1,\xi} \f$
  Double DyShapeFnc1(Double x, Double y) { return 4.*(x+y)-3.; }         //!< \f$ N_{1,\eta} \f$

  Double ShapeFnc2(Double x, Double y)   { return x*(2.*x-1.); }         //!< \f$ N_2 \f$
  Double DxShapeFnc2(Double x, Double y) { return 4.*x-1.; }             //!< \f$ N_{2,\xi} \f$
  Double DyShapeFnc2(Double x, Double y) { return 0.; }                  //!< \f$ N_{2,\eta} \f$

  Double ShapeFnc3(Double x, Double y)   { return y*(2.*y-1.); }         //!< \f$ N_3 \f$
  Double DxShapeFnc3(Double x, Double y) { return 0.; }                  //!< \f$ N_{3,\xi} \f$ 
  Double DyShapeFnc3(Double x, Double y) { return 4.*y-1.; }             //!< \f$ N_{3,\eta} \f$

  Double ShapeFnc4(Double x, Double y)   { return 4.*x*y; }              //!< \f$ N_4 \f$
  Double DxShapeFnc4(Double x, Double y) { return 4.*y; }                //!< \f$ N_{4,\xi} \f$ 
  Double DyShapeFnc4(Double x, Double y) { return 4.*x; }                //!< \f$ N_{4,\eta} \f$ 

  Double ShapeFnc5(Double x, Double y)   { return 4.*y*(1.-x-y); }       //!< \f$ N_5 \f$
  Double DxShapeFnc5(Double x, Double y) { return -4.*y; }               //!< \f$ N_{5,\xi} \f$ 
  Double DyShapeFnc5(Double x, Double y) { return 4.*(1-x-2.*y); }       //!< \f$ N_{5,\eta} \f$ 

  Double ShapeFnc6(Double x, Double y)   { return 4.*x*(1.-x-y); }       //!< \f$ N_6 \f$
  Double DxShapeFnc6(Double x, Double y) { return 4.*(1-2.*x-y); }       //!< \f$ N_{6,\xi} \f$ 
  Double DyShapeFnc6(Double x, Double y) { return -4.*x; }               //!< \f$ N_{6,\xi} \f$ 
  
  Vector<Double> ShapeFncAtIP1; //!< contains \f$ N_1 \f$ at all integration points
  Vector<Double> ShapeFncAtIP2; //!< contains \f$ N_2 \f$ at all integration points
  Vector<Double> ShapeFncAtIP3; //!< contains \f$ N_3 \f$ at all integration points
  Vector<Double> ShapeFncAtIP4; //!< contains \f$ N_4 \f$ at all integration points
  Vector<Double> ShapeFncAtIP5; //!< contains \f$ N_5 \f$ at all integration points
  Vector<Double> ShapeFncAtIP6; //!< contains \f$ N_6 \f$ at all integration points

  Vector<Double> DxShapeFncAtIP1; //!< contains \f$ N_{1,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP2; //!< contains \f$ N_{2,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP3; //!< contains \f$ N_{3,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP4; //!< contains \f$ N_{4,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP5; //!< contains \f$ N_{5,\xi} \f$ at all integration points
  Vector<Double> DxShapeFncAtIP6; //!< contains \f$ N_{6,\xi} \f$ at all integration points

  Vector<Double> DyShapeFncAtIP1; //!< contains \f$ N_{1,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP2; //!< contains \f$ N_{2,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP3; //!< contains \f$ N_{3,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP4; //!< contains \f$ N_{4,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP5; //!< contains \f$ N_{5,\eta} \f$ at all integration points
  Vector<Double> DyShapeFncAtIP6; //!< contains \f$ N_{6,\eta} \f$ at all integration points
 
};
 
}
 
#endif // FILE_TRIANGLE2_2003
