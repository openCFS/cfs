#ifndef FILE_GEHEXAHEDRAL_2003
#define FILE_GEHEXAHEDRAL_2003

#include <Elements/baseelem.hh>

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

  //! return FE-Type for LAS++
  virtual Integer feType() { return HEX;}
	
///////////////////////////////////////////////////////////////////////

virtual Vector<Double> *  GetDXiShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDEtShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDZeShFncAtIP(const Integer iShFnc) ;

protected:

   void test(){ std::cout << " hexahedral " << std::endl;}

   //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate shape functions at these integration points
  virtual void SetTransformFncAtIntPoints();

  //! Calculate derivatives of shape functions at these integration points
  virtual void SetDerTransformFncAtIntPoints();
  
  // Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<2> & J, const Integer ip,
                 Point<2> *  ptCoord, const Boolean NeedJinv=TRUE); 
    
  //! Calculation of Jacobian, inverse Jacobian and detJacobian
  /*! 
      \param J (output) Jacobian matrix (3,3)
      \param ip (input) number of integration point
      \param ptCoord (input) contains global coordinates of element nodes
      \param NeedJinv (input) TRUE: inverse Jacobian computed; FALSE: no inverse Jacobian computed  
  */
  void CalcJacobian(Jacobian<3> & J, const Integer ip,
                  Point<3> *  ptCoord, const Boolean NeedJinv=TRUE);

  //! Calculation of Jacobian, inverse Jacobian and detJacobian at center of element
  /*! 
      \param J (output) Jacobian matrix (3,3)
        \param ptCoord (input) contains global coordinates of element nodes
      \param NeedJinv (input) TRUE: inverse Jacobian computed; FALSE: no inverse Jacobian computed  
  */
  void CalcJacobianAtCenterOfElem(Jacobian<3> & J,
                 Point<3> *  ptCoord, const Boolean NeedJinv=TRUE);

  Vector<Double> TransFnc1AtIP; //!< contains \f$ N_1 \f$ at all integration points
  Vector<Double> TransFnc2AtIP; //!< contains \f$ N_2 \f$ at all integration points
  Vector<Double> TransFnc3AtIP; //!< contains \f$ N_3 \f$ at all integration points
  Vector<Double> TransFnc4AtIP; //!< contains \f$ N_4 \f$ at all integration points
  Vector<Double> TransFnc5AtIP; //!< contains \f$ N_5 \f$ at all integration points
  Vector<Double> TransFnc6AtIP; //!< contains \f$ N_6 \f$ at all integration points
  Vector<Double> TransFnc7AtIP; //!< contains \f$ N_7 \f$ at all integration points
  Vector<Double> TransFnc8AtIP; //!< contains \f$ N_8 \f$ at all integration points

  Vector<Double> DXiTransFnc1AtIP; //!< contains \f$ N_{1,\xi} \f$ at all integration points
  Vector<Double> DXiTransFnc2AtIP; //!< contains \f$ N_{2,\xi} \f$ at all integration points
  Vector<Double> DXiTransFnc3AtIP; //!< contains \f$ N_{3,\xi} \f$ at all integration points
  Vector<Double> DXiTransFnc4AtIP; //!< contains \f$ N_{4,\xi} \f$ at all integration points
  Vector<Double> DXiTransFnc5AtIP; //!< contains \f$ N_{5,\xi} \f$ at all integration points
  Vector<Double> DXiTransFnc6AtIP; //!< contains \f$ N_{6,\xi} \f$ at all integration points
  Vector<Double> DXiTransFnc7AtIP; //!< contains \f$ N_{7,\xi} \f$ at all integration points
  Vector<Double> DXiTransFnc8AtIP; //!< contains \f$ N_{8,\xi} \f$ at all integration points

  Vector<Double> DEtTransFnc1AtIP; //!< contains \f$ N_{1,\eta} \f$ at all integration points
  Vector<Double> DEtTransFnc2AtIP; //!< contains \f$ N_{2,\eta} \f$ at all integration points
  Vector<Double> DEtTransFnc3AtIP; //!< contains \f$ N_{3,\eta} \f$ at all integration points
  Vector<Double> DEtTransFnc4AtIP; //!< contains \f$ N_{4,\eta} \f$ at all integration points
  Vector<Double> DEtTransFnc5AtIP; //!< contains \f$ N_{5,\eta} \f$ at all integration points
  Vector<Double> DEtTransFnc6AtIP; //!< contains \f$ N_{6,\eta} \f$ at all integration points
  Vector<Double> DEtTransFnc7AtIP; //!< contains \f$ N_{7,\eta} \f$ at all integration points
  Vector<Double> DEtTransFnc8AtIP; //!< contains \f$ N_{8,\eta} \f$ at all integration points

  Vector<Double> DZeTransFnc1AtIP; //!< contains \f$ N_{1,\zeta} \f$ at all integration points
  Vector<Double> DZeTransFnc2AtIP; //!< contains \f$ N_{2,\zeta} \f$ at all integration points
  Vector<Double> DZeTransFnc3AtIP; //!< contains \f$ N_{3,\zeta} \f$ at all integration points
  Vector<Double> DZeTransFnc4AtIP; //!< contains \f$ N_{4,\zeta} \f$ at all integration points
  Vector<Double> DZeTransFnc5AtIP; //!< contains \f$ N_{5,\zeta} \f$ at all integration points
  Vector<Double> DZeTransFnc6AtIP; //!< contains \f$ N_{6,\zeta} \f$ at all integration points
  Vector<Double> DZeTransFnc7AtIP; //!< contains \f$ N_{7,\zeta} \f$ at all integration points
  Vector<Double> DZeTransFnc8AtIP; //!< contains \f$ N_{8,\zeta} \f$ at all integration points

  Boolean IsSet;

private:

  Double TransFnc1(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1-et)*(1-ze)); } //!< \f$ N_1 \f$
  Double TransFnc2(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1-et)*(1-ze)); } //!< \f$ N_2 \f$
  Double TransFnc3(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1+et)*(1-ze)); } //!< \f$ N_3 \f$
  Double TransFnc4(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1+et)*(1-ze)); } //!< \f$ N_4 \f$
  Double TransFnc5(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1-et)*(1+ze)); } //!< \f$ N_5 \f$
  Double TransFnc6(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1-et)*(1+ze)); } //!< \f$ N_6 \f$ 
  Double TransFnc7(Double xi,Double et, Double ze) { return 0.125*((1+xi)*(1+et)*(1+ze)); } //!< \f$ N_7 \f$ 
  Double TransFnc8(Double xi,Double et, Double ze) { return 0.125*((1-xi)*(1+et)*(1+ze)); } //!< \f$ N_8 \f$ 

  Double TransFnc1dxi (Double xi,Double et, Double ze) { return -0.125*(1-et)*(1-ze); }   //!< \f$ N_{1,\xi} \f$
  Double TransFnc1deta (Double xi,Double et, Double ze) { return -0.125*(1-xi)*(1-ze); }  //!< \f$ N_{1,\eta} \f$
  Double TransFnc1dzeta (Double xi,Double et, Double ze) { return -0.125*(1-xi)*(1-et); } //!< \f$ N_{1,\zeta} \f$

  Double TransFnc2dxi (Double xi,Double et, Double ze) { return 0.125*(1-et)*(1-ze); }    //!< \f$ N_{2,\xi} \f$
  Double TransFnc2deta (Double xi,Double et, Double ze) { return -0.125*(1+xi)*(1-ze); }  //!< \f$ N_{2,\eta} \f$
  Double TransFnc2dzeta (Double xi,Double et, Double ze) { return -0.125*(1+xi)*(1-et); } //!< \f$ N_{2,\zeta} \f$

  Double TransFnc3dxi (Double xi,Double et, Double ze)  { return 0.125*(1+et)*(1-ze); }    //!< \f$ N_{3,\xi} \f$
  Double TransFnc3deta (Double xi,Double et, Double ze)  { return 0.125*(1+xi)*(1-ze); }   //!< \f$ N_{3,\eta} \f$
  Double TransFnc3dzeta (Double xi,Double et, Double ze)  { return -0.125*(1+xi)*(1+et); } //!< \f$ N_{3,\zeta} \f$ 

  Double TransFnc4dxi (Double xi,Double et, Double ze)  { return -0.125*(1+et)*(1-ze); }   //!< \f$ N_{4,\xi} \f$
  Double TransFnc4deta (Double xi,Double et, Double ze)  { return 0.125*(1-xi)*(1-ze); }   //!< \f$ N_{4,\eta} \f$
  Double TransFnc4dzeta (Double xi,Double et, Double ze)  { return -0.125*(1-xi)*(1+et); } //!< \f$ N_{4,\zeta} \f$

  Double TransFnc5dxi (Double xi,Double et, Double ze) { return -0.125*(1-et)*(1+ze); }    //!< \f$ N_{5,\xi} \f$
  Double TransFnc5deta (Double xi,Double et, Double ze) { return -0.125*(1-xi)*(1+ze); }   //!< \f$ N_{5,\eta} \f$
  Double TransFnc5dzeta (Double xi,Double et, Double ze) { return 0.125*(1-xi)*(1-et); }   //!< \f$ N_{5,\zeta} \f$

  Double TransFnc6dxi (Double xi,Double et, Double ze) { return 0.125*(1-et)*(1+ze); }     //!< \f$ N_{6,\xi} \f$
  Double TransFnc6deta (Double xi,Double et, Double ze) { return -0.125*(1+xi)*(1+ze); }   //!< \f$ N_{6,\eta} \f$
  Double TransFnc6dzeta (Double xi,Double et, Double ze) { return 0.125*(1+xi)*(1-et); }   //!< \f$ N_{6,\zeta} \f$

  Double TransFnc7dxi (Double xi,Double et, Double ze)  { return 0.125*(1+et)*(1+ze); }    //!< \f$ N_{7,\xi} \f$
  Double TransFnc7deta (Double xi,Double et, Double ze)  { return 0.125*(1+xi)*(1+ze); }   //!< \f$ N_{7,\eta} \f$
  Double TransFnc7dzeta (Double xi,Double et, Double ze)  { return 0.125*(1+xi)*(1+et); }  //!< \f$ N_{7,\zeta} \f$
 
  Double TransFnc8dxi (Double xi,Double et, Double ze)  { return -0.125*(1+et)*(1+ze); }   //!< \f$ N_{8,\xi} \f$
  Double TransFnc8deta (Double xi,Double et, Double ze)  { return 0.125*(1-xi)*(1+ze); }   //!< \f$ N_{8,\eta} \f$
  Double TransFnc8dzeta (Double xi,Double et, Double ze)  { return 0.125*(1-xi)*(1+et); }  //!< \f$ N_{8,\zeta} \f$
};

}
#endif //
