#ifndef FILE_RECTANGLE_2002
#define FILE_RECTANGLE_2002

#include <Elements/baseelem.hh>
#include <Elements/jacobian.hh>

namespace CoupledField
{
//! Class with general procedures for quadrilateral finite elements
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights */
  
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

  //! Calculate shape functions at intergration points
  virtual void SetTransformFncAtIntPoints();

     //! Calculate derivatives of shape functions at these intergration points
  virtual void SetDerTransformFncAtIntPoints();

     //! Calculate derivatives of shape functions at these intergration points
  virtual void SetDerTransformFncAtIntPoints3D();
 
  //! Calculate shape functions at the center of element
  virtual void SetTransformFncAtCenter();

     //! Calculate derivatives of shape functions at the center of element
  virtual void SetDerTransformFncAtCenter();
  
  //! Calculation of Jacobian, inverse Jacobian and detJacobian
  /*! 
      \param J (output) Jacobian matrix (2,2)
      \param ip (input) number of integration point
      \param ptCoord (input) contains global coordinates of element nodes
      \param NeedJinv (input) TRUE: inverse Jacobian computed; FALSE: no inverse Jacobian computed  
  */
  virtual void CalcJacobian(Jacobian<2> & J, const Integer ip,
                 Point<2> * ptCoord, const Boolean NeedJinv=TRUE);

  // Calculation of Jacobian, inverse Jacobian and detJacobian
  virtual void CalcJacobian(Jacobian<3> & J, const Integer ip,
                  Point<3> * ptCoord, const Boolean NeedJinv=TRUE);

  //! Calculation of Jacobian for center point in 2D
  /*! 
      \param J (output) Jacobian matrix (2,2)
       \param ptCoord (input) contains global coordinates of element nodes
      \param NeedJinv (input) TRUE: inverse Jacobian computed; FALSE: no inverse Jacobian computed
  */
  virtual void CalcJacobianAtCenter(Jacobian<2> & J, Point<2> * ptCoord, const Boolean NeedJinv=TRUE);
  
  Vector<Double> TransFncAtIP1; //!< contains \f$ N_1 \f$ at all integration points
  Vector<Double> TransFncAtIP2; //!< contains \f$ N_2 \f$ at all integration points
  Vector<Double> TransFncAtIP3; //!< contains \f$ N_3 \f$ at all integration points
  Vector<Double> TransFncAtIP4; //!< contains \f$ N_4 \f$ at all integration points
 
  Vector<Double> DxTransFncAtIP1; //!< contains \f$ N_{1,\xi} \f$ at all integration points
  Vector<Double> DxTransFncAtIP2; //!< contains \f$ N_{2,\xi} \f$ at all integration points
  Vector<Double> DxTransFncAtIP3; //!< contains \f$ N_{3,\xi} \f$ at all integration points
  Vector<Double> DxTransFncAtIP4; //!< contains \f$ N_{4,\xi} \f$ at all integration points
 
  Vector<Double> DyTransFncAtIP1; //!< contains \f$ N_{1,\eta} \f$ at all integration points
  Vector<Double> DyTransFncAtIP2; //!< contains \f$ N_{2,\eta} \f$ at all integration points
  Vector<Double> DyTransFncAtIP3; //!< contains \f$ N_{3,\eta} \f$ at all integration points
  Vector<Double> DyTransFncAtIP4; //!< contains \f$ N_{4,\eta} \f$ at all integration points

  // These derivatives are used for the computation of Jacobian 3d. 
  Vector<Double> DxiTransFnc1AtIP; //!< contains 3D \f$ N_{1,\xi} \f$ at all integration points
  Vector<Double> DxiTransFnc2AtIP; //!< contains 3D \f$ N_{2,\xi} \f$ at all integration points
  Vector<Double> DxiTransFnc3AtIP; //!< contains 3D \f$ N_{3,\xi} \f$ at all integration points
  Vector<Double> DxiTransFnc4AtIP; //!< contains 3D \f$ N_{4,\xi} \f$ at all integration points
 
  Vector<Double> DetaTransFnc1AtIP; //!< contains 3D \f$ N_{1,\eta} \f$ at all integration points
  Vector<Double> DetaTransFnc2AtIP; //!< contains 3D \f$ N_{2,\eta} \f$ at all integration points
  Vector<Double> DetaTransFnc3AtIP; //!< contains 3D \f$ N_{3,\eta} \f$ at all integration points
  Vector<Double> DetaTransFnc4AtIP; //!< contains 3D \f$ N_{4,\eta} \f$ at all integration points

  Double TransFncAtCenter[4]; //!< contains \f$ N_1\f$ - \f$N_4\f$ at center 
  Double DxTransFncAtCenter[4]; //!< contains \f$ N_{1,\xi}\f$ - \f$N_{4,\xi}\f$ at center
  Double DyTransFncAtCenter[4]; //!< contains \f$ N_{1,\eta}\f$ - \f$N_{4,\eta}\f$ at center

  Boolean IsSet; //!< parameter is TRUE if we have already computed transformation function at integration points
  Boolean isSetAtCenter_; //!< parameter is TRUE if we calculated tranformation function at integration points

private:
  Double TransFnc1(Double x,Double y) {return 0.25*(1+x)*(1+y); } //!< \f$ N_1 \f$
  Double TransFnc2(Double x,Double y) {return 0.25*(1-x)*(1+y); } //!< \f$ N_2 \f$
  Double TransFnc3(Double x,Double y) {return 0.25*(1-x)*(1-y); } //!< \f$ N_3 \f$
  Double TransFnc4(Double x,Double y) {return 0.25*(1+x)*(1-y); } //!< \f$ N_4 \f$

  Double TransFnc1dx (Double x,Double y)  { return 0.25*(1+y);}  //!< \f$ N_{1,\xi} \f$
  Double TransFnc1dy (Double x,Double y)  { return 0.25*(1+x);}  //!< \f$ N_{1,\eta} \f$
  Double TransFnc2dx (Double x,Double y)  { return 0.25*(-1-y);} //!< \f$ N_{2,\xi} \f$
  Double TransFnc2dy (Double x,Double y)  { return 0.25*(1-x);}  //!< \f$ N_{2,\eta} \f$
  Double TransFnc3dx (Double x,Double y)  { return 0.25*(y-1);}  //!< \f$ N_{3,\xi} \f$
  Double TransFnc3dy (Double x,Double y)  { return 0.25*(x-1);}  //!< \f$ N_{3,\eta} \f$
  Double TransFnc4dx (Double x,Double y)  { return 0.25*(1-y);}  //!< \f$ N_{4,\xi} \f$
  Double TransFnc4dy (Double x,Double y)  { return 0.25*(-1-x);} //!< \f$ N_{4,\eta} \f$



  Double TransFnc1dxi (Double xi,Double et, Double ze) { return -0.125*(1-et)*(1-ze); }   //!< \f$ N_{1,\xi} \f$
  Double TransFnc1deta (Double xi,Double et, Double ze) { return -0.125*(1-xi)*(1-ze); }  //!< \f$ N_{1,\eta} \f$

  Double TransFnc2dxi (Double xi,Double et, Double ze) { return 0.125*(1-et)*(1-ze); }    //!< \f$ N_{2,\xi} \f$
  Double TransFnc2deta (Double xi,Double et, Double ze) { return -0.125*(1+xi)*(1-ze); }  //!< \f$ N_{2,\eta} \f$

  Double TransFnc3dxi (Double xi,Double et, Double ze)  { return 0.125*(1+et)*(1-ze); }    //!< \f$ N_{3,\xi} \f$
  Double TransFnc3deta (Double xi,Double et, Double ze)  { return 0.125*(1+xi)*(1-ze); }   //!< \f$ N_{3,\eta} \f$

  Double TransFnc4dxi (Double xi,Double et, Double ze)  { return -0.125*(1+et)*(1-ze); }   //!< \f$ N_{4,\xi} \f$
  Double TransFnc4deta (Double xi,Double et, Double ze)  { return 0.125*(1-xi)*(1-ze); }   //!< \f$ N_{4,\eta} \f$

};

}

#endif // FILE_RECTANGLE
