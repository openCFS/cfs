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

  //! return FE-Type for LAS++
  virtual Integer feType() { return TET;}
  
///////////////////////////////////////////////////////////////////////

virtual Vector<Double> *  GetDxShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDyShFncAtIP(const Integer iShFnc) ;
virtual Vector<Double> *  GetDzShFncAtIP(const Integer iShFnc) ;

protected:

   void test(){ std::cout << " tetrahedral " << std::endl;}

   //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate shape functions at these integration points
  virtual void SetTransformFncAtIntPoints();

  //! Calculate derivatives of shape functions at these integration points
  virtual void SetDerTransformFncAtIntPoints();
  
  //! Calculate shape functions at the center of element
  virtual void SetTransformFncAtCenter();

  //! Calculate derivatives of shape functions at the center of element
  virtual void SetDerTransformFncAtCenter();

  // Calculation of Jacobian, inverse Jacobian and detJacobian
  void CalcJacobian(Jacobian<2> & J, const Integer ip,
                 Point<2> * ptCoord, const Boolean NeedJinv=TRUE); 
    
  //! Calculation of Jacobian, inverse Jacobian and detJacobian
  /*! 
      \param J (output) Jacobian matrix (3,3)
      \param ip (input) number of integration point
      \param ptCoord (input) contains global coordinates of element nodes
      \param NeedJinv (input) TRUE: inverse Jacobian computed; FALSE: no inverse Jacobian computed  
  */
  void CalcJacobian(Jacobian<3> & J, const Integer ip,
                   Point<3> * ptCoord, const Boolean NeedJinv=TRUE);

 //! Calculation of Jacobian for center point in 3D
  /*! 
      \param J (output) Jacobian matrix (3,3)
        \param ptCoord (input) contains global coordinates of element nodes
      \param NeedJinv (input) TRUE: inverse Jacobian computed; FALSE: no inverse Jacobian computed  
  */
  virtual void CalcJacobianAtCenter(Jacobian<3> & J,  Point<3> *  ptCoord, const Boolean NeedJinv=TRUE);

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

  Vector<Double> DzTransFncAtIP1; //!< contains \f$ N_{1,\zeta} \f$ at all integration points
  Vector<Double> DzTransFncAtIP2; //!< contains \f$ N_{2,\zeta} \f$ at all integration points
  Vector<Double> DzTransFncAtIP3; //!< contains \f$ N_{3,\zeta} \f$ at all integration points
  Vector<Double> DzTransFncAtIP4; //!< contains \f$ N_{4,\zeta} \f$ at all integration points

  Double TransFncAtCenter[4];   //!< contains \f$ N_1\f$ - \f$N_4\f$ at center 
  Double DxTransFncAtCenter[4]; //!< contains \f$ N_{1,\xi}\f$ - \f$N_{4,\xi}\f$ at center
  Double DyTransFncAtCenter[4]; //!< contains \f$ N_{1,\eta}\f$ - \f$N_{4,\eta}\f$ at center
  Double DzTransFncAtCenter[4]; //!< contains \f$ N_{1,\zeta}\f$ - \f$N_{4,\zeta}\f$ at center

  Boolean IsSet; //!< parameter is TRUE if we have already computed transformation function at integration points
  Boolean isSetAtCenter_; //!< parameter is TRUE if we calculated tranformation function at integration points

private:

  Double TransFnc1(Double x,Double y, Double z) { return 1-x-y-z; } //!< \f$ N_1 \f$
  Double TransFnc2(Double x,Double y, Double z) { return x; }       //!< \f$ N_2 \f$
  Double TransFnc3(Double x,Double y, Double z) { return y; }       //!< \f$ N_3 \f$
  Double TransFnc4(Double x,Double y, Double z) { return z; }       //!< \f$ N_4 \f$

  Double TransFnc1dx (Double x,Double y, Double z) { return -1.0; }  //!< \f$ N_{1,\xi} \f$
  Double TransFnc1dy (Double x,Double y, Double z) { return -1.0; }  //!< \f$ N_{1,\eta} \f$
  Double TransFnc1dz (Double x,Double y, Double z) { return -1.0; }  //!< \f$ N_{1,\zeta} \f$

  Double TransFnc2dx (Double x,Double y, Double z) { return 1.0; }   //!< \f$ N_{2,\xi} \f$
  Double TransFnc2dy (Double x,Double y, Double z) { return 0.0; }   //!< \f$ N_{2,\eta} \f$
  Double TransFnc2dz (Double x,Double y, Double z) { return 0.0; }   //!< \f$ N_{2,\zeta} \f$

  Double TransFnc3dx (Double x,Double y, Double z)  { return 0.0; }  //!< \f$ N_{3,\xi} \f$
  Double TransFnc3dy (Double x,Double y, Double z)  { return 1.0; }  //!< \f$ N_{3,\eta} \f$
  Double TransFnc3dz (Double x,Double y, Double z)  { return 0.0; }  //!< \f$ N_{3,\zeta} \f$

  Double TransFnc4dx (Double x,Double y, Double z)  { return 0.0; }  //!< \f$ N_{4,\xi} \f$ 
  Double TransFnc4dy (Double x,Double y, Double z)  { return 0.0; }  //!< \f$ N_{4,\eta} \f$ 
  Double TransFnc4dz (Double x,Double y, Double z)  { return 1.0; }  //!< \f$ N_{4,\zeta} \f$
};

}
#endif //
