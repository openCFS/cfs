#ifndef FILE_GEOMETRY_LINE_2002
#define FILE_GEOMETRY_LINE_2002

#include <Elements/baseelem.hh>
#include <Elements/jacobian.hh>

namespace CoupledField
{
//! Class with general procedures for finite elements: line
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights
*/

class Line : public BaseElem
{
public:

  //! Constructor with type of integration rule
  Line();

  //! Deconstructor
  virtual ~Line();

  //! return FE-Type for CLA++
  virtual Integer feType() { return 1;}

    //!
  virtual void Init();


   //!
  virtual  void GetGradientShFnc(Vector<Double> & ,const Integer i, const Integer j)
    { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

 //!
  virtual Vector<Double> &  GetShFncAtIP(const Integer iShFnc);
//  { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

  //! Define integration points for this type of integration
  virtual void SetIntPoints();

  //! Calculate transformation function at these intergration points
  //virtual void SetTransformFncAtIntPoints()  { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}
	virtual void SetTransformFncAtIntPoints();

     //! Calculate derivative of transformation function at these intergration points
  //virtual void SetDerTransformFncAtIntPoints()  { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}
	virtual void SetDerTransformFncAtIntPoints();

   //! Calculation of Jacobian in 2D
  virtual void CalcJacobian(Jacobian<2> & J, const Integer ip,
          Point<2> *  ptCoord, const Boolean NeedJinv=TRUE)
 { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

  //! Calculation of Jacobian in 3D
 virtual void CalcJacobian(Jacobian<3> & J, const Integer ip,
              Point<3> *  ptCoord, const Boolean NeedJinv=TRUE)
 { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}


  virtual Double getIntVal(const Point<2> * const ptCoord);

protected:

  Vector<Double> TransFnc1AtIP;
  Vector<Double> TransFnc2AtIP;

  Vector<Double> DxiTransFnc1AtIP;
  Vector<Double> DxiTransFnc2AtIP;


	Boolean IsSet; //!< parameter is TRUE if we calculated transformation function at integration points
private:

	Double TransFnc1(Double xi) {return 0.5*(1-xi); }
  Double TransFnc2(Double xi) {return 0.5*(1+xi); }

  Double TransFnc1dxi () { return -0.5;}
  Double TransFnc2dxi () { return  0.5;}

};

}

#endif // FILE_LINE
