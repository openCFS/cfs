#ifndef FILE_GEOMETRY_LINE_2002
#define FILE_GEOMETRY_LINE_2002

#include "baseelem.hh"

namespace CoupledField
{
//! Class with general procedures for finite elements: line
/*! This class is derived from BaseElem. It stores general procedures for each type of finite element on quadralateral, such as calculation of Jacobian of transformation in standart and information about integration points and integration weights 
*/
  
class Line : public BaseElem
{
public:

  //! Constructor with type of integration rule
  Line(){;} 
  
  //! Deconstructor
  virtual ~Line(){;}

  //! return FE-Type for CLA++
  virtual Integer feType() { return NOFETYPE;}

    //! 
  virtual void Init()
   { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

   //!
  virtual  void GetGradientShFnc(Vector<Double> & ,const Integer i, const Integer j)
    { Error("This func is not implemented in 1d class",__FILE__,__LINE__);} 

 //!
  virtual Vector<Double> &  GetShFncAtIP(const Integer iShFnc)
  { Error("This func is not implemented in 1d class",__FILE__,__LINE__);} 

protected:

  //! Define integration points for this type of integration
  virtual void SetIntPoints() { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

  //! Calculate transformation function at these intergration points
  virtual void SetTransformFncAtIntPoints()  { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

     //! Calculate derivative of transformation function at these intergration points
  virtual void SetDerTransformFncAtIntPoints()  { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}
 
  // Calculation of Jacobian, inverse Jacobian and detJacobian
  virtual void CalcJacobian(Jacobian<2> & J, const Integer ip,
                 Point<2> * ptCoord, const Boolean NeedJinv=TRUE)
 { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

    // Calculation of Jacobian, inverse Jacobian and detJacobian
  virtual void CalcJacobian(Jacobian<3> & J, const Integer ip,
                  Point<3> * ptCoord, const Boolean NeedJinv=TRUE)
 { Error("This func is not implemented in 1d class",__FILE__,__LINE__);}

  virtual Double getIntVal(Point<2> * ptCoord);

private:

};

}

#endif // FILE_LINE
