#ifndef FILE_TRIANGLE1_2002
#define FILE_TRIANGLE1_2002

#include "getriangle.hh"
 
namespace CoupledField
{

//! Description of triangle element with 1 order interpolation 

  class Triangle1 : public GeTriangle
{
public:

  //! Constructor with type of integration rule
  Triangle1();
 
  //! Deconstructor
  virtual ~Triangle1();
 
  //! Define variables of this class
  virtual void Init();

  //! Return Vector<Double> of gradient of shape func with number i at integration point ip
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

  //! get gradient of shape fnc for node i at center of element
  virtual  void GetGradientShFncAtCenter(Vector<Double> & ,const Integer i);

  //! Return pointer to array of value shape function at integration points
  virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

  //! test type of element
  virtual enum ElementType type(){ return Triang1;}

protected:

private:

};
 
}
 
#endif // FILE_TRIANGLE1_2001
