#ifndef FILE_DQUADRILATERAL1
#define FILE_DQUADRILATERAL1

#include "rectangle.hh"

namespace CoupledField
{
//! Quadrilateral finite element with four nodes (linear interpolation function)

  class Quad1 : public Rectangle
{
public:

  //! Constructor with type of integration rule
  Quad1();
  
  //! Deconstructor
  virtual ~Quad1();

  //! Define variables of this class
  virtual void Init();

  //! Return Vector<Double> of gradient of shape func with number i at integration point ip 
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);  
 
  //! Return pointer to array of value shape function at integration points
  virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

protected:

private:
  ShortInt ElemType;  //!< element type
};

}

#endif // FILE_DQUADRILATERAL1
