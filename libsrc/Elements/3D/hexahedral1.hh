#ifndef FILE_HEXAHEDRAL_2003
#define FILE_HEXAHEDRAL_2003

#include "gehexahedral.hh"

namespace CoupledField
{

//! Class with  description of hexahedral element

class Hexahedral1 : public GeHexahedral
{
public:
   //! Constructor with type of integration rule
   Hexahedral1();

   //! Deconstructor
   virtual ~Hexahedral1();

   //! Define variables of this class
   virtual void Init();

   //! Store gradient of shape function with number i at integration point ip
   virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

   //! Store gradient of shape function with number i at center of element (0,0,0)
   virtual void GetGradientShFncAtCenterOfElem(Vector<Double> & grad, const Integer i);

   //! Return pointer to array of value shape function at integration points
   virtual Vector<Double> & GetShFncAtIP(const Integer iShFnc);

private:
   
};

}
#endif //
