#ifndef FILE_TETRAHEDRAL_2002
#define FILE_TETRAHEDRAL_2002

#include "getetrahedral.hh"

namespace CoupledField
{

//! Class with  description of tetrahedral element

class Tetrahedral1 : public GeTetrahedral
{
public:
   //! Constructor with type of integration rule
   Tetrahedral1();

   //! Deconstructor
   virtual ~Tetrahedral1();

   //! Define variables of this class
   virtual void Init();

   //! Store gradient of shape function with number i at integration point ip
   virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

   //! Return pointer to array of value shape function at integration points
   virtual Vector<Double> & GetShFncAtIP(const Integer iShFnc);

private:
   ShortInt ElemType; //!< element type
   
};

}
#endif //
