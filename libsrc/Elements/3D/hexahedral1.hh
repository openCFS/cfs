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
   /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(ip)\ , \ N_{i,\eta}(ip)\ , \ N_{i,zeta})\f$
     \param i (input) shape function number
     \param ip (input) integration point
   */
   virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

   //! Store gradient of shape function with number i at center of element (0,0,0)
   /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(0,0,0)\ , \ N_{i,\eta}(0,0,0))\f$
     \param i (input) shape function number
   */
   virtual void GetGradientShFncAtCenterOfElem(Vector<Double> & grad, const Integer i);

   //! Return pointer to array of value shape function at integration points
   /*!
     \param iShFnc (input) shape function number
   */ 
   virtual Vector<Double> & GetShFncAtIP(const Integer iShFnc);

private:
   
};

}
#endif //
