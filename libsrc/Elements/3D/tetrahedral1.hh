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
   /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(ip)\ , \ N_{i,\eta}(ip)\ , \ N_{i,zeta})\f$
     \param i (input) shape function number
     \param ip (input) integration point
   */
   virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);

   //! Return pointer to array of value shape function at integration points
   /*!
     \param iShFnc (input) shape function number
   */ 
   virtual Vector<Double> & GetShFncAtIP(const Integer iShFnc);

   //! return values of gradient at center point
   /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(0,0,0)\ , \ N_{i,\eta}(0,0,0))\f$
     \param iShFnc (input) shape function number
   */
   void GetGradientShFncAtCenter(Vector<Double>&grad,const Integer iShFnc);

private:
   
};

}
#endif //
