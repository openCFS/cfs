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
  /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(ip)\ , \ N_{i,\eta}(ip))\f$
     \param i (input) shape function number
     \param ip (input) integration point
  */
  virtual void GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip);  

  //! get gradient of shape fnc for node i at center of element
  /*!
     \param grad (output) contains local derivatives \f$ (N_{i,\xi}(0,0)\ , \ N_{i,\eta}(0,0))\f$
     \param ip (input) integration point
  */
  virtual  void GetGradientShFncAtCenter(Vector<Double> & ,const Integer i);

  //! Return pointer to array of value shape function at integration points
  /*!
     \param ShFnc (input) shape function number
  */
  virtual Vector<Double> & GetShFncAtIP(const Integer ShFnc) ;

  //! test type of element
  virtual enum ElementType type(){ return Quadrilateral1;}

protected:

private:
};

}

#endif // FILE_DQUADRILATERAL1
