#ifndef FILE_PDE_2001
#define FILE_PDE_2001

#include "material.hh"
#include "filetype.hh"
 
namespace CoupledField
{

  //! Base class for partial differential equation
  /*! Class PDE is base class, from which different type of PDE are derived. This class contains general things for different types of equations, for example, data about intergration parameters, method to calculate parameters for Newmark method.
   */

class PDE
{
public:

  //! Constructor( read integration parameters, define class Material)
  PDE(FileType * , Material *);

  //! Calculation parameteres for Newmark method for the time step
  void CalcParameters(const Double dt);

  //! Deconstructor
  virtual ~PDE(){;}

  //! Return solution
  virtual Vector<Double> & getS() { Error("Not implemented");}

 //! Return first derivatives of solution
  virtual Vector<Double> & getS1() { Error("Not implemented"); }

 //! Return second derivatives of solution
  virtual Vector<Double> & getS2() { Error("Not implemented"); }
 
protected:
  
  //! pointer to class Material
  Material * ptMaterial;

  //! Intergration parameters
  Double alpha,gamma, beta;

  //!
  Double a0,a1,a2,a3,a4,a5,a6,a7;
};

} // end of namespace

#endif
