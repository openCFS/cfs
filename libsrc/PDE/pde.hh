#ifndef FILE_PDE_2001
#define FILE_PDE_2001
 
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
  void CalcParamForNewmarkMethod(const Double dt);

  //! Deconstructor
  virtual ~PDE(){;}
 
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
