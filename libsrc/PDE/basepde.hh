#ifndef FILE_BASEPDE_2001
#define FILE_BASEPDE_2001

#include "material.hh"
#include "filetype.hh"
#include "abstractAlgSys.hh"
#include "element_header.hh"

namespace CoupledField
{

  //! Base class for partial differential equation
  /*! Class BasePDE is base class, from which different type of PDEs are derived. 
   */

class BasePDE
{
public:

  //! Constructor( read integration parameters, define class Material)
  BasePDE(FileType * , Material *);

  //! Deconstructor
  virtual ~BasePDE(){;}

  //!
  virtual void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps, Double &adampiter, 
                     Integer &amaxnumit)=0;
  //!
  virtual void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, 
                               Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)=0;
  //!
  virtual void CalcElemStiff()=0;

  //!
  virtual void AssembleGlobal(AbstractAlgebraicSys * algsys)=0;

  //! Calculation of integration parameteres 
  void CalcIntegrationParam(const Double dt);
 
protected:
  
  //! pointer to class Material
  Material * ptMaterial;

  //! Intergration parameters
  Double alpha,gamma_hyperbolic, gamma_parabolic, beta;

  //!
  Double cinteg[13];

};

} // end of namespace

#endif
