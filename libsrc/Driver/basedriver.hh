#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "General/environment.hh"
#include "Domain/domain.hh"

namespace CoupledField
{

//! a base class for driving classes where we implemented time-stepping
class BaseDriver
{
public:
  //! constructor
  /*!
    \param adomain pointer to class Domain
  */
  BaseDriver(Domain * adomain);

   //! deconstructor
  virtual ~BaseDriver();

  //! set the pdes, which have to be solved
  void SetPDEs(StdVector<BasePDE*> & pdes);
  
  //! main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem()=0;

  //! to setup matrices of PDE. we call according method of class PDE for setup matrices of PDE in assembling procedure.
  /*!
    \param pdenumber number of PDE
    \param matrixtype type of matrix
  */
  virtual void SetupMatricesPDE(Integer pdenumber, const Integer matrixtype)
  { Error("SetupMatricesPDE not implemented in base class!",__FILE__,__LINE__); };
  
 protected:
  //! pointer to class Domain
  Domain * ptdomain_;

  //! vector of PDEs to solve
  StdVector<BasePDE*> pdes_; 

  //! --------------------- stuff for computation with adaptivity
  //! for printing a sequence of files in dir meshes in gmv-format
  WriteResults * ptMeshes_;

  //! counter of meshes for printing meshes
  Integer nummeshes_;  

  //! print mesh in special file. this method is used in adaptive procedure for space.
  void PrintSeqMeshes();

  //! auxiliary function for computation with adaptivity: open files for printing sequence of refined meshes with error map 
  Boolean printMeshesOrNot();
  
private:
  
};

}

#endif // FILE_DRIVER
