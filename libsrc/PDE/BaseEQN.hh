#ifndef FILE_BASEEQN_2003
#define FILE_BASEEQN_2003

#include <General/environment.hh>
#include <Domain/grid.hh>
#include <Domain/bcs.hh>

namespace CoupledField
{

//! a base class for quation data handling

class BaseEQN
{
public:
  //! constructor
  BaseEQN(Grid * aptgrid, BCs *aptBCs, std::vector<std::string>& asubdoms, Integer actlevel);

  //! deconstructor
  virtual ~BaseEQN();
  
  //! initilization
  virtual void InitEQN()=0;

  //! account for homogeneous Dirichlet b.c.
  virtual void SetHomoDBCs()=0;

  //! applys the equation numbers
  virtual void ApplyEqnNrs()=0;

  //!retruns the equations numbers in Eqns for given node numbers nodes
  virtual void Mesh2Eqn(Vector<Integer>& Eqns, Vector<Integer>& nodes)=0;

  //! Print Equation-Array
  virtual void Print()=0;

protected:

  Grid * ptgrid_;  //!< pointer to Grid
  BCs * ptbcs_;   //!< pointer to boundary condition object
  std::vector<std::string> subdoms_;  //!< subdomain-levels belongig to PDE

  Integer actlevel_;
  Integer numPDENodes_;
  Integer numEqns_;

private:
   


};

}

#endif // FILE_BASEEQN
