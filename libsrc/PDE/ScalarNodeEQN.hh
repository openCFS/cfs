#ifndef FILE_SCALARNODEEQN_2003
#define FILE_SCALARNODEEQN_2003

#include "BaseEQN.hh"
#include <Utils/vector.hh>

namespace CoupledField
{

//! a base class for quation data handling

class ScalarNodeEQN : public BaseEQN
{

public:
  //! constructor
  ScalarNodeEQN(Grid * aptgrid, BCs *aptBCs, std::vector<std::string>& asubdoms, 
		std::vector<std::string>& bcs_hd, Integer actlevel);

  //! deconstructor
  virtual ~ScalarNodeEQN();
  
  //! initilization
  virtual void InitEQN();

  //! account for homogeneous Dirichlet b.c.
  virtual void SetHomoDBCs();

  //! applys the equation numbers
  virtual void ApplyEqnNrs();

  //!retruns the equations numbers in Eqns for given node numbers nodes
  virtual void Mesh2Eqn(Vector<Integer>& Eqns, Vector<Integer>& nodes);

  //! Print Equation-Array
  virtual void Print();

protected:


private:

  Vector<Integer> EQN_;  //!< contains the euqation numbers   
  Vector<Integer> EQN2MeshNode_; //!< mapping from equation number to node number
  std::vector<std::string> bcs_hd_; //!< homogeneous Dirichlet b.c.
};

}

#endif // FILE_SCALARNODEEQN
