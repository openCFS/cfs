#ifndef FILE_BASEEQN_2004
#define FILE_BASEEQN_2004

#include "General/environment.hh"
#include "Domain/grid.hh"
#include "Domain/bcs.hh"
#include <vector>

namespace CoupledField
{

//! a base class for equation data handling

class BaseEQN
{
public:

  //! Constructor
  BaseEQN(Grid * aptGrid,
	  BCs * aptBCs,
	  std::vector<std::string>& asubdoms, 
	  Integer actlevel, 
	  Integer dofsPerNode);

  //! Destructor
  virtual ~BaseEQN();
  
  //! Return the number of equations
  //! (= number of nodes * total number Dofs <br>
  //! - hom. Dirichlet nodes <br>
  //! - constraint nodes)
  inline Integer GetNumEQNs() const {return numEqns_;}
  
  //! Returns true, if nodal mapping is applied
  inline Boolean IsNodalMapped() const {return isNodalMapped_;}

  //! Returns true, if block mapping is applied
  inline Boolean IsBlockMapped() const {return isBlockMapped_;}

  //! Set the node numbers and dofs 
  //! with homogeneous Dirchlet BC
  virtual void SetHomoDirichletBCs(const std::vector<std::string> & nodeNrs,
				   const std::vector<std::string> & dofs);

  //! Set the node numbers and dofs which are
  //! slave nodes w.r.t. to constraints
  virtual void SetConstraints(const std::vector<Integer> & slaveNodeNrs,
			      const std::vector<Integer> & masterodeNrs,
			      const std::vector<std::string> & dofs);

  //! Calculate the mapping after Dirichlet and
  //! constraint nodes were set
  virtual void CalcMapping() = 0;

  //! Print the mapping nodes <->EQNs
  virtual void Print(std::ostream & out) const = 0;
  
protected:

  //! Flag for initialisation
  Boolean isInitialized_;
  
  //! Flag for indicating blockwise numbering
  Boolean isBlockMapped_;

  //! Flag for indicating nodal numbering vs. edge numbering
  Boolean isNodalMapped_;
  
  //! Pointer to Grid
  Grid * ptGrid_;  

  //! Pointer to BCs
  BCs * ptBCs_;

  //! Current refinement level (adaptivity)
  Integer actlevel_;

  //! Number of nodes in PDE
  Integer numPDENodes_;

  //! Number of dofs per node
  Integer dofsPerNode_;
 
  //! Number of equations in PDE
  Integer numEqns_;
  
   //! Vector containing subdomain names
  //! of according PDE
  std::vector<std::string> subdoms_;  

  //! Vector containing dirichlet Nodes
  std::vector<Integer> homoDirichletNodes_;

  //! Vector containing dofs associated
  //! with hom. Dirichlet nodes
  std::vector<Integer> homoDirichletDofs_;
  
  //! Vector containing constraint master nodes
  std::vector<Integer> constraintMasterNodes_;

  //! Vector containing according slave nodes
  std::vector<Integer> constraintSlaveNodes_;
  
  //! Vecotr containing dofs associated
  //! with slave nodes
  std::vector<Integer> constraintDofs_;
  
  //! Mapping of string-dof to integer (1-based)
  Integer GetBCDof(const std::string dofString) const;

private:
   


};

}

#endif // FILE_BASEEQN
