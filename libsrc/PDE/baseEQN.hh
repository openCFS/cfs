#ifndef FILE_BASEEQN_2004
#define FILE_BASEEQN_2004

#include "General/environment.hh"
#include "Domain/grid.hh"
#include "Domain/bcs.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

//! a base class for equation data handling

class BaseEQN
{
public:

  //! Constructor
  BaseEQN(Grid * aptGrid,
	  BCs * aptBCs,
	  StdVector<std::string>& asubdoms, 
	  Integer actlevel, 
	  Integer dofsPerNode);

  //! Destructor
  virtual ~BaseEQN();
  
  //! Return the number of equations
  //! (= number of nodes * total number Dofs <br>
  //! - hom. Dirichlet nodes <br>
  //! - constraint nodes)
  inline Integer GetNumEQNs() const {return numEqns_;}
  
  //! Return number of degree of freedoms per node
  inline Integer GetNumDofs() const {return dofsPerNode_;}

  //! Return number of dofs per equation
  inline Integer GetNumDofsPerEQN() const {return dofsPerEQN_;}
  
  //! Return number of build in Dirichlet values

  //! Returns the number of Dirichlet values that were assigned
  //! a 0 and therefore have not to set explicitly to 
  //! the algebraic system
  inline Integer GetNumBuildInDirichletEQNs() const
  {return numBuildInDirichletEQNs_;}

  //! Returns true, if nodal mapping is applied
  inline Boolean IsNodalMapped() const {return isNodalMapped_;}

  //! Returns true, if block mapping is applied
  inline Boolean IsBlockMapped() const {return isBlockMapped_;}

  //! Set the node numbers and dofs 
  //! with homogeneous Dirchlet BC
  virtual void SetHomoDirichletBCs(const StdVector<std::string> & nodeNrs,
				   const StdVector<std::string> & dofs);

  //! Set the node numbers and dofs which are
  //! slave nodes w.r.t. to constraints
  virtual void SetConstraints(const StdVector<Integer> & slaveNodeNrs,
			      const StdVector<Integer> & masterodeNrs,
			      const StdVector<std::string> & dofs);
			      

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
 
  //! Number of dofs per eqn
  Integer dofsPerEQN_;

  //! Number of equations in PDE
  Integer numEqns_;
  
  //! Number of Dirichlet values
  //! which have been eliminated
  //! by equation class
  Integer numBuildInDirichletEQNs_;

   //! Vector containing subdomain names
  //! of according PDE
  StdVector<std::string> subdoms_;  

  //! Vector containing dirichlet Nodes
  StdVector<Integer> homoDirichletNodes_;

  //! Vector containing dofs associated
  //! with hom. Dirichlet nodes
  StdVector<Integer> homoDirichletDofs_;
  
  //! Vector containing constraint master nodes
  StdVector<Integer> constraintMasterNodes_;

  //! Vector containing according slave nodes
  StdVector<Integer> constraintSlaveNodes_;
  
  //! Vector containing dofs associated
  //! with slave nodes
  StdVector<Integer> constraintDofs_;

  //! Mapping of string-dof to integer (1-based)
  Integer GetBCDof(const std::string dofString) const;

private:
   


};

}

#endif // FILE_BASEEQN
