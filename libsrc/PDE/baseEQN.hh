#ifndef FILE_BASEEQN_2003
#define FILE_BASEEQN_2003

#include "General/environment.hh"
#include "Domain/grid.hh"
#include <vector>

namespace CoupledField
{

//! a base class for equation data handling

class BaseEQN
{
public:

  //! Constructor
  BaseEQN(Grid * aptgrid, 
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
  
  //! Returns true, if block mapping is applied
  inline Boolean IsBlockMapped() const {return isBlockMapped_;}

  //! Set the node numbers and dofs 
  //! with homogeneous Dirchlet BC
  virtual void SetHomDirichletBCs(const std::vector<Integer> &nodeNrs,
				  const std::vector<std::string> &dofs);

  //! Set the node numbers and dofs which are
  //! slave nodes w.r.t. to constraints
  virtual void SetConstraints(const std::vector<Integer> &nodeNrs,
			      const std::vector<std::string> &dofs);

  //! Calculate the mapping after Dirichlet and
  //! constraint nodes were set
  /*!
    \param isBlock (input) If false, each block gets dofsPerNode
    equation numbers instead of only one
    \param isSuperBlock (input) If true, at each dof gets continuous 
    equation numbers, which results in a supberBlock system
  */
  //! \deprecated 'isSuperBlock' won't be needed when the complete
  //! Matrixcoupling structure is implemented
  virtual void CalcMapping() = 0;

  //! Print the mapping nodes <->EQNs
  virtual void Print(std::ostream & out) const = 0;
  
protected:

  //! Flag for initialisation
  Boolean isInitialized_;
  
  //! Flag for indicating blockwise numbering
  Boolean isBlockMapped_;
  
  //! Pointer to Grid
  Grid * ptgrid_;  

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
  std::vector<Integer> dirichletNodes_;

  //! Vector containing dofs associated
  //! with hom. Dirichlet nodes
  std::vector<Integer> dirichletDofs_;
  
  //! Vector containing constraint-slave nodes
  std::vector<Integer> constraintNodes_;

  //! Vecotr containing dofs associated
  //! with slave nodes
  std::vector<std::string> constraintDofs_;
  
  //! Mapping of string-dof to integer (0-based)
  Integer GetBCDof(const std::string dofString) const;

private:
   


};

}

#endif // FILE_BASEEQN
