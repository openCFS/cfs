#include "blocknodeEQN.hh"
#include <iomanip>

namespace CoupledField
{
  
BlockNodeEQN::BlockNodeEQN(Grid * aptGrid, 
			   BCs * aptBCs,
			   std::vector<std::string>& asubdoms, 
			   Integer actlevel, 
			   Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "BlockNodeEQN::BlockNodeEQN" );

  isBlockMapped_ = TRUE;
}

BlockNodeEQN::~BlockNodeEQN()
{
 
  ENTER_FCN( "BlockNodeEQN::BlockNodeEQN" );
}


void BlockNodeEQN::CalcMapping()
{

  ENTER_FCN( "BlockNodeEQN::CalcMapping" );
 
  // First apply Mapping from global to
  // local node numbers and back
  CalcMesh2PDENode(mesh2PDENode_, pde2MeshNode_);
 
 // Idea of the algorithm:
  // Step 1: Initialize pdeNode2eqn_ with 1
  // Step 2: Count, how many Dofs of every node are a hom. BC
  //         -> only if all Dofs are zero, the whole block can
  //         be omitted
  // Step 3: Count, how many Dofs of every node are a constraint slave node
  //         -> only if all Dofs are zero, the whole block can
  //
  // Step 4: Loop over all entries in pde2Meshnode
  //         and assign only a a equation number if not all dofs of one node
  //         are a hom. Dirichlet BC or a constraint slave node
  // Step 5: Afterwards loop again over all nodes in constraintSlaveNodes_
  //         and set the according entry in pdeNode2EQN_ to the
  //         value of constraintMasterNode

  Integer eqnCounter = 0;

  // STEP 1

  // Each local PDE-node number is assigned a 
  // corresponding EQN number:
  //
  // 0      : homog. Dirichlet BC
  // +number: normal Equation
  // -number: slave node (constraint)
  //
  // Therefore both pdeNode2EQN has the same size
  // as number of PDE nodes.
  // Furthermore, eqn2Pos has the length of
  // total number of equations, which is the total
  // number of pdenodes minus the number of hom.
  // Dirichlet nodes - the number of constraintNodes
  pdeNode2EQN_.clear();
  eqn2Pos_.clear();
  pdeNode2EQN_.resize(numPDENodes_,0);
  std::vector<Integer> eqn2Pos_Temp;
  eqn2Pos_Temp.reserve(numPDENodes_);

  // STEP 2
  // Check if there exist nodes, which only have
  // hom. Dirichlet BC dof
  std::vector<Integer> numDirichletDofsPerNode;
  numDirichletDofsPerNode.resize(numPDENodes_,0);
  for (Integer i=0; i<homoDirichletNodes_.size(); i++)
    numDirichletDofsPerNode[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]++;

  // STEP 3
  // Check if there are constraint slavenodes, where
  // all dofs depend on the same dofs of the same
  // master node
  std::vector<Integer> numConstraintDofsPerNode;
  std::vector<Integer> masterNodes;
  numConstraintDofsPerNode.resize(numPDENodes_);
  masterNodes.resize(numPDENodes_,0);

  for (Integer i=0; i<constraintSlaveNodes_.size(); i++)
    if (masterNodes[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] == 0)
      {
	// If masternodes are still empty
	// -> Set masternode 
	masterNodes[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] =
	  mesh2PDENode_[constraintMasterNodes_[i]-1];
	// increment number of constraint dofs for this node
	numConstraintDofsPerNode[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]++;
	// set the eqn-number of this node to -masternode number
	// NOTE: this will be overwritten in STEP 4, if not all
	// dofs of this node depend on the same master node
	pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] =
	  - mesh2PDENode_[constraintMasterNodes_[i]-1]; 
      }
    else if (masterNodes[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] ==
	     mesh2PDENode_[constraintMasterNodes_[i]-1])
      // If old masternode and new are the same
      // -> only increment number of constraint dofs for this node
      numConstraintDofsPerNode[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]++;
			       
  
  // STEP 4
  for (Integer i=0; i<pde2MeshNode_.size(); i++)
    if (numDirichletDofsPerNode[i] != dofsPerNode_ &&
	numConstraintDofsPerNode[i] != dofsPerNode_)
      {
	eqnCounter ++;
	pdeNode2EQN_[i] = eqnCounter;
	  eqn2Pos_Temp.push_back(pde2MeshNode_[i]);
      }	
  
  
  // Now object is initialized
  numDirichletDofsPerNode.clear();
  numConstraintDofsPerNode.clear();
  masterNodes.clear();
  eqn2Pos_ = eqn2Pos_Temp;
  eqn2Pos_Temp.clear();
  isInitialized_ = TRUE;
  numEqns_ = eqn2Pos_.size();
}

void BlockNodeEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "BlockNodeEQN::Print" );

  out << "Equation numbering - Information" << std::endl;
  out << "================================" << std::endl;
  out << "DOFs per Node: " << dofsPerNode_ << std::endl;
  out << "Using BLOCK numbering of equations" << std::endl;
  out << std::endl;

  // Print pde2MeshNode_ and pdeNode2EQN_
  out << std::setw(10) << "PDE NodeNr" << " | ";
  out << std::setw(13) << "Global NodeNr" << " | ";
  out << std::setw(10) << "EQ number";
  out << std::endl;
  out << std::setfill('-') << std::setw(39) << "-" << std::endl;
  out << std::setfill(' ');
  

  for (Integer i=0; i<pde2MeshNode_.size(); i++)
    {
      out << std::setw(10) << i+1  << " | ";
      out << std::setw(13) << pde2MeshNode_[i] << " | ";
      out << std::setw(10) << pdeNode2EQN_[i] << std::endl;
    }
}


void BlockNodeEQN::EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
				     std::vector<Integer> &pos) const
{
  ENTER_FCN( "BlockNodeEQN::EQN2SolVectorPos" );
  Info->Error( "Not implemented" );
}


void BlockNodeEQN::Node2EQN(const Integer nodeNr, std::vector<Integer> &eqns) const
{
  ENTER_FCN( "BlockNodeEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}
  

void BlockNodeEQN::Node2EQN(const std::vector<Integer> &nodeNr,
			     std::vector<Integer> &eqnNr) const
{
  ENTER_FCN( "BlockNodeEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}

} // end of namespace
