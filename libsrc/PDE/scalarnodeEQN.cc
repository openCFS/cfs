#include "scalarnodeEQN.hh"
#include <iomanip>

namespace CoupledField
{
  
ScalarNodeEQN::ScalarNodeEQN(Grid * aptGrid, 
			     BCs * aptBCs,
			     std::vector<std::string>& asubdoms, 
			     Integer actlevel, 
			     Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "NodeEQN::NodeEQN" );
  
  isBlockMapped_ = FALSE;
 
}

ScalarNodeEQN::~ScalarNodeEQN()
{
 
  ENTER_FCN( "NodeEQN::NodeEQN" );
 
}


void ScalarNodeEQN::CalcMapping()
{

  ENTER_FCN( "ScalarNodeEQN::CalcMapping" );
 
  // First apply Mapping from global to
  // local node numbers and back
  CalcMesh2PDENode(mesh2PDENode_, pde2MeshNode_);
 
  // Idea of the algorithm:
  // Step 1: Initialize pdeNode2eqn_ with 1
  // Step 2: For each entry in homoDirichletNodes_ set the according
  //         entry in pdeNode2eqn_ to 0
  // Step 3: For each entry in constraintSlaveNodes_ set the according
  //         entry in pdeNode2eqn_ to 0
  // Step 4: Loop over all entries in pde2Meshnode
  //         and assign each non-zero entry a equation number
  // Step 5: Afterwards loop again over all nodes in constraintSlaveNodes_
  //         and set the according entry in pdeNode2EQN_ to the
  //         value of constraintMasterNode

  Integer eqnCounter = 0;;

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
  pdeNode2EQN_.resize(numPDENodes_,1);
  eqn2Pos_.resize(numPDENodes_ 
		  - homoDirichletNodes_.size()
		  - constraintSlaveNodes_.size()); 

  // STEP 2
  for (Integer i=0; i<homoDirichletNodes_.size(); i++)
    pdeNode2EQN_[mesh2PDENode_[homoDirichletNodes_[i]-1]-1] = 0;
  
  // STEP 3
  for (Integer i=0; i<constraintSlaveNodes_.size(); i++)
    pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] = 0;

 

  // STEP 4
  for (Integer i=0; i<pde2MeshNode_.size(); i++)
      if (pdeNode2EQN_[i] != 0)
	{
	  eqnCounter ++;
	  pdeNode2EQN_[i] = eqnCounter;
	  eqn2Pos_[eqnCounter-1] = pde2MeshNode_[i];
	}
  

  // STEP 5
  for (Integer i=0; i<constraintSlaveNodes_.size(); i++)
    pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] =
      pdeNode2EQN_[mesh2PDENode_[constraintMasterNodes_[i]-1]-1];
  

  // Now object is initialized
  isInitialized_ = TRUE;
  numEqns_ = eqn2Pos_.size();
}

void ScalarNodeEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "ScalarNodeEQN::Print" );
  
  out << "Equation numbering - Information" << std::endl;
  out << "================================" << std::endl;
  out << "DOFs per Node: 1 " << std::endl;
  out << "Using SCALAR numbering of equations" << std::endl;
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


void ScalarNodeEQN::EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
				     std::vector<Integer> &pos) const
{
  ENTER_FCN( "ScalarNodeEQN::EQN2SolVectorPos" );
  Info->Error( "Not implemented" );
}


void ScalarNodeEQN::Node2EQN(const Integer nodeNr, std::vector<Integer> &eqnNr) const
{
  ENTER_FCN( "ScalarNodeEQN::Node2EQN" );
  
  eqnNr.resize(dofsPerNode_);

  eqnNr[0] = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1];
}
  

void ScalarNodeEQN::Node2EQN(const std::vector<Integer> &nodeNr,
			     std::vector<Integer> &eqnNr) const
{
  ENTER_FCN( "ScalarNodeEQN::Node2EQN" );
  
  eqnNr.resize(nodeNr.size());

  for (Integer i=0; i<nodeNr.size(); i++)
       eqnNr[i] =  pdeNode2EQN_[mesh2PDENode_[nodeNr[i]-1]-1];
}

} // end of namespace
