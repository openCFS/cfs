#include "scalarnodeEQN.hh"

#include <iomanip>

namespace CoupledField
{
  
ScalarNodeEQN::ScalarNodeEQN(Grid * aptGrid, 
			     BCs * aptBCs,
			     StdVector<std::string>& asubdoms, 
			     Integer actlevel, 
			     Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "NodeEQN::NodeEQN" );
  
  isBlockMapped_ = FALSE;
  dofsPerEQN_ = 1;
 
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
  CalcLocalGlobalMapping(mesh2PDENode_,
			 pde2MeshNode_,
			 mesh2PDEElem_,
			 pde2MeshElem_);
 
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
  std::string warnMsg, errMsg;
  Integer notIncludedBCs = 0;
  pdeNode2EQN_.Clear();
  eqn2Pos_.Clear();
  pdeNode2EQN_.Resize(numPDENodes_);
  pdeNode2EQN_.Init(1);
 

  // STEP 2
  StdVector<Integer> countNodes;
  countNodes.Resize(numPDENodes_);
  for (Integer i=0; i<homoDirichletNodes_.GetSize(); i++)
    {
      if (mesh2PDENode_[homoDirichletNodes_[i]-1]-1 < 0) {
	  warnMsg  = "ScalarNodeEQN::CalcMapping: Homogen. Dirichlet node nr. ";
	  warnMsg += Info->GenStr(homoDirichletNodes_[i]);
	  warnMsg += " is not contained in any of the regions for this PDE";
	  Info->Warning(warnMsg, __FILE__, __LINE__);
	  notIncludedBCs++;
	} 
      else if (countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1] != 0) {
	errMsg  = "ScalarNodeEQN::CalcMapping: HomDirihchletNode Nr. ";
	errMsg += Info->GenStr(homoDirichletNodes_[i]);
	errMsg += "\n occured already at least one time in the list of boundary nodes ";
	errMsg += "for this PDE!\n Please check, if this node is defined in";
	errMsg += " more than one level of boundary nodes!";
	Error(errMsg.c_str(), __FILE__, __LINE__);
      }
      else {
	pdeNode2EQN_[mesh2PDENode_[homoDirichletNodes_[i]-1]-1] = 0;
        countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]++;
      }
    }
  
//   std::cerr << "after step2 " << std::endl;
  
  // STEP 3
  for (Integer i=0; i<constraintSlaveNodes_.GetSize(); i++)
    pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] = 0;


  eqn2Pos_.Resize(numPDENodes_ 
		  - homoDirichletNodes_.GetSize()
		  - constraintSlaveNodes_.GetSize()
		  + notIncludedBCs); 


//   std::cerr << "after step3" << std::endl;
  
  // STEP 4
  for (Integer i=0; i<pde2MeshNode_.GetSize(); i++)
      if (pdeNode2EQN_[i] != 0)
	{
	  eqnCounter ++;
	  pdeNode2EQN_[i] = eqnCounter;
	  eqn2Pos_[eqnCounter-1] = pde2MeshNode_[i]-1;
	}
  
//   std::cerr << "after step4" << std::endl;

  // STEP 5
  for (Integer i=0; i<constraintSlaveNodes_.GetSize(); i++)
    pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] =
      pdeNode2EQN_[mesh2PDENode_[constraintMasterNodes_[i]-1]-1];
  
//   std::cerr << "after step5" << std::endl;
  // Now object is initialized
  isInitialized_ = TRUE;
  numEqns_ = eqn2Pos_.GetSize();

  numBuildInDirichletEQNs_ = numPDENodes_ - numEqns_;
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
  

  for (Integer i=0; i<pde2MeshNode_.GetSize(); i++)
    {
      out << std::setw(10) << i+1  << " | ";
      out << std::setw(13) << pde2MeshNode_[i] << " | ";
      out << std::setw(10) << pdeNode2EQN_[i] << std::endl;
    }
}


void ScalarNodeEQN::EQN2SolVectorPos(const StdVector<Integer> &eqnNr, 
				     StdVector<Integer> &pos) const
{
  ENTER_FCN( "ScalarNodeEQN::EQN2SolVectorPos" );
  Error( "Not implemented" );
}


void ScalarNodeEQN::Node2EQN(const Integer nodeNr, 
			     const Integer dof,
			     Integer & eqnNr,
			     Integer & eqnDof) const
{ 
  ENTER_FCN( "ScalarNodeEQN::Node2EQN" );
#ifdef CHECK_INDEX
  if (nodeNr > mesh2PDENode_.GetSize())
    Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
	  __FILE__, __LINE__);
#endif
  eqnNr = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1];
  eqnDof = 1;
}


void ScalarNodeEQN::Node2EQN(const Integer nodeNr, StdVector<Integer> &eqnNr) const
{
  ENTER_FCN( "ScalarNodeEQN::Node2EQN" );

#ifdef CHECK_INDEX
  if (nodeNr > mesh2PDENode_.GetSize())
    Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
	  __FILE__, __LINE__);
#endif 

  eqnNr.Resize(dofsPerNode_);

  eqnNr[0] = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1];
}
  

void ScalarNodeEQN::Node2EQN(const StdVector<Integer> &nodeNr,
			     StdVector<Integer> &eqnNr) const
{
  ENTER_FCN( "ScalarNodeEQN::Node2EQN" );
  
  eqnNr.Resize(nodeNr.GetSize());

  for (Integer i=0; i<nodeNr.GetSize(); i++)
    {
#ifdef CHECK_INDEX
      if (nodeNr[i] > mesh2PDENode_.GetSize())
	Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
	      __FILE__, __LINE__);
#endif
      eqnNr[i] =  pdeNode2EQN_[mesh2PDENode_[nodeNr[i]-1]-1];
    }
}

} // end of namespace
