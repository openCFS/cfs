#include "superblockEQN.hh"

#include <iomanip>

namespace CoupledField
{
  
SuperBlockEQN::SuperBlockEQN(Grid * aptGrid, 
			     BCs * aptBCs,
			     StdVector<std::string>& asubdoms, 
			     Integer actlevel, 
			     Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "SuperBlockEQN::SuperBlockEQN" );
  isBlockMapped_ = FALSE;
  dofsPerEQN_ = 1;
}

SuperBlockEQN::~SuperBlockEQN()
{
 
  ENTER_FCN( "SuperBlockEQN::SuperBlockEQN" );
}


void SuperBlockEQN::CalcMapping()
{

  ENTER_FCN( "SuperBlockEQN::CalcMapping" );

  // First apply Mapping from global to
  // local node/elem numbers and back
  CalcLocalGlobalMapping(mesh2PDENode_,
			 pde2MeshNode_,
			 mesh2PDEElem_,
			 pde2MeshElem_); 

  Integer numMechBCs = 0;


  // 1. STEP: Count number of mechanic eqns
  for (Integer i=0; i<homoDirichletNodes_.GetSize(); i++)
    // check if node has a mechanic dof
    if(homoDirichletDofs_[i] < dofsPerNode_)
      numMechBCs++;

  
  numMechEQNs_ = pde2MeshNode_.GetSize() * (dofsPerNode_-1) - numMechBCs;
  numElecEQNs_ = pde2MeshNode_.GetSize() - 
    (homoDirichletNodes_.GetSize() - numMechBCs);

  std::cerr << "We habe " << pde2MeshNode_.GetSize() << " local nodes " << std::endl;
  std::cerr << "We have " << numMechEQNs_ << " mechanic EQNs" << std::endl;
  std::cerr << "We have " << numElecEQNs_ << " electric EQNs" << std::endl;
  
  mechNode2EQN_.Resize(pde2MeshNode_.GetSize(), dofsPerNode_-1);
  mechNode2EQN_.Init(1);
  elecNode2EQN_.Resize(pde2MeshNode_.GetSize());
  elecNode2EQN_.Init(1);
  eqn2Pos_.Resize(numPDENodes_ * dofsPerNode_
		  - homoDirichletNodes_.GetSize()
		  - constraintSlaveNodes_.GetSize());
  
  std::cerr << "size of eqn2Pos_ = " << eqn2Pos_.GetSize() << std::endl;
   
  // STEP 2
  for (Integer i=0; i<homoDirichletNodes_.GetSize(); i++)
    {
      if (homoDirichletDofs_[i] == dofsPerNode_)
	elecNode2EQN_[mesh2PDENode_[homoDirichletNodes_[i]-1]-1] = 0;
      else
	mechNode2EQN_[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
	  [homoDirichletDofs_[i]-1] = 0;
    }
  std::cerr << "After setting homoDirichletNodes" << std::endl;

  // STEP 3
  for (Integer i=0; i<constraintSlaveNodes_.GetSize(); i++)
    if (homoDirichletDofs_[i] == dofsPerNode_)
       elecNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] = 0;
    else
      mechNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]
	[constraintDofs_[i]-1] = 0;
  
  std::cerr << "After setting constraints" << std::endl;

  // STEP 4

  Integer mechEQNCounter = 0;
  Integer elecEQNCounter = numMechEQNs_;
  for (Integer iNode=0; iNode<pde2MeshNode_.GetSize(); iNode++)
    {
      std::cerr << "Checking Local Node " << iNode << std::endl;
      // Assign mechanic equation numbers
      for (Integer iDof=0; iDof<dofsPerNode_-1; iDof++)
	{
	  if (mechNode2EQN_[iNode][iDof] != 0)
	    {
	      mechEQNCounter++;
	      mechNode2EQN_[iNode][iDof] = mechEQNCounter;
	      std::cerr << "mechNode2EQN_[" << iNode << "][" << iDof << " = " << mechEQNCounter << std::endl;
	      eqn2Pos_[mechEQNCounter-1] = 
	      (pde2MeshNode_[iNode]-1)*dofsPerNode_ + iDof;
	      std::cerr << "eqnPos_[" << mechEQNCounter-1 << "] = " << (pde2MeshNode_[iNode]-1)*dofsPerNode_ + iDof << std::endl;
	      
	    }
	}
	  // Assign electric equation numbers
	  if (elecNode2EQN_[iNode] != 0)
	    {
	      elecEQNCounter++;
	      elecNode2EQN_[iNode] = elecEQNCounter;
	       std::cerr << "elecNode2EQN_[" << iNode << "] = " << elecEQNCounter << std::endl;
	      eqn2Pos_[elecEQNCounter-1] = 
		(pde2MeshNode_[iNode]-1)*dofsPerNode_ + dofsPerNode_ -1;
	      std::cerr << "eqnPos_[" << elecEQNCounter-1 << "] = " << (pde2MeshNode_[iNode]-1)*dofsPerNode_  + dofsPerNode_ -1  << std::endl;
	    }
    }
  
      
  // STEP 5
  for (Integer i=0; i<constraintSlaveNodes_.GetSize(); i++)
    if (constraintDofs_[i] == dofsPerNode_)
      elecNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] =
	elecNode2EQN_[mesh2PDENode_[constraintMasterNodes_[i]-1]-1];
    else
      mechNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]
	[constraintDofs_[i]-1] =
	mechNode2EQN_[mesh2PDENode_[constraintMasterNodes_[i]-1]-1]
	[constraintDofs_[i]-1];
  
  // Now object is initialized
  isInitialized_ = TRUE;
  numEqns_ = numMechEQNs_ + numElecEQNs_;

  numBuildInDirichletEQNs_ = numPDENodes_ * dofsPerNode_ - numEqns_;
  std::cerr << "size of homoDirichletNodes = " << homoDirichletNodes_.GetSize() << std::endl;
  std::cerr << "Size of numBuildInDirichletEQNs_ = " << numBuildInDirichletEQNs_ << std::endl;
}

void SuperBlockEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "SuperBlockEQN::Print" );
  out << "Equation numbering - Information" << std::endl;
  out << "================================" << std::endl;
  out << "DOFs per Node: " << dofsPerNode_ << std::endl;
  out << "Using SUPER-BLOCK numbering of equations" << std::endl;
  out << "Number of MECHANIC EQNs: " << numMechEQNs_ << std::endl;
  out << "Number of ELECTRIC EQNs: " << numElecEQNs_ << std::endl;
  out << std::endl;

  // Print pde2MeshNode_ and pdeNode2EQN_
  out << std::setw(10) << "PDE NodeNr" << " | ";
  out << std::setw(13) << "Global NodeNr" << " | ";
  out << std::setw(10) << "EQ number";
  out << std::endl;
  out << std::setfill('-') << std::setw(39) << "-" << std::endl;
  out << std::setfill(' ');
  
  for (Integer iNode=0; iNode<pde2MeshNode_.GetSize(); iNode++)
    {
      // Print out first dof
      out << std::setw(10) << iNode+1  << " | ";
      out << std::setw(13) << pde2MeshNode_[iNode] << " | ";
      out << std::setw(10) << mechNode2EQN_[iNode][0] << std::endl;
      
      // Print out rest of mechanic dofs
      for (Integer iDof = 1; iDof < dofsPerNode_-1; iDof++)
	{
	  out << std::setw(10) << " " << " | ";
	  out << std::setw(13) << " " << " | ";
	  out << std::setw(10) << mechNode2EQN_[iNode][iDof] << std::endl;
	}
      // Print out electric dof
      out << std::setw(10) << " " << " | ";
      out << std::setw(13) << " " << " | ";
      out << std::setw(10) << elecNode2EQN_[iNode] << std::endl;
    }
  
}


void SuperBlockEQN::EQN2SolVectorPos(const StdVector<Integer> &eqnNr, 
				     StdVector<Integer> &pos) const
{
  ENTER_FCN( "SuperBlockEQN::EQN2SolVectorPos" );
  Error( "Not implemented" );
}

void SuperBlockEQN::Node2EQN(const Integer nodeNr, 
			     const Integer dof,
			     Integer & eqnNr,
			     Integer & eqnDof) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );
#ifdef CHECK_INDEX
  if (nodeNr > mesh2PDENode_.GetSize())
    Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
	  __FILE__, __LINE__);
#endif
 
  if (dof == dofsPerNode_)
    eqnNr = elecNode2EQN_[mesh2PDENode_[nodeNr-1]-1];
  else
    eqnNr = mechNode2EQN_[mesh2PDENode_[nodeNr-1]-1][dof-1];
  
  eqnDof = 1;
    
}

void SuperBlockEQN::Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );
  Error( "Not implemented" );
}
  

void SuperBlockEQN::Node2EQN(const StdVector<Integer> &nodeNr,
			     StdVector<Integer> &eqnNr) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );

  eqnNr.Resize(nodeNr.GetSize()*dofsPerNode_);

  for(Integer iNode=0; iNode<nodeNr.GetSize(); iNode++)
    {
#ifdef CHECK_INDEX
      if (nodeNr[iNode] > mesh2PDENode_.GetSize())
	Error("SuperBlockEQN::Node2EQN: Index out of bounds", 
	      __FILE__, __LINE__);
#endif

      // Set mechanic eqnNrs
      for (Integer iDof=0; iDof<dofsPerNode_-1; iDof++)
	eqnNr[iNode*dofsPerNode_ + iDof] =
	  mechNode2EQN_[mesh2PDENode_[nodeNr[iNode]-1]-1][iDof];

      // Set electric eqnNrs
      
      eqnNr[iNode*dofsPerNode_ + dofsPerNode_ -1] =
	elecNode2EQN_[mesh2PDENode_[nodeNr[iNode]-1]-1];
    }
}

} // end of namespace
