#include "scalarblockEQN.hh"

#include <iomanip>

namespace CoupledField
{
  
ScalarBlockEQN::ScalarBlockEQN(Grid * aptGrid, 
			       BCs * aptBCs,
			       StdVector<std::string>& asubdoms, 
			       Integer actlevel, 
			       Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
  isBlockMapped_ = FALSE;
  dofsPerEQN_ = 1;
  
}

ScalarBlockEQN::~ScalarBlockEQN()
{
 
  ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
}


void ScalarBlockEQN::CalcMapping()
{

  ENTER_FCN( "ScalarBlockEQN::CalcMapping" );

  // First apply Mapping from global to
  // local node/elem numbers and back
  CalcLocalGlobalMapping(mesh2PDENode_,
			 pde2MeshNode_,
			 mesh2PDEElem_,
			 pde2MeshElem_);
 
  Integer eqnCounter = 0;
  std::string warnMsg, errMsg;
//  std::cerr << "Content of homoDirichletNodes: " << std::endl;
//   std::cerr << "-------------------------------" << std::endl;
//   std::cerr << homoDirichletNodes_ << std::endl;

//   std::cerr << "Content of homoDirichletDofs : " << std::endl;
//   std::cerr << "-------------------------------" << std::endl;
//   std::cerr << homoDirichletDofs_ << std::endl;

//   std::cerr << "we have " << numPDENodes_  << " PDE nodes" << std::endl;
//   std::cerr << "size of Mesh2PDENode " << mesh2PDENode_.GetSize() << std::endl;
//   std::cerr << "size of PDE2MeshNode " << pde2MeshNode_.GetSize() << std::endl;
//   std::cerr << "We have " << homoDirichletNodes_.GetSize() << " homo Dirichlet nodes" << std::endl;
//   std::cerr << "pdeNode2EQN.RowSize = " << numPDENodes_ << std::endl;
//   std::cerr << "pdeNode2EQN.ColSize = " << numPDENodes_ << std::endl;

  // STEP 1
  Integer notIncludedBCs = 0;
  Integer multipleBCs = 0;
  pdeNode2EQN_.Resize(numPDENodes_,dofsPerNode_);
  pdeNode2EQN_.Init(1);
 

  
//   std::cerr << "size of mesh2PDENode_" << mesh2PDENode_.GetSize() << std::endl;
//   std::cerr << "size of pdeNode2EQN_" << pdeNode2EQN_.GetSizeRow() << std::endl;
  
  // STEP 2

  Matrix<Integer> countNodes;
  countNodes.Resize(numPDENodes_,dofsPerNode_);
  countNodes.Init(0);
  for (Integer i=0; i<homoDirichletNodes_.GetSize(); i++)
    {
//       std::cerr << "homoDirichletNodes_[i]-1 = " << homoDirichletNodes_[i]-1 << std::endl;
//       std::cerr << "mesh2PDENode_[homoDirichletNodes_[i]-1]-1 = " << mesh2PDENode_[homoDirichletNodes_[i]-1]-1 << std::endl;
      // Check if homDirichletNode belongs to one of my 
      // subdomains
      if (mesh2PDENode_[homoDirichletNodes_[i]-1]-1 < 0)
	{
	  warnMsg  = "ScalarBlockEQN::CalcMapping: Homogen. Dirichlet node nr. ";
	  warnMsg += Info->GenStr(homoDirichletNodes_[i]);
	  warnMsg += " is not contained in any of the regions for this PDE";
	  Info->Warning(warnMsg, __FILE__, __LINE__);
	  notIncludedBCs++;
	} 
      else if (countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
	       [homoDirichletDofs_[i]-1] != 0)
	{
	  multipleBCs++;
	}
      else {
	pdeNode2EQN_[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
	  [homoDirichletDofs_[i]-1] = 0;
	countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
	  [homoDirichletDofs_[i]-1]++;
      }
    }


  if (multipleBCs > 0) {
    errMsg  = "ScalarBlockEQN::CalcMapping: Some hom dirichlet nodes ";
    errMsg += " occured already at least two times in the list of homDirichlet boundary nodes ";
    errMsg += "for this PDE!\n Please check, if this node is defined in";
    errMsg += " more than one level of boundary nodes!";
    Warning(errMsg.c_str(), __FILE__, __LINE__);
  }
  
  
  // STEP 3
  for (Integer i=0; i<constraintSlaveNodes_.GetSize(); i++)
    pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]
      [constraintDofs_[i]-1] = 0;
  
//    std::cerr << "after step3" << std::endl;
    // STEP 4
  for (Integer iNode=0; iNode<pde2MeshNode_.GetSize(); iNode++)
    for (Integer iDof=0; iDof<dofsPerNode_; iDof++)
      if (pdeNode2EQN_[iNode][iDof] != 0)
	{
	  eqnCounter++;
// 	  std::cerr << "size of pdeNode2EQN = " << pdeNode2EQN_.GetSizeRow() << std::endl;
// 	  std::cerr << "pdeNode2EQN_[" << iNode << "][" << iDof << "] = " << eqnCounter << std::endl;
	  pdeNode2EQN_[iNode][iDof] = eqnCounter;

	}
//     std::cerr << "after step4" << std::endl;
      
  // STEP 5
  for (Integer i=0; i<constraintSlaveNodes_.GetSize(); i++)
    pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]
      [constraintDofs_[i]-1] =
      pdeNode2EQN_[mesh2PDENode_[constraintMasterNodes_[i]-1]-1]
      [constraintDofs_[i]-1];

  // Now object is initialized
  isInitialized_ = TRUE;
  numEqns_ = eqnCounter;

  numBuildInDirichletEQNs_ = numPDENodes_ * dofsPerNode_ 
                           - numEqns_ + multipleBCs;

}

void ScalarBlockEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "ScalarBlockEQN::Print" );
  out << "Equation numbering - Information" << std::endl;
  out << "================================" << std::endl;
  out << "DOFs per Node: " << dofsPerNode_ << std::endl;
  out << "Using SCALAR-BLOCK numbering of equations" << std::endl;
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
      out << std::setw(10) << pdeNode2EQN_[iNode][0] << std::endl;
      
      // Print out further dofs if needed
      for (Integer iDof = 1; iDof < dofsPerNode_; iDof++)
	{
	  out << std::setw(10) << " " << " | ";
	  out << std::setw(13) << " " << " | ";
	  out << std::setw(10) << pdeNode2EQN_[iNode][iDof] << std::endl;
	}
    }
  
}


void ScalarBlockEQN::Node2EQN(const Integer nodeNr, 
			      const Integer dof,
			      Integer & eqnNr,
			      Integer & eqnDof) const
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );
 #ifdef CHECK_INDEX
  if (nodeNr > mesh2PDENode_.GetSize())
    Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
	  __FILE__, __LINE__);
#endif
  eqnNr = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1][dof-1];
  eqnDof = 1;
}

void ScalarBlockEQN::Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );

#ifdef CHECK_INDEX
  if (nodeNr > mesh2PDENode_.GetSize())
    Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
	  __FILE__, __LINE__);
#endif

  eqns.Resize(dofsPerNode_);
  
  for (Integer i=0; i<dofsPerNode_; i++)
    eqns[i] = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1][i];
}


void ScalarBlockEQN::Node2EQN(const StdVector<Integer> &nodeNr,
			      StdVector<Integer> &eqnNr) const
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );

  eqnNr.Resize(nodeNr.GetSize()*dofsPerNode_);

  for(Integer iNode=0; iNode<nodeNr.GetSize(); iNode++)
    {
#ifdef CHECK_INDEX
      if (nodeNr[iNode] > mesh2PDENode_.GetSize())
	Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
	      __FILE__, __LINE__);
#endif
      for (Integer iDof=0; iDof<dofsPerNode_; iDof++)
	eqnNr[iNode*dofsPerNode_ + iDof] =
	  pdeNode2EQN_[mesh2PDENode_[nodeNr[iNode]-1]-1][iDof];
    }
}

  // ==========================================
  //   Equation mapping according to reordering
  // ==========================================
  void ScalarBlockEQN::ReorderMapping(Integer *order) {

    ENTER_FCN( "ScalarBlockEQN::ReorderMapping" );

    for ( Integer i = 0; i < pdeNode2EQN_.GetSizeRow(); i++ ) {
      for ( Integer j = 0; j < pdeNode2EQN_.GetSizeCol(); j++ ) {
	if ( pdeNode2EQN_[i][j] > 0 ) {
	  pdeNode2EQN_[i][j] = order[pdeNode2EQN_[i][j]-1];
	}
	else if(pdeNode2EQN_[i][j] < 0 ) {
	  //due to constraints
	  pdeNode2EQN_[i][j] = -order[-pdeNode2EQN_[i][j]-1];	
	}
      }
    }
    
  }

} // end of namespace
