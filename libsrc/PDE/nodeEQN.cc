#include "nodeEQN.hh"

#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField
{
  
NodeEQN::NodeEQN(Grid * aptGrid, 
		 BCs * aptBCs,
		 std::vector<std::string>& asubdoms, 
		 Integer actlevel, 
		 Integer dofsPerNode)
  : BaseEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "NodeEQN::NodeEQN" );

  isNodalMapped_ = TRUE;
}
  

NodeEQN::~NodeEQN()
{
 
  ENTER_FCN( "NodeEQN::NodeEQN" );
}
  

void NodeEQN::CalcMesh2PDENode(std::vector<Integer> & mesh2PDENode,
			       std::vector<Integer> & pde2MeshNode)
{
  ENTER_FCN( "NodeEQN::CalcMesh2PDENode" );
  
  mesh2PDENode.clear();
  mesh2PDENode.resize(ptGrid_->GetMaxnumnodes(actlevel_),0);
  pde2MeshNode.clear();
  
  Integer nodeCounter = 0;
 
  std::vector<Elem*> subdom;
 
  // iterate over all subdomains
  for (Integer iSD=0; iSD<subdoms_.size(); iSD++)
    {
      ptGrid_->GetElemSD(subdom,subdoms_[iSD],actlevel_); 

      // iterate over all elems in subdomain
      for (Integer iElem=0; iElem<subdom.size(); iElem++)
	// iterate over all nodes in elem
	for (Integer iNode=0; iNode<subdom[iElem]->connect.GetSize(); iNode++)
	  // Check if node was already assigned
	  if (mesh2PDENode[subdom[iElem]->connect[iNode]-1] == 0)
	    {
	      mesh2PDENode[subdom[iElem]->connect[iNode]-1] = ++nodeCounter;
	      pde2MeshNode.push_back(subdom[iElem]->connect[iNode]);
	    }
    }

  numPDENodes_ = nodeCounter;
}



} // end of namespace
