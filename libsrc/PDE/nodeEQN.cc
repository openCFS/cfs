#include "nodeEQN.hh"

#include "Domain/grid.hh"

namespace CoupledField
{
  
NodeEQN::NodeEQN(Grid * aptGrid, 
		 BCs * aptBCs,
		 StdVector<std::string>& asubdoms, 
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
  
void NodeEQN::Mesh2PDENode(StdVector<Integer> & PDENodes,
			   const StdVector<Integer> & MeshNodes) const
{
  ENTER_FCN( "NodeEQN::Mesh2PDENode" );
  
  PDENodes.Resize(MeshNodes.GetSize());
   
  for (Integer i=0; i<MeshNodes.GetSize(); i++) 
    PDENodes[i] = mesh2PDENode_[MeshNodes[i]-1];
  
}


void NodeEQN::CalcMesh2PDENode(StdVector<Integer> & mesh2PDENode,
			       StdVector<Integer> & pde2MeshNode)
{
  ENTER_FCN( "NodeEQN::CalcMesh2PDENode" );
  
  mesh2PDENode.Clear();
  mesh2PDENode.Resize(ptGrid_->GetMaxnumnodes(actlevel_));
  pde2MeshNode.Clear();
  
  Integer nodeCounter = 0;
 
  StdVector<Elem*> subdom;
 
  // iterate over all subdomains
  for (Integer iSD=0; iSD<subdoms_.GetSize(); iSD++)
    {
      ptGrid_->GetElemSD(subdom,subdoms_[iSD],actlevel_); 

      // iterate over all elems in subdomain
      for (Integer iElem=0; iElem<subdom.GetSize(); iElem++)
	// iterate over all nodes in elem
	for (Integer iNode=0; iNode<subdom[iElem]->connect.GetSize(); iNode++)
	  // Check if node was already assigned
	  if (mesh2PDENode[subdom[iElem]->connect[iNode]-1] == 0)
	    {
	      mesh2PDENode[subdom[iElem]->connect[iNode]-1] = ++nodeCounter;
	      pde2MeshNode.Push_back(subdom[iElem]->connect[iNode]);
	    }
    }

  numPDENodes_ = nodeCounter;
}



} // end of namespace
