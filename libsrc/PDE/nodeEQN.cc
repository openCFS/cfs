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


void NodeEQN::CalcLocalGlobalMapping(StdVector<Integer> & mesh2PDENode,
				     StdVector<Integer> & pde2MeshNode,
				     StdVector<Integer> & mesh2PDEElem,
				     StdVector<Integer> & pde2MeshElem)
{
  ENTER_FCN( "NodeEQN::CalcLocalGlobalMapping" );
  
  mesh2PDENode.Resize(ptGrid_->GetMaxnumnodes(actlevel_));
  mesh2PDENode.Init(-1);
  pde2MeshNode.Clear();
  //std::cerr << "Number of global elems " << ptGrid_->GetMaxnumElem(actlevel_);
  mesh2PDEElem.Resize(ptGrid_->GetMaxnumElem(actlevel_));
  //std::cerr << "number of local elems " <<ptGrid_->GetMaxnumElem(actlevel_,subdoms_); 
  mesh2PDEElem.Init(-1);
  pde2MeshElem.Resize(ptGrid_->GetMaxnumElem(actlevel_,subdoms_));
  //std::cerr << "Size of pde2MeshElem" << pde2MeshElem.GetSize() << std::endl;
  pde2MeshElem.Init(-1);
  //std::cerr << "After init of pde2MeshElem" << std::endl;
  //std::cerr << "Size of pde2MeshEl

  Integer nodeCounter = 0;
  Integer elemCounter = 1;
 
  StdVector<Elem*> subdom;
 
  // iterate over all subdomains
  for (Integer iSD=0; iSD<subdoms_.GetSize(); iSD++)
    {
      ptGrid_->GetElemSD(subdom,subdoms_[iSD],actlevel_); 

      // iterate over all elems in subdomain
      for (Integer iElem=0; iElem<subdom.GetSize(); iElem++)
	{
	  // *** Mapping of Elements ***
	  mesh2PDEElem[subdom[iElem]->elemNum - 1 ] = elemCounter;
	  pde2MeshElem[elemCounter-1] = subdom[iElem]->elemNum;
	  elemCounter++;

	  
	  // *** Mapping of Nodes ***
	  
	  // iterate over all nodes in elem
	  for (Integer iNode=0; iNode<subdom[iElem]->connect.GetSize(); iNode++)
	    // Check if node was already assigned
	    if (mesh2PDENode[subdom[iElem]->connect[iNode]-1] == -1)
	      {
		mesh2PDENode[subdom[iElem]->connect[iNode]-1] = ++nodeCounter;
		pde2MeshNode.Push_back(subdom[iElem]->connect[iNode]);
	      }
	}
    }

  numPDENodes_ = nodeCounter;
}



} // end of namespace
