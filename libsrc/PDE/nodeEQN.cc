#include "nodeEQN.hh"

#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField
{
  
NodeEQN::NodeEQN(Grid * aptgrid, 
		 std::vector<std::string>& asubdoms, 
		 Integer actlevel, 
		 Integer dofsPerNode)
  : BaseEQN(aptgrid, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "NodeEQN::NodeEQN" );
  Info->Error( "Not implemented" );
}
  

NodeEQN::~NodeEQN()
{
 
  ENTER_FCN( "NodeEQN::NodeEQN" );
  Info->Error( "Not implemented" );
}
  

void NodeEQN::CalcMesh2PDENode(std::vector<Integer> & mesh2PDENode)
{
  ENTER_FCN( "NodeEQN::CalcMesh2PDENode" );
  
  mesh2PDENode.resize(ptgrid_->GetMaxnumnodes(actlevel_));
  
  Integer nodeCounter = 0;
 
  std::vector<Elem*> subdom;
 
  // iterate over all subdomains
  for (Integer iSD=0; iSD<subdoms_.size(); iSD++)
    {
      ptgrid_->GetElemSD(subdom,subdoms_[iSD],actlevel_); 

      // iterate over all elems in subdomain
      for (Integer iElem=0; iElem<subdom.size(); iElem++)
	// iterate over all nodes in elem
	for (Integer iNode=0; iNode<subdom[iElem]->connect.GetSize(); iNode++)
	  // Check if node was already assigned
	  if (mesh2PDENode[subdom[iElem]->connect[iNode]-1] = 0)
	    mesh2PDENode[subdom[iElem]->connect[iNode]-1] = ++nodeCounter;
    }
  numPDENodes_ = nodeCounter;
}



} // end of namespace
