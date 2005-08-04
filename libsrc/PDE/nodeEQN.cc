#include "nodeEQN.hh"

#include "Domain/grid.hh"

namespace CoupledField
{
  
  // ***************
  //   Constructor
  // ***************
  NodeEQN::NodeEQN( Grid *aptGrid, StdVector<RegionIdType> &asubdoms, 
                    UInt dofsPerNode, Boolean sortEQNs )
    : BaseEQN( aptGrid, asubdoms, dofsPerNode, sortEQNs ) {

    ENTER_FCN( "NodeEQN::NodeEQN" );

    isNodalMapped_ = TRUE;
  }
  

  // **************
  //   Destructor
  // **************
  NodeEQN::~NodeEQN() {
 
    ENTER_FCN( "NodeEQN::NodeEQN" );

  }
  
  void NodeEQN::Mesh2PDENode(StdVector<UInt> & PDENodes,
                             const StdVector<UInt> & MeshNodes) const
  {
    ENTER_FCN( "NodeEQN::Mesh2PDENode" );
  
    PDENodes.Resize(MeshNodes.GetSize());
   
    for (UInt i=0; i<MeshNodes.GetSize(); i++) 
      PDENodes[i] = mesh2PDENode_[MeshNodes[i]-1];
  
  }

  void NodeEQN::PDE2MeshNode(StdVector<UInt> & meshNodes,
                             const StdVector<UInt> & pdeNodes) const
  {
    ENTER_FCN( "NodeEQN::PDE2MeshNode" );
  
    meshNodes.Resize(pdeNodes.GetSize());
   
    for (UInt i=0; i<pdeNodes.GetSize(); i++) 
      meshNodes[i] = pde2MeshNode_[pdeNodes[i]-1];
  
  }


  void NodeEQN::CalcLocalGlobalMapping(StdVector<Integer> & mesh2PDENode,
                                       StdVector<UInt> & pde2MeshNode,
                                       StdVector<Integer> & mesh2PDEElem,
                                       StdVector<UInt> & pde2MeshElem)
  {
    ENTER_FCN( "NodeEQN::CalcLocalGlobalMapping" );
  
    mesh2PDENode.Resize(ptGrid_->GetNumNodes());
    mesh2PDENode.Init(-1);
    pde2MeshNode.Clear();

    mesh2PDEElem.Resize(ptGrid_->GetNumElems());
    mesh2PDEElem.Init(-1);
    pde2MeshElem.Resize(ptGrid_->GetNumElems(subdoms_));
    pde2MeshElem.Init(0);
    //std::cerr << "After init of pde2MeshElem" << std::endl;
    //std::cerr << "Size of pde2MeshEl


    UInt nodeCounter = 0;
    UInt elemCounter = 1;
 
    StdVector<Elem*> subdom;
 

    // iterate over all subdomains
    for (UInt iSD=0; iSD<subdoms_.GetSize(); iSD++)
      {
        ptGrid_->GetElems(subdom,subdoms_[iSD]); 

        // iterate over all elems in subdomain
        for (UInt iElem=0; iElem<subdom.GetSize(); iElem++)
          {
            //std::cerr << "before mapping of elements" << std::endl;
            // *** Mapping of Elements ***
            mesh2PDEElem[subdom[iElem]->elemNum - 1 ] = elemCounter;
            pde2MeshElem[elemCounter-1] = subdom[iElem]->elemNum;
            elemCounter++;
            //std::cerr << "after mapping of elements" << std::endl;

          
            // *** Mapping of Nodes ***
            // iterate over all nodes in elem
            for (UInt iNode=0; iNode<subdom[iElem]->connect.GetSize(); iNode++)
              // Check if node was already assigned
              if (mesh2PDENode[subdom[iElem]->connect[iNode]-1] == -1)
                {
                  // std::cerr << "mesh2PDENode[" << subdom[iElem]->connect[iNode]-1 << "] = ";
                  //              std::cerr << nodeCounter +1 << std::endl;
                  mesh2PDENode[subdom[iElem]->connect[iNode]-1] = ++nodeCounter;
                  pde2MeshNode.Push_back(subdom[iElem]->connect[iNode]);
                }
          }
      }

    numPDENodes_ = nodeCounter;
  }



} // end of namespace
