#include "blocknodeEQN.hh"
#include <iomanip>

namespace CoupledField
{
  
  BlockNodeEQN::BlockNodeEQN(Grid * aptGrid, 
                             BCs * aptBCs,
                             StdVector<std::string>& asubdoms, 
                             Integer actlevel, 
                             Integer dofsPerNode)
    : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
  {
    ENTER_FCN( "BlockNodeEQN::BlockNodeEQN" );

    isBlockMapped_ = TRUE;
    dofsPerEQN_ = dofsPerNode;
  }

  BlockNodeEQN::~BlockNodeEQN()
  {
 
    ENTER_FCN( "BlockNodeEQN::BlockNodeEQN" );
  }


  void BlockNodeEQN::CalcMapping()
  {

    ENTER_FCN( "BlockNodeEQN::CalcMapping" );
 
    // First apply Mapping from global to
    // local node/elem numbers and back
    CalcLocalGlobalMapping(mesh2PDENode_,
                           pde2MeshNode_,
                           mesh2PDEElem_,
                           pde2MeshElem_);
 
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
    std::string warnMsg, errMsg;
    pdeNode2EQN_.Clear();
    pdeNode2EQN_.Resize(numPDENodes_);

    // STEP 2
    // Check if there exist nodes, which only have
    // hom. Dirichlet BC dof
    StdVector<Integer> numDirichletDofsPerNode;
    numDirichletDofsPerNode.Resize(numPDENodes_);
    for (Integer i=0; i<homoDirichletNodes_.GetSize(); i++)
      {
        if (mesh2PDENode_[homoDirichletNodes_[i]-1]-1 < 0)
          {
            warnMsg  = "BlockNodeEQN::CalcMapping: Homogen. Dirichlet node nr. ";
            warnMsg += Info->GenStr(homoDirichletNodes_[i]);
            warnMsg += " is not contained in any of the regions for this PDE";
            Info->Warning(warnMsg, __FILE__, __LINE__);
          }
        else if (numDirichletDofsPerNode[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
                 > dofsPerNode_) {
          errMsg  = "BlockNodeEQN::CalcMapping: Homogen. Dirichlet node nr. ";
          errMsg += Info->GenStr(homoDirichletNodes_[i]);
          errMsg += " is not contained in any of the regions for this PDE";
          Warning(errMsg.c_str(), __FILE__, __LINE__);
        
        }
        else
          numDirichletDofsPerNode[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]++;
      }

    // REMOVE LATER !!!!!!!!!!!!!!!!!!!
    //numDirichletDofsPerNode.Init(0);

    // STEP 3
    // Check if there are constraint slavenodes, where
    // all dofs depend on the same dofs of the same
    // master node
    StdVector<Integer> numConstraintDofsPerNode;
    StdVector<Integer> masterNodes;
    numConstraintDofsPerNode.Resize(numPDENodes_);
    masterNodes.Resize(numPDENodes_);

    for (Integer i=0; i<constraintSlaveNodes_.GetSize(); i++)
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
                               
    //std::cerr << "********************************" << std::endl;
    //std::cerr << "** CALCULATING MAPPING *********" << std::endl;
    //std::cerr << std::endl;
    //std::cerr << "pde2MeshNode" << std::endl << "----------" << std::endl;
    //std::cerr << pde2MeshNode_ << std::endl;

    // STEP 4
    for (Integer i=0; i<pde2MeshNode_.GetSize(); i++)
      if (numDirichletDofsPerNode[i] != dofsPerNode_ &&
          numConstraintDofsPerNode[i] != dofsPerNode_)
        {
          eqnCounter ++;
          pdeNode2EQN_[i] = eqnCounter;
          //std::cerr << "Pushing back" << (pde2MeshNode_[eqnCounter-1]-1)*dofsPerNode_ << std::endl;
        } 
  
  
    // Now count number of dirichlet BCs, which were not 
    // thrown out
    numBuildInDirichletEQNs_ = 0;
    for (Integer i=0; i<numDirichletDofsPerNode.GetSize(); i++)
      if (numDirichletDofsPerNode[i] == dofsPerNode_)
        numBuildInDirichletEQNs_ += dofsPerNode_;

   
    //std::cerr << "NumBuildInDirichletEQNs = " <<numBuildInDirichletEQNs_ << std::endl; 
  

    // Now object is initialized
    numDirichletDofsPerNode.Clear();
    numConstraintDofsPerNode.Clear();
    masterNodes.Clear();
    isInitialized_ = TRUE;
    numEqns_ = eqnCounter;
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
  

    for (Integer i=0; i<pde2MeshNode_.GetSize(); i++)
      {
        out << std::setw(10) << i+1  << " | ";
        out << std::setw(13) << pde2MeshNode_[i] << " | ";
        out << std::setw(10) << pdeNode2EQN_[i] << std::endl;
      }
  }


  void BlockNodeEQN::Node2EQN(const Integer nodeNr, 
                              const Integer dof,
                              Integer & eqnNr,
                              Integer & eqnDof) const 
  {
    ENTER_FCN( "BlockNodeEQN::Node2EQN" );
#ifdef CHECK_INDEX
    if (nodeNr > mesh2PDENode_.GetSize())
      Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
            __FILE__, __LINE__);
#endif
    eqnNr = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1];
    eqnDof = dof;
  }

  void BlockNodeEQN::Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const
  {
    ENTER_FCN( "BlockNodeEQN::Node2EQN" );
    Error( "Not implemented",__FILE__,__LINE__);
  }
  

  void BlockNodeEQN::Node2EQN(const StdVector<Integer> &nodeNr,
                              StdVector<Integer> &eqnNr) const
  {
    ENTER_FCN( "BlockNodeEQN::Node2EQN" );

    eqnNr.Resize(nodeNr.GetSize());

    for (Integer i=0; i<nodeNr.GetSize(); i++)
      eqnNr[i] =  pdeNode2EQN_[mesh2PDENode_[nodeNr[i]-1]-1];
  }

  // ==========================================
  //   Equation mapping according to reordering
  // ==========================================
  void BlockNodeEQN::ReorderMapping(Integer *order) {

    ENTER_FCN( "BlockNodeEQN::ReorderMapping" );

    for ( Integer i = 0; i < pdeNode2EQN_.GetSize(); i++ ) {
      if ( pdeNode2EQN_[i] > 0 ) {
        pdeNode2EQN_[i] = order[pdeNode2EQN_[i]-1];
      }
      else if(pdeNode2EQN_[i] < 0 ) {
        //due to constraints
        pdeNode2EQN_[i] = -order[-pdeNode2EQN_[i]-1];   
      }
    }
  }


} // end of namespace
