#include "baseEQN.hh"
#include <list>
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField
{

  BaseEQN :: BaseEQN(Grid * aptGrid, 
		     BCs * aptBCs,
		     StdVector<std::string>& asubdoms, 
		     Integer actlevel, 
		     Integer dofsPerNode)
{
  ENTER_FCN( "BaseEQN::BaseEQN" );

  ptGrid_   = aptGrid;
  ptBCs_    = aptBCs;
  subdoms_  = asubdoms;
  actlevel_ = actlevel;
  dofsPerNode_ = dofsPerNode;

  isInitialized_ = FALSE;
}


BaseEQN :: ~BaseEQN()
{
  ENTER_FCN( "~BaseEQN::~BaseEQN" );

}



void BaseEQN::SetHomoDirichletBCs(const StdVector<std::string> &nodeLevel,
				  const StdVector<std::string> &dofs)
{
  ENTER_FCN( "BaseEQN::SetHomoDirichletBCs" );

  std::list<Integer> tempNodeList;
  homoDirichletNodes_.Clear();

  for (Integer i=0; i<nodeLevel.GetSize(); i++)
    {
      
      tempNodeList = ptBCs_->GetNodesLevel(nodeLevel[i]);
      
      std::list<Integer>::iterator it;
      
      for (it=tempNodeList.begin(); it != tempNodeList.end(); it++)
	{
	  homoDirichletNodes_.Push_back(*it);
	  if (dofsPerNode_ > 1)
	    homoDirichletDofs_.Push_back(GetBCDof(dofs[i]));
	}
      
    }
  
}

void BaseEQN::SetConstraints(const StdVector<Integer> &slaveNodeNrs,
			     const StdVector<Integer> &masterNodeNrs,
			     const StdVector<std::string> &dofs)
{
  ENTER_FCN( "BaseEQN::SetConstraints" );

  constraintSlaveNodes_ = slaveNodeNrs;
  constraintMasterNodes_ = masterNodeNrs;

  if (dofsPerNode_ > 1)
    {
      constraintDofs_.Resize(dofs.GetSize());
      
      for (Integer i=0; i<dofs.GetSize(); i++)
	constraintDofs_[i] = GetBCDof(dofs[i]);
    }
}




Integer BaseEQN::GetBCDof(const std::string dofString) const
{
  ENTER_FCN( "BaseEQN::GetBCDof" );
  
  Integer retVal = 0;
  
  if (dofString == "ux") 
    retVal = 1;
  else if 
    (dofString == "uy") retVal = 2;
  else if 
    (dofString == "uz") retVal = 3;
  
  // hard-coded for Piezo PDE
  else if 
    (dofString == "ep") retVal = dofsPerNode_;

  else {
    // According to the Schema definition of the parameter file this cannot
    // happen. Did the parser not perform validation?
    std::string errmsg = "Direction should be one of ux, uy, uz\n";
    errmsg += "and not" + dofString;
    Info->Error( errmsg, __FILE__, __LINE__ );
  }
  
  return retVal;
}

} // end of namespace
