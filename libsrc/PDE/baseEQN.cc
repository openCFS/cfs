#include "baseEQN.hh"
#include <list>
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField
{

  BaseEQN :: BaseEQN(Grid * aptGrid, 
		     BCs * aptBCs,
		     std::vector<std::string>& asubdoms, 
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



void BaseEQN::SetHomoDirichletBCs(const std::vector<std::string> &nodeLevel,
				  const std::vector<std::string> &dofs)
{
  ENTER_FCN( "BaseEQN::SetHomoDirichletBCs" );

  std::list<Integer> tempNodeList;
  homoDirichletNodes_.clear();
  Integer oldSize = 0;
  Integer counter = 0;

  //std::cerr << "nodeLevel.size=" << nodeLevel.size() << std::endl;

  for (Integer i=0; i<nodeLevel.size(); i++)
    {
      tempNodeList = ptBCs_->GetNodesLevel(nodeLevel[i]);
      
      oldSize = homoDirichletNodes_.size();
      homoDirichletNodes_.resize(oldSize + tempNodeList.size());
      //std::cerr << "oldsize: " << oldSize << "tempNodeList.size= " << tempNodeList.size() << std::endl;

      std::list<Integer>::iterator it;
      
      for (it=tempNodeList.begin(); it != tempNodeList.end(); it++)
	homoDirichletNodes_[counter++] = (*it);
      
      if (dofsPerNode_ > 1)
	{
	  homoDirichletDofs_.insert(homoDirichletDofs_.end(),
				    GetBCDof(dofs[i]),
				    tempNodeList.size());
	}
	
    }
  
  // Debug
  //for (Integer i=0; i<homoDirichletNodes_.size(); i++)
    //   std::cerr << homoDirichletNodes_[i] << std::endl;

}


void BaseEQN::SetConstraints(const std::vector<Integer> &slaveNodeNrs,
			     const std::vector<Integer> &masterNodeNrs,
			     const std::vector<std::string> &dofs)
{
  ENTER_FCN( "BaseEQN::SetConstraints" );

  constraintSlaveNodes_ = slaveNodeNrs;
  constraintMasterNodes_ = masterNodeNrs;

  if (dofsPerNode_ > 1)
    {
      constraintDofs_.resize(dofs.size());
      
      for (Integer i=0; i<dofs.size(); i++)
	constraintDofs_[i] = GetBCDof(dofs[i]);
    }
}




Integer BaseEQN::GetBCDof(const std::string dofString) const
{
  ENTER_FCN( "BaseEQN::GetBCDof" );
  
  Integer retVal = 0;
  
  if (dofString == "ux") retVal = 1;
  else if (dofString == "uy") retVal = 2;
  else if (dofString == "uz") retVal = 3;
  else
    {
#ifndef XMLPARAMS
      std::string errmsg = "The direction '" + dofString;
      errmsg += "' mentioned in the config-file is not implemented";
      Info->Error( errmsg, __FILE__, __LINE__ );
#else
      // According to the Schema definition of the parameter file this cannot
      // happen. Did the parser not perform validation?
      std::string errmsg = "Direction should be one of ux, uy, uz\n";
      errmsg += "and not" + dofString;
      Info->Error( errmsg, __FILE__, __LINE__ );
#endif
    }
  
  return retVal;
}

} // end of namespace
