#include "baseEQN.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField
{

  BaseEQN :: BaseEQN(Grid * aptgrid, 
		     std::vector<std::string>& asubdoms, 
		     Integer actlevel, 
		     Integer dofsPerNode)
{
  ENTER_FCN( "BaseEQN::BaseEQN" );

  ptgrid_   = aptgrid;
  subdoms_  = asubdoms;
  actlevel_ = actlevel;
  dofsPerNode_ = dofsPerNode;

  isInitialized_ = FALSE;
}


BaseEQN :: ~BaseEQN()
{
  ENTER_FCN( "~BaseEQN::~BaseEQN" );

}



void BaseEQN::SetHomDirichletBCs(const std::vector<Integer> &nodeNrs,
			const std::vector<std::string> &dofs)
{
  ENTER_FCN( "BaseEQN::SetHomDirichletBCs" );

  dirichletNodes_ = nodeNrs;

  if (dofsPerNode_ > 1)
    {
      dirichletDofs_.resize(dofs.size());
      
      for (Integer i=0; i<dofs.size(); i++)
	{
	  dirichletDofs_[i] = GetBCDof(dofs[i]);
	}
    }			   
}

void BaseEQN::SetConstraints(const std::vector<Integer> &nodeNrs,
			     const std::vector<std::string> &dofs)
{
  ENTER_FCN( "BaseEQN::SetConstraints" );
  Info->Error( "Not implemented yet" );
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
