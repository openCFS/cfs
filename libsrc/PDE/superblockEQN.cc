#include "superblockEQN.hh"

namespace CoupledField
{
  
SuperBlockEQN::SuperBlockEQN(Grid * aptGrid, 
			     BCs * aptBCs,
			     std::vector<std::string>& asubdoms, 
			     Integer actlevel, 
			     Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "SuperBlockEQN::SuperBlockEQN" );
  Info->Error( "Not implemented" );
}

SuperBlockEQN::~SuperBlockEQN()
{
 
  ENTER_FCN( "SuperBlockEQN::SuperBlockEQN" );
  Info->Error( "Not implemented" );
}


void SuperBlockEQN::CalcMapping()
{

  ENTER_FCN( "SuperBlockEQN::CalcMapping" );
  Info->Error( "Not implemented" );
}

void SuperBlockEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "SuperBlockEQN::Print" );
  Info->Error( "Not implemented" );
}


void SuperBlockEQN::EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
				     std::vector<Integer> &pos) const
{
  ENTER_FCN( "SuperBlockEQN::EQN2SolVectorPos" );
  Info->Error( "Not implemented" );
}


void SuperBlockEQN::Node2EQN(const Integer nodeNr, std::vector<Integer> &eqns) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}
  

void SuperBlockEQN::Node2EQN(const std::vector<Integer> &nodeNr,
			     std::vector<Integer> &eqnNr) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}

} // end of namespace
