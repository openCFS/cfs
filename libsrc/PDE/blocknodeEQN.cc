#include "blocknodeEQN.hh"

namespace CoupledField
{
  
BlockNodeEQN::BlockNodeEQN(Grid * aptgrid, 
		 std::vector<std::string>& asubdoms, 
		 Integer actlevel, 
		 Integer dofsPerNode)
  : NodeEQN(aptgrid, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "BlockNodeEQN::BlockNodeEQN" );
  Info->Error( "Not implemented" );
}

BlockNodeEQN::~BlockNodeEQN()
{
 
  ENTER_FCN( "BlockNodeEQN::BlockNodeEQN" );
  Info->Error( "Not implemented" );
}


void BlockNodeEQN::CalcMapping()
{

  ENTER_FCN( "BlockNodeEQN::CalcMapping" );
  Info->Error( "Not implemented" );
}

void BlockNodeEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "BlockNodeEQN::Print" );
  Info->Error( "Not implemented" );
}


void BlockNodeEQN::EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
				     std::vector<Integer> &pos)
{
  ENTER_FCN( "BlockNodeEQN::EQN2SolVectorPos" );
  Info->Error( "Not implemented" );
}


void BlockNodeEQN::Node2EQN(const Integer nodeNr, std::vector<Integer> &eqns)
{
  ENTER_FCN( "BlockNodeEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}
  

void BlockNodeEQN::Node2EQN(const std::vector<Integer> &nodeNr,
			     std::vector<Integer> &eqnNr)
{
  ENTER_FCN( "BlockNodeEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}

} // end of namespace
