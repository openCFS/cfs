#include "scalarnodeEQN.hh"

namespace CoupledField
{
  
ScalarNodeEQN::ScalarNodeEQN(Grid * aptgrid, 
		 std::vector<std::string>& asubdoms, 
		 Integer actlevel, 
		 Integer dofsPerNode)
  : NodeEQN(aptgrid, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "NodeEQN::NodeEQN" );
  Info->Error( "Not implemented" );
}

ScalarNodeEQN::~ScalarNodeEQN()
{
 
  ENTER_FCN( "NodeEQN::NodeEQN" );
  Info->Error( "Not implemented" );
}


void ScalarNodeEQN::CalcMapping()
{

  ENTER_FCN( "ScalarNodeEQN::CalcMapping" );
  Info->Error( "Not implemented" );
}

void ScalarNodeEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "ScalarNodeEQN::Print" );
  Info->Error( "Not implemented" );
}


void ScalarNodeEQN::EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
				     std::vector<Integer> &pos)
{
  ENTER_FCN( "ScalarNodeEQN::EQN2SolVectorPos" );
  Info->Error( "Not implemented" );
}


void ScalarNodeEQN::Node2EQN(const Integer nodeNr, std::vector<Integer> &eqns)
{
  ENTER_FCN( "ScalarNodeEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}
  

void ScalarNodeEQN::Node2EQN(const std::vector<Integer> &nodeNr,
			     std::vector<Integer> &eqnNr)
{
  ENTER_FCN( "ScalarNodeEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}

} // end of namespace
