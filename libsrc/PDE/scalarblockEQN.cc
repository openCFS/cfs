#include "scalarblockEQN.hh"

namespace CoupledField
{
  
ScalarBlockEQN::ScalarBlockEQN(Grid * aptgrid, 
		 std::vector<std::string>& asubdoms, 
		 Integer actlevel, 
		 Integer dofsPerNode)
  : NodeEQN(aptgrid, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
  Info->Error( "Not implemented" );
}

ScalarBlockEQN::~ScalarBlockEQN()
{
 
  ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
  Info->Error( "Not implemented" );
}


void ScalarBlockEQN::CalcMapping()
{

  ENTER_FCN( "ScalarBlockEQN::CalcMapping" );
  Info->Error( "Not implemented" );
}

void ScalarBlockEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "ScalarBlockEQN::Print" );
  Info->Error( "Not implemented" );
}


void ScalarBlockEQN::EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
				     std::vector<Integer> &pos)
{
  ENTER_FCN( "ScalarBlockEQN::EQN2SolVectorPos" );
  Info->Error( "Not implemented" );
}


void ScalarBlockEQN::Node2EQN(const Integer nodeNr, std::vector<Integer> &eqns)
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}
  

void ScalarBlockEQN::Node2EQN(const std::vector<Integer> &nodeNr,
			     std::vector<Integer> &eqnNr)
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );
  Info->Error( "Not implemented" );
}

} // end of namespace
