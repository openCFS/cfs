#include "scalarblockEQN.hh"

namespace CoupledField
{
  
ScalarBlockEQN::ScalarBlockEQN(Grid * aptGrid, 
			       BCs * aptBCs,
			       StdVector<std::string>& asubdoms, 
			       Integer actlevel, 
			       Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
  Error( "Not implemented" );
}

ScalarBlockEQN::~ScalarBlockEQN()
{
 
  ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
  Error( "Not implemented" );
}


void ScalarBlockEQN::CalcMapping()
{

  ENTER_FCN( "ScalarBlockEQN::CalcMapping" );
  Error( "Not implemented" );
}

void ScalarBlockEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "ScalarBlockEQN::Print" );
  Error( "Not implemented" );
}


void ScalarBlockEQN::EQN2SolVectorPos(const StdVector<Integer> &eqnNr, 
				     StdVector<Integer> &pos) const
{
  ENTER_FCN( "ScalarBlockEQN::EQN2SolVectorPos" );
  Error( "Not implemented" );
}


Integer ScalarBlockEQN::Node2EQN(const Integer nodeNr, 
			      const Integer dof) const
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );
  Error( "Not implemented yet" );
  return 0;
}

void ScalarBlockEQN::Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );
  Error( "Not implemented" );
}
  

void ScalarBlockEQN::Node2EQN(const StdVector<Integer> &nodeNr,
			     StdVector<Integer> &eqnNr) const
{
  ENTER_FCN( "ScalarBlockEQN::Node2EQN" );
  Error( "Not implemented" );
}

} // end of namespace
