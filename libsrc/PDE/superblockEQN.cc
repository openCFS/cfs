#include "superblockEQN.hh"

namespace CoupledField
{
  
SuperBlockEQN::SuperBlockEQN(Grid * aptGrid, 
			     BCs * aptBCs,
			     StdVector<std::string>& asubdoms, 
			     Integer actlevel, 
			     Integer dofsPerNode)
  : NodeEQN(aptGrid, aptBCs, asubdoms, actlevel, dofsPerNode)
{
  ENTER_FCN( "SuperBlockEQN::SuperBlockEQN" );
  Error( "Not implemented" );
}

SuperBlockEQN::~SuperBlockEQN()
{
 
  ENTER_FCN( "SuperBlockEQN::SuperBlockEQN" );
  Error( "Not implemented" );
}


void SuperBlockEQN::CalcMapping()
{

  ENTER_FCN( "SuperBlockEQN::CalcMapping" );
  Error( "Not implemented" );
}

void SuperBlockEQN::Print(std::ostream & out) const
{
  ENTER_FCN( "SuperBlockEQN::Print" );
  Error( "Not implemented" );
}


void SuperBlockEQN::EQN2SolVectorPos(const StdVector<Integer> &eqnNr, 
				     StdVector<Integer> &pos) const
{
  ENTER_FCN( "SuperBlockEQN::EQN2SolVectorPos" );
  Error( "Not implemented" );
}

Integer SuperBlockEQN::Node2EQN(const Integer nodeNr, 
				const Integer dof) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );
  Error( "Not implemented yet" );
  return 0;
}

void SuperBlockEQN::Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );
  Error( "Not implemented" );
}
  

void SuperBlockEQN::Node2EQN(const StdVector<Integer> &nodeNr,
			     StdVector<Integer> &eqnNr) const
{
  ENTER_FCN( "SuperBlockEQN::Node2EQN" );
  Error( "Not implemented" );
}

} // end of namespace
