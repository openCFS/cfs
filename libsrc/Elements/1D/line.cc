#include "line.hh"

namespace CoupledField
{

Double Line::getIntVal(Point<2> * nodes)
{
	return dist(nodes[0],nodes[1])*0.5;
}

} // end of namespace
