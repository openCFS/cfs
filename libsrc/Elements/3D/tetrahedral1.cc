#include <iostream>
#include <fstream>

#include "tetrahedral1.hh"

namespace CoupledField
{

Tetrahedral1::Tetrahedral1(enum IntegrationType aintegtype):GeTetrahedral(aintegtype)
{ 
#ifdef TRACE
  (*trace) << "entering Tetrahedral1::Tetrahedral1" << std::endl;
#endif
   ElemType = TETRAHEDRAL1;
   Init();
}

Tetrahedral1::~Tetrahedral1()
{
#ifdef TRACE
  (*trace) << "entering Tetrahedral1::~Tetrahedral1" << std::endl;
#endif
 ; 
}

void Tetrahedral1::Init()
{
#ifdef TRACE
  (*trace) << "entering Tetrahedral1::Init" << std::endl;
#endif
  
  Dim = 3;
  NumNodes = 4;
  NumEdges = 6;
  NumFaces = 4;

  SetIntPoints();
  SetTransformFncAtIntPoints();
  SetDerTransformFncAtIntPoints();
  IsSet=TRUE;
}

} // end of namespace

