#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

#include <general_head.hh>
#include <utils_head.hh>
#include "baseelem.hh"
#include "triangle1.hh"

namespace CoupledField
{
                   
Triangle1 :: Triangle1(ShortInt aintegtype) : BaseElem()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::Triangle1" << endl;
#endif
  
  IntegType = aintegtype;
}
  
Triangle1 :: ~Triangle1()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::~Triangle1" << endl;
#endif

  ;
}

void Triangle1 :: Init()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::Init" << endl;
#endif

  NumNodes = 3;
  NumEdges = 3;
  NumFaces = 1;
  ;
}

void Triangle1 :: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::SetIntegration" << endl;
#endif

  ;
}

void Triangle1 :: SetShapeFnc()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::SetShapeFnc" << endl;
#endif

  ;
}

void Triangle1 :: SetDShapeFnc()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::SetDShapeFnc" << endl;
#endif

  ;
}

}
