#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include <Elements/1D/line.hh>

namespace CoupledField
{

Line :: Line()
{
#ifdef TRACE
  (*trace) << "entering Line::Line" << std::endl;
#endif

  std::string integtype;
  conf->get("line",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

  IsSet=FALSE;

	Init();
}

Line :: ~Line()
{
#ifdef TRACE
  (*trace) << "entering Line::~Line" << std::endl;
#endif

  if (IntWeights) delete IntWeights;
}

void Line :: Init()
{
#ifdef TRACE
  (*trace) << "entering Line::Init" << std::endl;
#endif

  Dim      = 1;
  NumNodes = 2;
  NumEdges = 1;
  NumFaces = 1;

  SetIntPoints();
  SetTransformFncAtIntPoints();
  SetDerTransformFncAtIntPoints();
  IsSet=TRUE;
}
void Line :: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Line::SetIntPoints" << std::endl;
#endif

  switch(IntegType)
    {
    case GaussOrder2:

      NumIntPoints=2;
      DegreeInteg=2;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=0;

      IntPoints[0][0] = -0.57735026919;
      IntPoints[1][0] =  0.57735026919;
      //IntPoints[2][0] =  0.57735026919;
      //IntPoints[3][0] = -0.57735026919;
      //IntPoints[0][1] = 1.0;
      //IntPoints[1][1] = -1.0;
      //IntPoints[2][1] =  0.57735026919;
      //IntPoints[3][1] =  0.57735026919;
      break;

    default:
      std::cerr << "Integration type " << IntegType
	   << " is not implemented\n" << std::endl; exit(-1);
    }
}

void Line :: SetTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Line::SetShapeFncAtIntPoints" << std::endl;
#endif
  Integer i;
  TransFnc1AtIP.Resize(NumIntPoints);
  TransFnc2AtIP.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      TransFnc1AtIP[i]=TransFnc1(IntPoints[i][0]);
      TransFnc2AtIP[i]=TransFnc2(IntPoints[i][0]);
    }
}

void Line :: SetDerTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Line::SetDerTransFncAtIntPoints" << std::endl;
#endif
  Integer i;

  DxiTransFnc1AtIP.Resize(NumIntPoints);
  DxiTransFnc2AtIP.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      DxiTransFnc1AtIP[i]=TransFnc1dxi();
      DxiTransFnc2AtIP[i]=TransFnc2dxi();
    }

}

Vector<Double> & Line::GetShFncAtIP(const Integer iShFnc)
{
#ifdef TRACE
  (*trace) << "entering Line::GetShFncAtIP" << std::endl;
#endif
  switch(iShFnc)
    {
    case 1:
      return TransFnc1AtIP;
      break;
    case 2:
      return TransFnc2AtIP;
      break;
    default:
      Error("Number of shape function is wrong", __FILE__, __LINE__);
    }
}

Double Line::getIntVal(const Point<2> * const nodes)
{
	return dist(nodes[0],nodes[1])*0.5;
}

} // end of namespace
