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
   ElemType  = TRIANGLE1;   
   IntegType = aintegtype;
   Init();
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

  Dim = 2;
  NumNodes = 3;
  NumEdges = 3;
  NumFaces = 1;

  SetIntPoints();
  SetShapeFncAtIntPoints();
  SetDShapeFncAtIntPoints();
  ;
}

void Triangle1 :: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::SetIntPoints" << endl;
#endif

 switch(IntegType)
    {
    case GaussOrder3:
      NumIntPoints=4;
      DegreeInteg=3;
      IntPoints=new  Matrix<Double>(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);
      (*IntPoints)[0][0] = 1.0/3.0;
      (*IntPoints)[1][0] = 3.0/5.0;
      (*IntPoints)[2][0] = 1.0/5.0;
      (*IntPoints)[3][0] = 1.0/5.0;
      (*IntPoints)[0][1] = 1.0/3.0;
      (*IntPoints)[1][1] = 1.0/5.0;
      (*IntPoints)[2][1] = 3.0/5.0;
      (*IntPoints)[3][1] = 1.0/5.0;
      
      IntWeights[0]=-0.5625;
      IntWeights[1]=0.520833333333333;
      IntWeights[2]=0.520833333333333;
      IntWeights[3]=0.520833333333333;

      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 4 nodes, degree of precision is 3 " << endl;
      break;

    case GaussOrder2:
      NumIntPoints=3;
      DegreeInteg=2;
      IntPoints=new  Matrix<Double>(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);
      (*IntPoints)[0][0] = 0.166666666666667;
      (*IntPoints)[1][0] = 0.666666666666667; 
      (*IntPoints)[2][0] = 0.166666666666667;
      (*IntPoints)[0][1] = 0.166666666666667;
      (*IntPoints)[1][1] = 0.166666666666667;
      (*IntPoints)[2][1] = 0.666666666666667;

      IntWeights[0]= 0.166666666666667 ;
      IntWeights[1]= 0.166666666666667 ;
      IntWeights[2]= 0.166666666666667 ;

      if (InfoPrint)
    (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 3 nodes, degree of precision is 2 " << endl;
      break;
 
    default:
      cerr << "Integration type " << IntegType
           << " is not implemented \n" << endl; exit(-1);
    }
  
}

void Triangle1 :: SetShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::SetShapeFncAtIntPoints" << endl;
#endif
 Integer i;
  ShFncAtIP1.Resize(NumIntPoints);
  ShFncAtIP2.Resize(NumIntPoints);
  ShFncAtIP3.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      ShFncAtIP1[i]=ShapeFnc1((*IntPoints)[i][0],(*IntPoints)[i][1]);
      ShFncAtIP2[i]=ShapeFnc2((*IntPoints)[i][0],(*IntPoints)[i][1]);
      ShFncAtIP3[i]=ShapeFnc3((*IntPoints)[i][0],(*IntPoints)[i][1]);
    }
  
}

void Triangle1 :: SetDShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::SetDShapeFnc" << endl;
#endif
 Integer i;

  DxShFncAtIP1.Resize(NumIntPoints);
  DxShFncAtIP2.Resize(NumIntPoints);
  DxShFncAtIP3.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      DxShFncAtIP1[i]=ShapeFnc1dx((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DxShFncAtIP2[i]=ShapeFnc2dx((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DxShFncAtIP3[i]=ShapeFnc3dx((*IntPoints)[i][0],(*IntPoints)[i][1]);
    }
 
  DyShFncAtIP1.Resize(NumIntPoints);
  DyShFncAtIP2.Resize(NumIntPoints);
  DyShFncAtIP3.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      DyShFncAtIP1[i]=ShapeFnc1dy((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DyShFncAtIP2[i]=ShapeFnc2dy((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DyShFncAtIP3[i]=ShapeFnc3dy((*IntPoints)[i][0],(*IntPoints)[i][1]);
    }
}

Vector<Double> & Triangle1::GetShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
    {
    case 1:
      return ShFncAtIP1;
    case 2:
      return ShFncAtIP2;
    case 3:
      return ShFncAtIP3;
    default:
    Error("Shape function does not exist with this number", __FILE__,__LINE__);
    }
}

}
