#include <iostream>
#include <fstream>

#include "baseelem.hh"
 
namespace CoupledField
{

template <>
Jacobian<Point2D>::Jacobian()
{
  J.Resize(2,2);
  Jinv.Resize(2,2); 
}

template <>
Jacobian<Point3D>::Jacobian()
{
  J.Resize(3,3);
  Jinv.Resize(3,3);
}

template<>
void Jacobian<Point3D>::GetJinvX(Vector<Double> & JinvX)
{
  JinvX.Resize(3);

  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];
  JinvX[2]=Jinv[2][0];
}

template<>
void Jacobian<Point2D>::GetJinvX(Vector<Double> & JinvX)
{
  JinvX.Resize(2);
  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];
}

template<>
void Jacobian<Point2D>::GetJinvY(Vector<Double> & JinvY)
{
  JinvY.Resize(2); 

  JinvY[0]=Jinv[0][1];
  JinvY[1]=Jinv[1][1]; 
}

template<>
void Jacobian<Point3D>::GetJinvY(Vector<Double> & JinvY)
{
  JinvY.Resize(3);

  JinvY[0]=Jinv[0][1];
  JinvY[1]=Jinv[1][1];
  JinvY[2]=Jinv[2][1];
}
//----------------------------------------------------------------------------
 BaseElem :: BaseElem()
{
#ifdef TRACE
  (*trace) << "entering BaseElem::BaseElem" << std::endl;
#endif
  ;
}
 
BaseElem :: ~BaseElem()
{
#ifdef TRACE
  (*trace) << "entering BaseElem::~BaseElem" << std::endl;
#endif
 
  ;
}

}
