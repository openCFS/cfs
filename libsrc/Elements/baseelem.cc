#include <iostream>
#include <fstream>

//#include "general_head.hh" 
//#include "utils_head.hh"
#include "baseelem.hh"
 
namespace CoupledField
{

/*
template<> 
Jacobian<Point2D>::Jacobian()
{
  J.Resize(2,2);
  Jinv.Resize(2,2);
}

template<> 
Jacobian<Point3D>::Jacobian()
{
  J.Resize(3,3);
  Jinv.Resize(3,3);
}
*/

template <class Dim>
Jacobian<Dim>::Jacobian()
{
  Dim a;
 if (a.is2D()) { J.Resize(2,2);
                   Jinv.Resize(2,2); 
                 }
 else { J.Resize(3,3);
        Jinv.Resize(3,3);
 }
}

template<class Dim>
void Jacobian<Dim>::GetJinvX(Vector<Double> & JinvX)
{
  Dim a;
  if (a.is2D()) { JinvX.Resize(2);

  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];}
  else
    { 
  JinvX.Resize(3);

  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];
  JinvX[2]=Jinv[2][0];
    } 
}

/*
template<>
void Jacobian<Point2D>::GetJinvX(Vector<Double> & JinvX)
{
  JinvX.Resize(2);

  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];
}

template<>
void Jacobian<Point3D>::GetJinvX(Vector<Double> & JinvX)
{
  JinvX.Resize(3);

  JinvX[0]=Jinv[0][0];
  JinvX[1]=Jinv[1][0];
  JinvX[2]=Jinv[2][0];
}
*/

template<class Dim>
void Jacobian<Dim>::GetJinvY(Vector<Double> & JinvY)
{
  Dim a;
  if (a.is2D()) { JinvY.Resize(2); 

  JinvY[0]=Jinv[0][1];
  JinvY[1]=Jinv[1][1]; }
  else {  JinvY.Resize(3);

  JinvY[0]=Jinv[0][1];
  JinvY[1]=Jinv[1][1];
  JinvY[2]=Jinv[2][1];}
}
/*
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
*/

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
