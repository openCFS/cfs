#include <iostream>
#include <fstream>

#include "linearform.hh"

namespace CoupledField
{

Double RhsFunc(Double x, Double y, Double t){ return 0;} 

template <class Dim>
LinearForm<Dim>:: LinearForm(BaseElem * aptelem) : BaseForm<Dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering LinearForm::LinearForm" << std::endl;
#endif
  ;
}
  
template <class Dim>
LinearForm<Dim> :: ~LinearForm()
{
#ifdef TRACE
  (*trace) << "entering LinearForm::~LinearForm" << std::endl;
#endif

  ;
}

template<class Dim>
void LinearForm<Dim> :: CalcElemMatrix(Dim * ptCoord, Vector<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering LinearForm::~LinearForm" << std::endl;
#endif
  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();

  Double det,help;
  Integer i,ii;

  Jacobian<Dim> J;

  Result.Resize(n);

  Vector<Double> * Sf=new Vector<Double> [n];
 
  for (i=0; i < n; i++)
    Sf[i]=ptelem->GetShFncAtIP(i+1);

  for (i=0; i<l; i++)
    {
      ptelem->CalcJacobian(J,i,ptCoord, FALSE);  
 
      help=FuncAtIP(i, RhsFunc );

      for (ii=0; ii < n; ii++) Result[ii]+=J.detJ*help*Sf[ii][i];
    }
}

//template void BilinearForm<Point2D>:: FuncAtIP<RHSfunction>(const ShortInt, RHSfunction();
//template void BilinearForm<Point3D>:: FuncAtIP<RHSfunction>(const ShortInt, RHSfunction());

} // end of namespace
