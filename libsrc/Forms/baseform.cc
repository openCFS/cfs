#include <iostream>
#include <fstream>

#include "baseform.hh"

namespace CoupledField
{

template <Integer dim>
BaseForm<dim> :: BaseForm(BaseElem * aptelem) 
{
#ifdef TRACE
  (*trace) << "entering BaseForm::BaseForm" << std::endl;
#endif

  ptelem=aptelem;

}
 
template <Integer dim> 
BaseForm<dim> :: ~BaseForm()
{
#ifdef TRACE
  (*trace) << "entering BaseForm::~BaseForm" << std::endl;
#endif
;
}

template <Integer dim>
void BaseForm<dim> :: CalcShapeFncDxShapeFncDx(Point<dim> * ptCoord, 
                                      Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering BaseForm::CalcShapeFncDyShapeFncDy" << std::endl;
#endif

  Integer l=ptelem->GetNumIntPoints(); // l - number of integration points
  Integer n=ptelem->GetNumNodes(); // n - number of nodes in element

  Jacobian<dim>  J;
  Vector<Double> JinvX;
 
  Vector<Double> * help=new Vector<Double>[n];
  Integer i,ii,iii;
 
  Result.Resize(n,n);
 

  Vector<Double> * intWeights=ptelem->GetIntWeights();  

  for (i=0; i<l; i++)
    {

      ptelem->CalcJacobian(J,i,ptCoord);

      J.GetJinvX(JinvX);

      for (ii=0; ii<n; ii++)
	{ 
          ptelem->GetGradientShFnc(help[ii],ii+1,i);
        }
 

      for (ii=0; ii < n; ii++)
	for (iii=0; iii<ii+1; iii++)
	  {
	    if (intWeights)  
	      Result[ii][iii]+=(help[ii]*JinvX)*(help[iii]*JinvX)*J.detJ*(*intWeights)[i];
	    else
	      Result[ii][iii]+=(help[ii]*JinvX)*(help[iii]*JinvX)*J.detJ;
	  }

    }

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] help;
}

template <Integer dim>
void BaseForm<dim> :: CalcShapeFncDyShapeFncDy(Point<dim> * ptCoord,
                           Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering BaseForm::CalcShapeFncDyShapeFncDy" << std::endl;
#endif
  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();

  Jacobian<dim>  J;
  Vector<Double> JinvY;
 
  Vector<Double> * help=new Vector<Double>[n];
  Integer i,ii,iii;

  Result.Resize(n,n);

  Vector<Double> * intWeights=ptelem->GetIntWeights(); 

  for (i=0; i<l; i++)
    {
  
      ptelem->CalcJacobian(J,i,ptCoord);  

      J.GetJinvY(JinvY);

      for (ii=0; ii<n; ii++)
        {
          ptelem->GetGradientShFnc(help[ii],ii+1,i);
        }

      for (ii=0; ii < n; ii++)
	for (iii=0; iii<ii+1; iii++)
	  {
	    if (intWeights) 
	      Result[ii][iii]+=(help[ii]*JinvY)*(help[iii]*JinvY)*J.detJ*(*intWeights)[i];
	    else
	      Result[ii][iii]+=(help[ii]*JinvY)*(help[iii]*JinvY)*J.detJ;
	  }

    }
 
  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] help;
}

template <Integer dim>
void BaseForm<dim> :: CalcShapeFncShapeFnc(Point<dim> * ptCoord,
                                 Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering BaseForm::CalcShapeFncShapeFnc" << std::endl;
#endif

  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();

  Jacobian<dim> J; 
  Integer i,ii,iii;
 
  Result.Resize(n,n);
 
  Vector<Double> * Sf=new Vector<Double> [n];
  for (i=0; i < n; i++)
    Sf[i]=ptelem->GetShFncAtIP(i+1);

  Vector<Double> * intWeights=ptelem->GetIntWeights();  
 
  for (i=0; i<l; i++)
    {
      ptelem->CalcJacobian(J, i, ptCoord, FALSE);
 
      for (ii=0; ii < n; ii++)
	for (iii=0; iii<ii+1; iii++)
	  { 
	    if (intWeights)	  
	      Result[ii][iii]+=J.detJ*Sf[ii][i]*Sf[iii][i]*(*intWeights)[i];
	    else
	      Result[ii][iii]+=J.detJ*Sf[ii][i]*Sf[iii][i];
	  }

    }

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] Sf;
}


template <Integer dim>
Double BaseForm<dim> :: FuncAtIP(const ShortInt iIP, RHS f)
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::FuncAtIP" << std::endl;
#endif

  Double result;

  switch(dim)
    { 
    case 2:
      Double x,y;

      x=ptelem->GetIntPointsX(iIP);
      y=ptelem->GetIntPointsY(iIP);
      result=f(x,y,0);

      break;

    case 3:
      Error("3D has not implemented yet",__FILE__,__LINE__);
    }

  return result;
}

template <Integer dim>
void BaseForm<dim>::CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> & StiffMat) 
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::CalcElemMatrix" << std::endl;
#endif
  Error(" Function CalcElemMatrix of BilinieaForm is virtual. You can use it for derived classes LaplaceInt or MassInt.",__FILE__,__LINE__);
}

template <Integer dim>
void BaseForm<dim>::CalcElemVector(Point<dim> * ptCoord, Vector<Double> & StiffMat) 
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::CalcElemVector" << std::endl;
#endif
  Error(" Function CalcElemVector of BaseForm is virtual. You can use it for derived class LinearForm.",__FILE__,__LINE__);
}

template <Integer dim>
void BaseForm<dim>::CalcElemVector4InterpolatedFnc(Point<dim> * ptCoord, const Integer aComponent, Vector<Double> & aValueAtNodePoints, Vector<Double> & Result)
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::CalcElemVector" << std::endl;
#endif
  Error(" Function CalcElemVector4InterpolatedFnc of BaseForm is virtual. You can use it for derived class LinearForm.",__FILE__,__LINE__);
}

template <Integer dim>
void BaseForm<dim>::Print(std::ostream * out, const Matrix<Double> Result) const
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::Print" << std::endl;
#endif
 Error(" Function Print of BilinieaForm is virtual. You can use it for derived  classes LaplaceInt or MassInt.",__FILE__,__LINE__);
}



}
