#include <iostream.h>
#include <fstream.h>

#include <general_head.hh>
#include <utils_head.hh>
#include <elements_head.hh>
#include "baseform.hh"

namespace CoupledField
{

template <class Dim>
BaseForm<Dim> :: BaseForm(BaseElem * aptelem) 
{
#ifdef TRACE
  (*trace) << "entering BaseForm::BaseForm" << endl;
#endif

  ptelem=aptelem;

}
 
template <class Dim> 
BaseForm<Dim> :: ~BaseForm()
{
#ifdef TRACE
  (*trace) << "entering BaseForm::~BaseForm" << endl;
#endif
;
}

template <class Dim>
void BaseForm<Dim> :: CalcShapeFncDxShapeFncDx(Dim * ptCoord, 
                                      Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering BaseForm::CalcShapeFncDyShapeFncDy" << endl;
#endif

  Integer l=ptelem->GetNumIntPoints(); // l - number of integration points
  Integer n=ptelem->GetNumNodes(); // n - number of nodes in element

  Jacobian<Dim>  J;
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
          Result[ii][iii]+=(help[ii]*JinvX)*(help[iii]*JinvX)*J.detJ;

      if (intWeights)
	{
	  Result[ii][iii]*=(*intWeights)[i];
	}

    }

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] help;
}

template <class Dim>
void BaseForm<Dim> :: CalcShapeFncDyShapeFncDy(Dim * ptCoord,
                           Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering BaseForm::CalcShapeFncDyShapeFncDy" << endl;
#endif
  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();

  Jacobian<Dim>  J;
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
	  Result[ii][iii]+=(help[ii]*JinvY)*(help[iii]*JinvY)*J.detJ;

      if (intWeights)
	{
	  Result[ii][iii]*=(*intWeights)[i];
	}

    }
 
  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] help;
}

template <class Dim>
void BaseForm<Dim> :: CalcShapeFncShapeFnc(Dim * ptCoord,
                                 Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering BaseForm::CalcShapeFncShapeFnc" << endl;
#endif

  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();

  Jacobian<Dim> J; 
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
	  Result[ii][iii]+=J.detJ*Sf[ii][i]*Sf[iii][i];

      if (intWeights)
	{
	  Result[ii][iii]*=(*intWeights)[i];
	}

    }

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] Sf;
}


template <class Dim>
Double BaseForm<Dim> :: FuncAtIP(const ShortInt iIP, RHS f)
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::FuncAtIP" << endl;
#endif

  Double result;

  Integer dim=ptelem->GetDim();

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

template <class Dim>
void BaseForm<Dim>::CalcElemMatrix(Dim * ptCoord, Matrix<Double> & StiffMat) 
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::CalcElemMatrix" << endl;
#endif
  Error(" Function CalcElemMatrix of BilinieaForm is virtual. You can use it for derived classes LaplaceInt or MassInt.",__FILE__,__LINE__);
}

template <class Dim>
void BaseForm<Dim>::Print(ostream * out, const Matrix<Double> Result) const
{
#ifdef TRACE
  (*trace) <<  "entering BaseForm::Print" << endl;
#endif
 Error(" Function Print of BilinieaForm is virtual. You can use it for derived  classes LaplaceInt or MassInt.",__FILE__,__LINE__);
}
}
