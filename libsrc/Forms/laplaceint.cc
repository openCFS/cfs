#include <iostream.h>
#include <fstream.h>

#include <general_head.hh>
#include <utils_head.hh>
#include "baseelem.hh"
#include "baseform.hh"
#include "laplaceint.hh"

namespace CoupledField
{

template <class Dim>
LaplaceInt<Dim> :: LaplaceInt(BaseElem * aptelem, const ShortInt ndofs)
              : BaseForm<Dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::LaplaceInt" << endl;
#endif

  DofsPerNode = ndofs;
}
 
template <class Dim> 
LaplaceInt<Dim> :: ~LaplaceInt()
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::~LaplaceInt" << endl;
#endif
 

  ;
}

template <class Dim>
void LaplaceInt<Dim> :: CalcElemMatrix(Dim * ptCoord, Matrix<Double> & Result)   
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::CalcElemMatrix" << endl;
#endif
  
  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();
 
  Jacobian<Dim>  J;
  Vector<Double> JinvX, JinvY;
 
  Vector<Double> * help=new Vector<Double>[n];
  Integer i,ii,iii;
 
  Result.Resize(n,n);

  Vector<Double> * intWeights=ptelem->GetIntWeights();  

  for (i=0; i<l; i++)
    {
 
      ptelem->CalcJacobian(J,i,ptCoord);

      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);
   
      for (ii=0; ii<n; ii++)
        {
          ptelem->GetGradientShFnc(help[ii],ii+1,i);
        }
 
      for (ii=0; ii < n; ii++)
        for (iii=0; iii<ii+1; iii++)
        {
          if (intWeights) 
	  Result[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX)
            +(help[ii]*JinvY)*(help[iii]*JinvY))*J.detJ*(*intWeights)[i];
          else  
         Result[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX)
                 +(help[ii]*JinvY)*(help[iii]*JinvY))*J.detJ;
        }
    }

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];
 
  delete [] help;
}

template <class Dim>
void LaplaceInt<Dim>::Print(ostream * out, const Matrix<Double> Result) const
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::Print" << endl;
#endif
  (*out)<< "Stiffness matrix" << endl;
  (*out) << Result;

}
}
