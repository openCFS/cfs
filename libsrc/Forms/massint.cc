#include <iostream.h>
#include <fstream.h>

#include <general_head.hh>
#include <utils_head.hh>
#include "baseelem.hh"
#include "baseform.hh"
#include "massint.hh"

namespace CoupledField
{

template <class Dim>
MassInt<Dim> :: MassInt(BaseElem * aptelem, const ShortInt ndofs) : BaseForm<Dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering MassInt::MassInt" << endl;
#endif
  DofsPerNode = ndofs;
}
  
template <class Dim>
MassInt<Dim> :: ~MassInt()
{
#ifdef TRACE
  (*trace) << "entering MassInt::~MassInt" << endl;
#endif

  ;
}

template <class Dim>
void MassInt<Dim> :: CalcElemMatrix(Dim * ptCoord, Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering MassInt::CalcElemMatrix" << endl;
#endif

  Integer l=ptelem->GetNumIntPoints(); // l - number of integration points
  Integer n=ptelem->GetNumNodes();   // n - number of nodes
 
  Jacobian<Dim> J;
  Integer i,ii,iii;
 
  Result.Resize(n,n);
 
  Vector<Double> * Sf=new Vector<Double>[n];
  for (i=0; i < n; i++)
    Sf[i]=ptelem->GetShFncAtIP(i+1);

  Vector<Double> * intWeights=ptelem->GetIntWeights();
 
  for (i=0; i<l; i++)
    {
      ptelem->CalcJacobian(J, i, ptCoord, FALSE);
 
      for (ii=0; ii < n; ii++)
	for (iii=0; iii<ii+1; iii++)
       {
        if  (intWeights)
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


}
