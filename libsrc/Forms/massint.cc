#include <iostream>
#include <fstream>

#include "massint.hh"
#include <Elements/jacobian.hh>
#include <Elements/basefe.hh>

namespace CoupledField
{

template <Integer dim>
MassInt<dim> :: MassInt(BaseFE * aptelem, const ShortInt ndofs) : BaseForm<dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering MassInt::MassInt" << std::endl;
#endif
  DofsPerNode = ndofs;
}
  
template <Integer dim>
MassInt<dim> :: ~MassInt()
{
#ifdef TRACE
  (*trace) << "entering MassInt::~MassInt" << std::endl;
#endif

  ;
}

template <Integer dim>
void MassInt<dim> :: CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> & Result)
{

  Integer l=ptelem->GetNumIntPoints(); // l - number of integration points
  Integer n=ptelem->GetNumNodes();   // n - number of nodes
 
  Jacobian<dim> Jii;
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
