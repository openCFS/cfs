#include <iostream>
#include <fstream>

#include "dampint.hh"

namespace CoupledField
{

template <Integer Dim>
DampInt<Dim> :: DampInt(BaseElem * aptelem, const ShortInt ndofs) : BaseForm<Dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering DampInt::DampInt" << std::endl;
#endif
  DofsPerNode = ndofs;
}
  
template <Integer Dim>
DampInt<Dim> :: ~DampInt()
{
#ifdef TRACE
  (*trace) << "entering DampInt::~DampInt" << std::endl;
#endif

  ;
}

template <Integer Dim>
void DampInt<Dim> :: CalcElemMatrix(Point<Dim> * ptCoord, Matrix<Double> & Result)
{

  Integer l=ptelem->GetNumIntPoints(); // l - number of integration points
  Integer n=ptelem->GetNumNodes();   // n - number of nodes
  
  Jacobian<Dim> J;
  Integer i,ii,iii;
  Double length; 
  Result.Resize(n,n);

  Vector<Double> * Sf=new Vector<Double>[n];
  for (i=0; i < n; i++)
    Sf[i]=ptelem->GetShFncAtIP(i+1);

  Vector<Double> * intWeights=ptelem->GetIntWeights();
 
  for (i=0; i<l; i++)
    {
      // ptelem->CalcJacobian(J, i, ptCoord, FALSE);
      length=dist(ptCoord[0],ptCoord[1]);
      for (ii=0; ii < n; ii++)
	for (iii=0; iii<ii+1; iii++)
       {
        if  (intWeights)
	  Result[ii][iii]+=length/2.0*Sf[ii][i]*Sf[iii][i]*(*intWeights)[i];
         else
           Result[ii][iii]+=length/2.0*Sf[ii][i]*Sf[iii][i];
       }

    }

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] Sf;
}


}
