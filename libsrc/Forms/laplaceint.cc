#include <iostream>
#include <fstream>

#include "laplaceint.hh"

namespace CoupledField
{

template <Integer dim>
LaplaceInt<dim> :: LaplaceInt(BaseFE * aptelem, const ShortInt ndofs)
              : BaseForm<dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::LaplaceInt" << std::endl;
#endif

  DofsPerNode = ndofs;
}
 
template <Integer dim> 
LaplaceInt<dim> :: ~LaplaceInt()
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::~LaplaceInt" << std::endl;
#endif
 
  ;
}

template <Integer dim>
void LaplaceInt<dim> :: CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> & Result)  
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::CalcElemMatrix" << std::endl;
#endif
  
  Integer l=ptelem->GetNumIntPoints();

  Integer n=ptelem->GetNumNodes();

  Jacobian<dim>  J;
  Vector<Double> JinvX, JinvY, JinvZ;
 
  Vector<Double> * help=new Vector<Double>[n];
  Integer i,ii,iii;
 
  Result.Resize(n,n);

  Vector<Double> * intWeights=ptelem->GetIntWeights();  

  for (i=0; i<l; i++)
    {

      ptelem->CalcJacobian(J,i,ptCoord);

      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);
      if (dim==3) J.GetJinvZ(JinvZ);
   
      for (ii=0; ii<n; ii++)
        {
          ptelem->GetGradientShFnc(help[ii],ii+1,i);
        }
 
      if (dim==2) {
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

      if (dim==3) {      
      for (ii=0; ii < n; ii++)
        for (iii=0; iii<ii+1; iii++)
        {
          if (intWeights)
          Result[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX)
            +(help[ii]*JinvY)*(help[iii]*JinvY)+(help[ii]*JinvZ)*(help[iii]*JinvZ))*J.detJ*(*intWeights)[i];
          else
         Result[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX)
                 +(help[ii]*JinvY)*(help[iii]*JinvY)+(help[ii]*JinvZ)*(help[iii]*JinvZ))*J.detJ;
        }
      }

    }

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];
 
  delete [] help;

#ifdef TRACE
  (*trace) << "leaving LaplaceInt::CalcElemMatrix" << std::endl;
#endif
}

// template <>
// void LaplaceInt<3> :: CalcElemMatrix(Point<3> * ptCoord, Matrix<Double> & Result)
// {
// #ifdef TRACE
//   (*trace) << "entering LaplaceInt::CalcElemMatrix" << std::endl;
// #endif

//   Integer l=ptelem->GetNumIntPoints();

//   Integer n=ptelem->GetNumNodes();

//   Jacobian<3>  J;
//   Vector<Double> JinvX, JinvY, JinvZ;

//   Vector<Double> * help=new Vector<Double>[n];
//   Integer i,ii,iii;

//   Result.Resize(n,n);

//   Vector<Double> * intWeights=ptelem->GetIntWeights();

//   for (i=0; i<l; i++)
//     {
//       ptelem->CalcJacobian(J,i,ptCoord);

//       J.GetJinvX(JinvX);
//       J.GetJinvY(JinvY);
//       J.GetJinvZ(JinvZ);

//       for (ii=0; ii<n; ii++)
//         {
//           ptelem->GetGradientShFnc(help[ii],ii+1,i);
//         }

//       for (ii=0; ii < n; ii++)
//         for (iii=0; iii<ii+1; iii++)
//         {
//           if (intWeights)
//           Result[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX)
//             +(help[ii]*JinvY)*(help[iii]*JinvY)+(help[ii]*JinvZ)*(help[iii]*JinvZ))*J.detJ*(*intWeights)[i];
//           else
//          Result[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX)
//                  +(help[ii]*JinvY)*(help[iii]*JinvY)+(help[ii]*JinvZ)*(help[iii]*JinvZ))*J.detJ;
//         }
//     }

//   for (ii=0; ii<n; ii++)
//     for (iii=0; iii<ii; iii++)
//       Result[iii][ii]=Result[ii][iii];

//   delete [] help;
// }

/*
template <>
void LaplaceInt<Point<3>> :: CalcElemMatrix(Point3D * ptCoord, Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::CalcElemMatrix" << endl;
#endif

  Integer l=ptelem->GetNumIntPoints();
  Integer n=ptelem->GetNumNodes();

  Double * help1=new Double[n];
  Double * help2=new Double[n];
  Double * help3=new Double[n];

  Double det;
  Integer i,ii,iii;

  Result.Resize(n,n);

  Vector<Double> ** Vx=new Vector<Double>*[n];
  Vector<Double> ** Vy=new Vector<Double>*[n];
  Vector<Double> ** Vz=new Vector<Double>*[n];

  for (i=0; i < n; i++)
  Vx[i]=ptelem->GetDxShFncAtIP(i+1);

  for (i=0; i < n; i++)
  Vy[i]=ptelem->GetDyShFncAtIP(i+1);

  for (i=0; i<n; i++)
  Vz[i]=ptelem->GetDzShFncAtIP(i+1);

  for (i=0; i<l; i++)
{

   for (ii=0; ii<n; ii++)
{
help1[ii] =
(*Vx[ii])[i]*(TransYDnu(i,ptCoord)*TransZDgam(i,ptCoord)-TransYDgam(i,ptCoord)*TransZDnu(i,ptCoord))+
(*Vy[ii])[i]*(-TransYDxi(i,ptCoord)*TransZDgam(i,ptCoord)+TransYDgam(i,ptCoord)*TransZDxi(i,ptCoord))+
(*Vz[ii])[i]*(TransYDxi(i,ptCoord)*TransZDnu(i,ptCoord)-TransYDnu(i,ptCoord)*TransZDxi(i,ptCoord));

help2[ii] =
(*Vx[ii])[i]*(TransXDgam(i,ptCoord)*TransZDnu(i,ptCoord)-TransXDnu(i,ptCoord)*TransZDgam(i,ptCoord))+
(*Vy[ii])[i]*(TransXDxi(i,ptCoord)*TransZDgam(i,ptCoord)-TransXDgam(i,ptCoord)*TransZDxi(i,ptCoord))+
(*Vz[ii])[i]*(-TransXDxi(i,ptCoord)*TransZDnu(i,ptCoord)+TransXDnu(i,ptCoord)*TransZDxi(i,ptCoord));

help3[ii] =
(*Vx[ii])[i]*(TransXDnu(i,ptCoord)*TransYDgam(i,ptCoord)-TransXDgam(i,ptCoord)*TransYDnu(i,ptCoord))+
(*Vy[ii])[i]*(-TransXDxi(i,ptCoord)*TransYDgam(i,ptCoord)+TransXDgam(i,ptCoord)*TransYDxi(i,ptCoord))+
(*Vz[ii])[i]*(TransXDxi(i,ptCoord)*TransYDnu(i,ptCoord)-TransXDnu(i,ptCoord)*TransYDxi(i,ptCoord));
}

   det=CalcDet(i,ptCoord);
   for (ii=0; ii < n; ii++)
     for (iii=0; iii<ii+1; iii++)
        Result[ii][iii]+=(1.0/det)*(help1[ii]*help1[iii]+help2[ii]*help2[iii]+help3[ii]*help3[iii]);
   for (ii=0; ii<n; ii++)
     for (iii=0; iii<ii; iii++)
        Result[iii][ii]=Result[ii][iii];
 }
}
*/

template <Integer dim>
void LaplaceInt<dim>::Print(std::ostream * out, const Matrix<Double> Result) const
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::Print" << std::endl;
#endif
  (*out)<< "Stiffness matrix" << std::endl;
  (*out) << Result;

}
}
