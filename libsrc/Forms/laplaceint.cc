#include <iostream>
#include <fstream>

#include "laplaceint.hh"

namespace CoupledField
{

template <class Dim>
LaplaceInt<Dim> :: LaplaceInt(BaseElem * aptelem, const ShortInt ndofs)
              : BaseForm<Dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::LaplaceInt" << std::endl;
#endif

  DofsPerNode = ndofs;
}
 
template <class Dim> 
LaplaceInt<Dim> :: ~LaplaceInt()
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::~LaplaceInt" << std::endl;
#endif
 
  ;
}

template <>
void LaplaceInt<Point2D> :: CalcElemMatrix(Point2D * ptCoord, Matrix<Double> & Result)  
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::CalcElemMatrix" << std::endl;
#endif
  
  Integer l=ptelem->GetNumIntPoints();

  Integer n=ptelem->GetNumNodes();

  Jacobian<Point2D>  J;
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

template <>
void LaplaceInt<Point3D> :: CalcElemMatrix(Point3D * ptCoord, Matrix<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::CalcElemMatrix" << std::endl;
#endif

  Integer l=ptelem->GetNumIntPoints();

  Integer n=ptelem->GetNumNodes();

  Jacobian<Point3D>  J;
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
      J.GetJinvZ(JinvZ);

      std::cout << JinvX << std::endl;
      std::cout << JinvY << std::endl;
      std::cout << JinvZ << std::endl;

      for (ii=0; ii<n; ii++)
        {
          ptelem->GetGradientShFnc(help[ii],ii+1,i);
        }

       Matrix<Double> Test;
       Test.Resize(n,n);

       for (ii=0; ii < n; ii++)
       for (iii=0; iii<ii+1; iii++)
       {
          if (intWeights)
          Test[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX))*J.detJ*(*intWeights)[i];
          else
          Test[ii][iii]+=((help[ii]*JinvX)*(help[iii]*JinvX))*J.detJ;
       }

       std::cout << " Test " << Test << std::endl;

       for (ii=0; ii < n; ii++)
       for (iii=0; iii<ii+1; iii++)
       {
          if (intWeights)
          Test[ii][iii]+=((help[ii]*JinvY)*(help[iii]*JinvY))*J.detJ*(*intWeights)[i];
          else
          Test[ii][iii]+=((help[ii]*JinvY)*(help[iii]*JinvY))*J.detJ;
       }

       std::cout << " Test " << Test << std::endl;

       for (ii=0; ii < n; ii++)
       for (iii=0; iii<ii+1; iii++)
       {
          if (intWeights)
          Test[ii][iii]+=((help[ii]*JinvZ)*(help[iii]*JinvZ))*J.detJ*(*intWeights)[i];
          else
          Test[ii][iii]+=((help[ii]*JinvZ)*(help[iii]*JinvZ))*J.detJ;
       }

       std::cout << " Test " << Test << std::endl;

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

  for (ii=0; ii<n; ii++)
    for (iii=0; iii<ii; iii++)
      Result[iii][ii]=Result[ii][iii];

  delete [] help;
}

/*
template <>
void LaplaceInt<Point3D> :: CalcElemMatrix(Point3D * ptCoord, Matrix<Double> & Result)
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

template <class Dim>
void LaplaceInt<Dim>::Print(std::ostream * out, const Matrix<Double> Result) const
{
#ifdef TRACE
  (*trace) << "entering LaplaceInt::Print" << std::endl;
#endif
  (*out)<< "Stiffness matrix" << std::endl;
  (*out) << Result;

}
}
