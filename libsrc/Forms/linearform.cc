#include <iostream>
#include <fstream>
#include <math.h>

#include "linearform.hh"

namespace CoupledField
{

// Double RhsFunc(Double x, Double y, Double t)
// {
//   Double alpha=0.8;
//   Double xi0=20;
//   Double xi=(x+y)/sqrt(2);
  
//   Double part1=-2*(y-y*y+x-x*x)/tan(alpha*(xi-xi0));
//   Double part2=-2*alpha*(1/sqrt(2))*((1-2*x)*(y-y*y)+(1-2*y)*(x-x*x))/sqr(sin(alpha*(xi-xi0)));
//   Double part3=2*sqr(alpha)*(x-sqr(x))*(y-sqr(y))*cos(alpha*(xi-xi0))/(sqr(sin(alpha*(xi-xi0)))*sin(alpha*(xi-xi0)));

//   Double res=part1+part2+part3;
//   return res;
//   // return 0;
// } 

 Double RhsFunc(Double x, Double y, Double t)
   {
     // return -2*(x+y-x*x-y*y);
     return 0;
 //     return x;
  }

template <Integer dim>
LinearForm<dim>:: LinearForm(BaseElem * aptelem) : BaseForm<dim>(aptelem)
{
#ifdef TRACE
  (*trace) << "entering LinearForm::LinearForm" << std::endl;
#endif
  ;
}
  
template <Integer dim>
LinearForm<dim> :: ~LinearForm()
{
#ifdef TRACE
  (*trace) << "entering LinearForm::~LinearForm" << std::endl;
#endif

  ;
}

template<Integer dim>
void LinearForm<dim> :: CalcElemVector(Point<dim> * ptCoord, Vector<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering LinearForm::CalcElemVector" << std::endl;
#endif

  // get number of integr-points
  Integer noIntPnt=ptelem->GetNumIntPoints();
  // get number of nodes in standard element
  Integer noNd=ptelem->GetNumNodes();

  Jacobian<dim> J;

  Result.Resize(noNd);

  // declaration of shape functions
  Vector<Double> * Sf=new Vector<Double> [noIntPnt];
 
  // get values of shape functions at integr-points
  Integer ish;
  for (ish=0; ish < noNd; ish++)
    Sf[ish]=ptelem->GetShFncAtIP(ish+1);

  // loop over integr-points
  Vector<Double> * intWeights; 
  intWeights=ptelem->GetIntWeights();
  Integer iIP;
  for (iIP=0; iIP<noIntPnt; iIP++)
    {
      // calculation of integral at integr-point
      ptelem->CalcJacobian(J,iIP,ptCoord, FALSE);  

      // get value of function RHSFunc at integr-point. this method is in baseform.hh
   Point<2> Point;
   Point[0] = 0;
   Point[1] = 0;
   for (ish=0; ish<noNd ; ish++) {
     Point[0] += Sf[ish][iIP]*ptCoord[ish][0];
     Point[1] += Sf[ish][iIP]*ptCoord[ish][1];
       }

   Double help=RhsFunc(Point[0],Point[1],0); // 0 - time
 
      // loop over shape fncs
   for (ish=0; ish < noNd; ish++) {

     if (intWeights) 
       Result[ish]+=J.detJ*help*Sf[ish][iIP]*(*intWeights)[iIP];
     else 
       Result[ish]+=J.detJ*help*Sf[ish][iIP];    
    }

    }
}

template<Integer dim>
void LinearForm<dim> :: CalcElemVector4InterpolatedFnc(Point<dim> * ptCoord, const Integer aComponent, Vector<Double> & aValFncAtNode, Vector<Double> & Result)
{
#ifdef TRACE
  (*trace) << "entering LinearForm::CalcElemVector4InterpolatedFnc" << std::endl;
#endif
  Integer noIntgrPnts=ptelem->GetNumIntPoints();
  Integer noNds=ptelem->GetNumNodes();

  Double det,fnc=0;
  Integer i,ii;

  Jacobian<dim> J;

  Result.Resize(noNds);

  Vector<Double> * Sf=new Vector<Double> [noNds];
 
  // get values of shape fncs at integr. points
  for (i=0; i < noNds; i++)
    Sf[i]=ptelem->GetShFncAtIP(i+1);
  // get integration weights
  Vector<Double> * intWeights; 
  intWeights=ptelem->GetIntWeights();

  Vector<Double> valFncRHS;
  valFncRHS.Resize(noIntgrPnts); 

  Vector<Double> vec,grad,glgrad;
  Integer ip,isf,idim;
  Vector<Double> JinvX,JinvY,JinvZ;

  for (ip=0; ip<noIntgrPnts; ip++) 
    { // loop over integration points
    
      ptelem->CalcJacobian(J,ip,ptCoord, TRUE);
  
      J.GetJinvX(JinvX);
      J.GetJinvY(JinvY);
      if (dim==3)  J.GetJinvZ(JinvZ);

      vec.Resize(dim);
      glgrad.Resize(dim);

      for (isf=0; isf<noNds; isf++) 
	{ // over shape fnc
	  ptelem->GetGradientShFnc(grad,isf+1,ip); // ii-fnc, i-ip

	  glgrad[0]=JinvX*grad;
	  glgrad[1]=JinvY*grad;
	  if (dim==3)  glgrad[2]=JinvZ*grad;

	  for (idim=0; idim<grad.size();idim++) 
	    {// loop over dimension
	      vec[idim]+=aValFncAtNode[isf]*glgrad[idim];
	    }
	}
      // valFncInRHS[ip]=vec.normL2(); // abs val
      if (aComponent==0)
	valFncRHS[ip]=vec[0]; // Dx(ShFnc)
      if (aComponent==1)
        valFncRHS[ip]=vec[1]; // Dy(ShFnc)
      if (aComponent==2)
	   valFncRHS[ip]=vec[2]; //Dz(ShFnc)

      for (ii=0; ii < noNds; ii++) {
	// loop over shape fnc
	if (intWeights) {
	 Result[ii]+=J.detJ*Sf[ii][ip]*valFncRHS[ip]*(*intWeights)[ip];
	}
	else {
	  Result[ii]+=J.detJ*Sf[ii][ip]*valFncRHS[ip];  
	}
      }
    }
}

// template<>
// void LinearForm<3> :: CalcElemVector4InterpolatedFnc(Point<3> * ptCoord, const Integer aComponent, Vector<Double> & aValFncAtNode, Vector<Double> & Result)
// {
// #ifdef TRACE
//   (*trace) << "entering LinearForm::CalcElemVector4InterpolatedFnc" << std::endl;
// #endif
//   Integer noIntgrPnts=ptelem->GetNumIntPoints();
//   Integer noNds=ptelem->GetNumNodes();

//   Double det,fnc=0;
//   Integer i,ii;

//   Jacobian<3> J;

//   Result.Resize(noNds);

//   Vector<Double> * Sf=new Vector<Double> [noNds];
 
//   // get values of shape fncs at integr. points
//   for (i=0; i < noNds; i++)
//     Sf[i]=ptelem->GetShFncAtIP(i+1);
//   // get integration weights
//   Vector<Double> * intWeights; 
//   intWeights=ptelem->GetIntWeights();

//   Vector<Double> valFncRHS;
//   valFncRHS.Resize(noIntgrPnts); 

//   Vector<Double> vec,grad,glgrad;
//   Integer ip,isf,idim;
//   Vector<Double> JinvX,JinvY,JinvZ;

//   for (ip=0; ip<noIntgrPnts; ip++) 
//     { // loop over integration points
    
//       ptelem->CalcJacobian(J,ip,ptCoord, TRUE);
  
//       J.GetJinvX(JinvX);
//       J.GetJinvY(JinvY);
//       J.GetJinvZ(JinvZ);

//       vec.Resize(3);
//       glgrad.Resize(3);

//       for (isf=0; isf<noNds; isf++) 
// 	{ // over shape fnc
// 	  ptelem->GetGradientShFnc(grad,isf+1,ip); // ii-fnc, i-ip

// 	  glgrad[0]=JinvX*grad;
// 	  glgrad[1]=JinvY*grad;
// 	  glgrad[2]=JinvZ*grad;

// 	  for (idim=0; idim<grad.size();idim++) 
// 	    {// loop over dimension
// 	      vec[idim]+=aValFncAtNode[isf]*glgrad[idim];
// 	    }
// 	}
//       // valFncInRHS[ip]=vec.normL2(); // abs val
//       if (aComponent==0)
// 	valFncRHS[ip]=vec[0]; // Dx(ShFnc)
//       if (aComponent==1)
//         valFncRHS[ip]=vec[1]; // Dy(ShFnc)
//       if (aComponent==2)
//         valFncRHS[ip]=vec[2]; // Dz(ShFnc)

//       for (ii=0; ii < noNds; ii++) {
// 	// loop over shape fnc
// 	if (intWeights) {
// 	 Result[ii]+=J.detJ*Sf[ii][ip]*valFncRHS[ip]*(*intWeights)[ip];
// 	}
// 	else {
// 	  Result[ii]+=J.detJ*Sf[ii][ip]*valFncRHS[ip];  
// 	}
//       }
//     }
// }

} // end of namespace
