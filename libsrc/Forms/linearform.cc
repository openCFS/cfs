#include <iostream>
#include <fstream>
#include <math.h>

#include "linearform.hh"
#include <Elements/jacobian.hh>
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
LinearForm<dim>:: LinearForm(BaseFE * aptelem) : BaseForm<dim>(aptelem)
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

// Doing it for acoustic RHS (dipole)
//! only for 2D
// template<>
// void LinearForm<2> :: CalcElemVector4FlowSrcDip(Point<2> * ptCoord,const Vector<Integer> & connecth, Vector<Double> & Result, const std::vector<Double> gradN_x_P)
// {
// #ifdef TRACE
//   (*trace) << "entering LinearForm::CalcElemVector4FlowSrc" << std::endl;
// #endif

//   Integer l=ptelem->GetNumIntPoints();
//   Integer n=ptelem->GetNumNodes();
//   Vector<Double> * help=new Vector<Double>[n];

       
//   Double det;
//   Integer i,ii;
//   Vector<Double> JinvX, JinvY;
//   Jacobian<2> J;
//   Double density=1.0;
//   Result.Resize(n);
//   Vector<Double> * Sf=new Vector<Double> [n];


//   // Dipole Source Term: (w,deltp)onbnd

//   std::vector<Double> normal;
//   normal.resize(2);
//   Double multiplier, dNdxNormalP, length;

//  calcNormal2Line(normal,ptCoord[0],ptCoord[1]);

//       for (i=0; i < n; i++)
//         {
//           Sf[i]=ptelem->GetShFncAtIP(i+1);
//         }

//       for (i=0; i<l; i++)
//         {
//           length=dist(ptCoord[0],ptCoord[1]);

//           for (ii=0; ii < n; ii++)
//             {
//                 multiplier=ScalarMult(normal,gradN_x_P); // n * gradShFnc
//                 Result[ii]+= -length/2.0*multiplier*Sf[ii][i];
//             }

//         }

// } // end of method

// Doing it for acoustic RHS (dipole)
template<>
void LinearForm<3> :: CalcElemVector4FlowSrcDip(Point<3> * ptCoord,const Vector<Integer> & connecth, Vector<Double> & Result, const std::vector<Double> gradN_x_P)
{
#ifdef TRACE
  (*trace) << "entering LinearForm::CalcElemVector4FlowSrc" << std::endl;
#endif

  Error("This function is not implemented for 3d",__FILE__,__LINE__);

} // end of method

// Doing it for acoustic RHS (quadrupole)
// template<Integer dim>
// void LinearForm<dim> :: CalcElemVector4FlowSrcQuad(Point<dim> * ptCoord,const Vector<Integer> & connecth,const Matrix<Double> & FlowData, Vector<Double> & Result)
// {
// #ifdef TRACE
//   (*trace) << "entering LinearForm::CalcElemVector4FlowSrc" << std::endl;
// #endif

//   Integer l=ptelem->GetNumIntPoints();
//   Integer n=ptelem->GetNumNodes();
//   Vector<Double> * help=new Vector<Double>[n];
       
//   Double det;
//   Integer i,ii;
//   Vector<Double> JinvX, JinvY, JinvZ;
//   Jacobian<dim> J;
//   Double density=1.0;
//   Result.Resize(n);
//   Vector<Double> * Sf=new Vector<Double> [n];

//   Double dTxxdx1part,dTxxdx2part,dTxydx1part,dTxydx2part;
//   Double dTyxdy1part,dTyxdy2part,dTyydy1part,dTyydy2part;
//   for (i=0; i < n; i++)
//     Sf[i]=ptelem->GetShFncAtIP(i+1);

//   if (dim==2)
//     {
//   for (i=0; i<l; i++)
//     {
//           ptelem->CalcJacobian(J,i,ptCoord);
//           J.GetJinvX(JinvX);
//           J.GetJinvY(JinvY);
//           // Resetting derivatives of tensor
//           dTxxdx1part=0;
//           dTxxdx2part=0;
//           dTxydx1part=0;
//           dTxydx2part=0;
//           dTyxdy1part=0;
//           dTyxdy2part=0;
//           dTyydy1part=0;
//           dTyydy2part=0;
//           // This work only for quad1 elements since we have values
//           // for ux and uy only at the corners!
//           // Here we compute the derivatives of the Lighthill's tensor needed in the quadrupole term
//           // at the ith integration point
//           for (ii=0; ii<n; ii++)
//             {
//               ptelem->GetGradientShFnc(help[ii],ii+1,i);

//               dTxxdx1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]); 
 
//               dTxxdx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];
//               dTxydx1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
//               dTxydx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];

//               dTyxdy1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
//               dTyxdy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
//               dTyydy1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
// 	      dTyydy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
//             }


//  // Setting up the RHS (dw/di,dTij/di) part
// 	  for (ii=0; ii < n; ii++) Result[ii]+= density*((help[ii]*JinvX)*(2.*dTxxdx1part*dTxxdx2part+dTxydx1part*dTxydx2part)+(help[ii]*JinvY)*(dTyxdy1part*dTyxdy2part+2.*dTyydy1part*dTyydy2part))*J.detJ;
//     }
//     } // end of if (dim==2)

//   if (dim==3)
//     {

//       Double dTxzdx1part, dTxzdx2part, dTyzdy1part, dTyzdy2part, dTzxdz1part, dTzxdz2part, dTzydz1part, dTzydz2part, dTzzdz1part, dTzzdz2part;
//       for (i=0; i<l; i++)
// 	{
// 	  ptelem->CalcJacobian(J,i,ptCoord);  
// 	  J.GetJinvX(JinvX);
// 	  J.GetJinvY(JinvY);
// 	  J.GetJinvZ(JinvZ);
// 	  // Resetting derivatives of tensor
// 	  dTxxdx1part=0;
// 	  dTxxdx2part=0;
// 	  dTxydx1part=0;
// 	  dTxydx2part=0;
// 	  dTxzdx1part=0;
// 	  dTxzdx2part=0;
	  
// 	  dTyxdy1part=0;
// 	  dTyxdy2part=0;
// 	  dTyydy1part=0;
// 	  dTyydy2part=0;
// 	  dTyzdy1part=0;
// 	  dTyzdy2part=0;

// 	  dTzxdz1part=0;
// 	  dTzxdz2part=0;
// 	  dTzydz1part=0;
// 	  dTzydz2part=0;
// 	  dTzzdz1part=0;
// 	  dTzzdz2part=0;

// 	  // This work only for trilinear hexahedral elements since we have values
// 	  // for ux uy and uz only at the corners!
// 	  // Here we compute the derivatives of the Lighthill's tensor needed in the quadrupole term
// 	  // at the ith integration point
// 	  for (ii=0; ii<n; ii++)
// 	    {
// 	      ptelem->GetGradientShFnc(help[ii],ii+1,i);

// 	      dTxxdx1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
// 	      dTxxdx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];
// 	      dTxydx1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
// 	      dTxydx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];
// 	      dTxzdx1part+=(FlowData[3][connecth[ii]-1]*Sf[ii][i]);
// 	      dTxzdx2part+=(help[ii]*JinvX)*FlowData[1][connecth[ii]-1];

// 	      dTyxdy1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
// 	      dTyxdy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
// 	      dTyydy1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
// 	      dTyydy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];
// 	      dTyzdy1part+=(FlowData[3][connecth[ii]-1]*Sf[ii][i]);
// 	      dTyzdy2part+=(help[ii]*JinvY)*FlowData[2][connecth[ii]-1];

// 	      dTzxdz1part+=(FlowData[1][connecth[ii]-1]*Sf[ii][i]);
// 	      dTzxdz2part+=(help[ii]*JinvZ)*FlowData[3][connecth[ii]-1];
// 	      dTzydz1part+=(FlowData[2][connecth[ii]-1]*Sf[ii][i]);
// 	      dTzydz2part+=(help[ii]*JinvZ)*FlowData[3][connecth[ii]-1];
// 	      dTzzdz1part+=(FlowData[3][connecth[ii]-1]*Sf[ii][i]);
// 	      dTzzdz2part+=(help[ii]*JinvZ)*FlowData[3][connecth[ii]-1];
// 	    }

// 	  // Setting up the RHS (dw/di,dTij/di) part
// 	  for (ii=0; ii < n; ii++) Result[ii]+= density*((help[ii]*JinvX)*(2.*dTxxdx1part*dTxxdx2part+dTxydx1part*dTxydx2part+dTxzdx1part*dTxzdx2part)+(help[ii]*JinvY)*(dTyxdy1part*dTyxdy2part+2.*dTyydy1part*dTyydy2part+dTyzdy1part*dTyzdy2part)+(help[ii]*JinvZ)*(dTzxdz1part*dTzxdz2part+dTzydz1part*dTzydz2part+2.*dTzzdz1part*dTzzdz2part))*J.detJ;
// 	}   
//     } // end of if (dim==3)

// } // end of method



} // end of namespace
