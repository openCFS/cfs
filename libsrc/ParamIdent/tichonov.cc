#include <iostream>
#include <fstream>
#include <string>
//#include "staticdriver.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "CoupledPDE/basecoupledpde.hh"
#include "General/environment.hh"
#include "PDE/basePDE.hh" 

#include "piezoParamIdent.hh"
#include "Forms/baseForm.hh"
#include "Utils/vector.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/MaterialData.hh"
#include "PDE/timestepping.hh"
#include "Utils/baseelemstoresol.hh"
#include "Driver/singleDriver.hh"
#include "PDE/nodeEQN.hh"
#include <Domain/elem.hh>
#include "Forms/forms_header.hh"
#include "Utils/mathfunctions.hh"


#ifdef __sgi
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#define POW pow
#else
#include <cstdarg>
#include <cstdio>
#include <cmath>
#define POW std::pow
#endif

#include <sstream>
#include <iomanip>



#include "Utils/tools.hh"
#include <PDE/pdes_header.hh>

#include <Driver/piezoParamIdent.hh>



//#include "/../OLAS/algsys/basesystem.hh"
//#include "DataInOut/piezoParameterData.hh"



namespace CoupledField
{

  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxx tichonov xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  void piezoParamIdent::tichonov(){   // here we scale in advance, just once before NewtonCG and afterwards
    ENTER_FCN("piezoParamIdent::tichonov");
    std::cout<<"\n Entering piezoParamIdent::tichonov ... "<< std::endl;

    Integer nrNewtonIterations=0;
    Integer backtrackIterator=0;
    Integer i;
    Double misfit, misfit_new, lin_misfit;

    MaterialData * ptMaterial=pdes_[0]->getPDEMaterialData();   // Pointer to MaterialData

    ptBCs = pdes_[0]->getPDE_BCs();
    
    Vector<Complex> y_hat_F_hat (nrMeasuredData);
    Vector<Complex> JacFS(nrMeasuredData);
    Vector<Complex> linres(nrMeasuredData);
    Vector<Complex> basC;
    bas.Resize(actNrParameter);
    basC.Resize(actNrParameter);
    Vector<Double> stepR(actNrParameter);

    Matrix<Double> testA(3,3);
    Vector<Double> eigen(3);

   //  testA[0][0]=1.0;
//     testA[0][1]=2.0;
//     testA[0][2]=1.0;
//     testA[1][0]=2.0;
//     testA[1][1]=5.0;
//     testA[1][2]=0.5;
//     testA[2][0]=1.0;
//     testA[2][1]=0.5;
//     testA[2][2]=1.0;

//     eigenValues(testA,0.0001,eigen);
//      for (Integer i=0;i<3;i++)
//        std::cout<<eigen[i]<<std::endl;

//      std::getchar();



    for(Integer i=0;i<actNrParameter;i++)
      bas[i]=1.0;

    Vector<Complex> step(actNrParameter);

    createF(ptMaterial, ptBCs, F_hat,TRUE);

    for (i=0; i<nrMeasuredData;i++)
      y_hat_F_hat[i]=y_hat[i]-F_hat[i];

    misfit=sqrt(realA2norm(y_hat_F_hat)); 
    misfit_new=misfit;
    lin_misfit;
    Double eta = 0.9;

    std::cout << "\n misfit = |y-F| =  " <<  misfit <<std::endl;

    while((misfit>tau*delta||nrNewtonIterations<25)&&nrNewtonIterations<20){ // Newton

    for (Integer i=0;i<parameter.GetSize();i++)
      scaling[i]=1.0;

      nrNewtonIterations++;
      createF(ptMaterial, ptBCs, y_hat,FALSE);

      createJacobiMatrix2(JacobiMatrix);
      std::cout<<"\n Tichonov 1"<<std::endl;

      createAdjointJacobiMatrix(JacobiMatrix,adjJacobiMatrix);

      std::cout<<"\n Tichonov 2"<<std::endl;
      //      adjJacobiMatrix.Mult(y_hat_F_hat,basC);
      basC = y_hat_F_hat*adjJacobiMatrix;

      std::cout<<"\n Tichonov 3"<<std::endl;
      
      std::cout<<"\n update bas"<<std::endl;
      for(Integer i=0;i<actNrParameter;i++)
	bas[i]=basC[i].real();

      std::cout<<"\n Tichonov 4"<<std::endl;
      // Choice of regularisation parameter alpha!
      Integer innerIter=0;

      while (lin_misfit>eta*misfit_new&&innerIter<1||innerIter<1){
	innerIter++;

      Double reg_alpha=5.0e-15;

      Matrix<Complex> JacobiMatrixNE(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
      Matrix<Double> JacobiMatrixNE_R(JacobiMatrix.GetSizeCol(), JacobiMatrix.GetSizeCol());
      Vector<Complex> y_hat_F_hatNE(JacobiMatrix.GetSizeCol());

      adjJacobiMatrix.Mult(JacobiMatrix,JacobiMatrixNE);
      adjJacobiMatrix.Mult(y_hat_F_hat, y_hat_F_hatNE);


	std::cout<<"\n --> Value of reg_alpha = " << reg_alpha << std::endl;
	std::cout<<"\n lin_misfit = " << lin_misfit << "; misfit_new = "<<misfit_new<<std::endl;

	std::cout<<"\n Build up normal equation"<<std::endl;

	for (Integer i=0;i<JacobiMatrixNE.GetSizeRow();i++)
	  for (Integer j=0;j<JacobiMatrixNE.GetSizeCol();j++){
	    JacobiMatrixNE[i][j]=JacobiMatrixNE[i][j].real();
	    JacobiMatrixNE_R[i][j]=JacobiMatrixNE[i][j].real();
	    if(i==j)
	      JacobiMatrixNE[i][i]=JacobiMatrixNE[i][i]+reg_alpha;
	     std::cout << std::setprecision(10);
	    std::cout<<JacobiMatrixNE_R[i][j]<<", ";
	    if(j==JacobiMatrix.GetSizeCol()-1)
	      std::cout<<"; \n";
	  }

	std::cout<<"\n"<<std::endl;

	for (Integer j=0;j<y_hat_F_hatNE.GetSize();j++){
	  y_hat_F_hatNE[j]=y_hat_F_hatNE[j].real();
	  //std::cout<<y_hat_F_hatNE[j].real()<<" ";
	}

	// 	Vector<Double> eigenvalues(JacobiMatrixNE_R.GetSizeRow());
	// 	eigenValues(JacobiMatrixNE_R,0.000001,eigenvalues);

	// 	for(Integer i=0;i<eigenvalues.GetSize();i++)
	// std::cout <<" eig ("<<i<<") = "<< eigenvalues[i] << ", "<<std::endl; 

	//	std::cout<<"\n Condition number of normalequation: "<< sqrt(eigenvalues[0]/eigenvalues[eigenvalues.GetSize()-1])<<std::endl;

	//	std::getchar();
	JacobiMatrixNE.DirectSolve(s,y_hat_F_hatNE);
      
        std::cout<< "\n s determined with direct solver ...\n " <<std::endl;
	for (int i=0;i<parameter.GetSize();i++){ 
	  std::cout<<"s("<<i<<")="<<s[i]<<"; \t";
	  //	  parameter_new[i]=parameter[i]; //+(1.0/scaling[i])*step[i].real();        
	  //		parameter_new=parameter;
	  //	std::cout<<"parameter("<<i<<")="<<parameter[i]<<"; \t";
	}	     
	parameter_new = parameter;
      Double theta = 1.0;
      Matrix<Double> *matMat = ptMaterial->GetMatrix();
      scaling[0]=1.0/((*matMat)[0][0]); 
      scaling[1]=1.0/((*matMat)[2][2]);
      scaling[2]=1.0/((*matMat)[1][0]);
      scaling[3]=1.0/((*matMat)[0][2]);
      scaling[4]=1.0/((*matMat)[3][3]); 
      scaling[5]=1.0/((*matMat)[6][4]);
      scaling[6]=std::abs(1.0/((*matMat)[8][0]));
      scaling[7]=1.0/((*matMat)[8][2]);
      scaling[8]=1.0/((*matMat)[6][6]); 
      scaling[9]=1.0/((*matMat)[8][8]);
      
    for (Integer i=0;i<nrParameter;i++)
      stepR[i]=s[i].real();

      setNewParameterSet(parameter, parameter_new, scaling, theta, stepR, whichParameterToUpdate);

// 	parameter_new[9]=parameter[9]+(1.0/scaling[9])*step[9].real();        // eps33
// 	parameter_new[1]=parameter[1]+(1.0/scaling[1])*step[1].real(); // c33
// 	parameter_new[7]=parameter[7]+(1.0/scaling[7])*step[7].real(); //e33
	      
        std::cout<<"\n"<<std::endl;
	updateMaterialData(parameter_new,ptMaterial);

	createF(ptMaterial, ptBCs, F_hat,FALSE);
	y_hat_F_hat=y_hat-F_hat;
	//	misfit = realA2norm(y_hat_F_hat);
	JacobiMatrix.Mult(step,JacFS);
	linres = JacFS-y_hat_F_hat;
	lin_misfit = sqrt(realA2norm(linres));
	misfit_new = sqrt(realA2norm(y_hat_F_hat));
	reg_alpha = 1.1 * reg_alpha;
	
	for (Integer i=0;i<nrParameter;i++){
	  std::cout<<"Paramter_new["<<i<<"]= " << parameter_new[i]<<std::endl; 
	}      

      }

    
      std::cout<<"\n"<<std::endl;

      // precautionary measure, if paramters tend to far away from initial values ...
      //       for(Integer i=0;i<nrParameter;i++)
      // 	if (std::abs(parameter_new[i]-parstart[i])>0.5*std::abs(parstart[i])){
      // 	  parameter_new[i]=parameter[i];
      // 	  std::cout<<"\n parameter("<<i<<") was reset to value " << parameter[i] <<std::endl;
      // 	}

      updateMaterialData(parameter_new, ptMaterial);

      createF(ptMaterial, ptBCs, F_hat,FALSE);
      y_hat_F_hat=y_hat-F_hat;
      misfit = realA2norm(y_hat_F_hat);
      parameter=parameter_new;
	

    } // end while
  }

  // void piezoParamIdent::jacobi(Matrix<Double>& a, Double eps, Integer l_sort, Integer l_print, Vector<Double> & d)
//   {
//     std::cout<<"\nEntering jacobi ... "<<std::endl;
//     Integer i, j;
//     Integer k, kmax;
//     Integer l_conv;
//     Double a2, eps2, dkmax;
//     Double n, n2;

//     Integer ndim = a.GetSizeRow();
  
//     a2 = 0.0;
//     for (i=0; i<ndim; ++i)
//       for (j=0; j<i-1; ++j)
// 	a2 += a[i][j] * a[i][j];
//     a2 *= 2.0;
//     std::cout<<"\nEntering jacobi 1... "<<std::endl;
  
//     n = (Double)ndim;
//     n2 = n * n;
//     eps2 = eps * eps;
//     dkmax = std::log(eps)/std::log((n2-n-2)/(n2-n));
//     kmax = (Integer)std::ceil(dkmax);

//     std::cout<<"\nEntering jacobi 2... "<<std::endl;
  
//     for (k=1; k<kmax; ++k)
//       {
// 	givens_rotation(ndim, a);
// 	test_termination(ndim, a, a2, eps2, &l_conv);

// 	if (l_conv==1)
// 	   break;
//         }
//     if (l_conv==0) std::cerr("\n Note: The Jacobi method did not converge!.\n");
  
//     for (i=0; i<ndim; ++i)
//       d[i] = a[i][i];
//     sort_array(ndim, l_sort, d);
//   }
//   /*
//     Givens rotation (Rutishauser's rule)
  
//     (input)
//     ndim int : the matrix size
//     a    double [][NDIM] : the square matrix A
//   */
//   void piezoParamIdent::givens_rotation(Integer ndim, Matrix<Double> & a)
//   {
//     std::cout<<"\nEntering givens rotation ... "<<std::endl;
//     Matrix<Double> b(ndim,ndim);
//     //Double b[NDIM][NDIM];
//     Double a2, max_a2;
//     Integer i, j, k, p, q;
//     Double z, t, c, s, u;
  
//     max_a2 = 0.0;
//     for (i=0; i<ndim; ++i)
//       {
// 	for (j=0; j<=i-1; ++j)
// 	  {
// 	    a2 = a[i][j];
// 	    a2 = a2 * a2;
// 	    if (a2 > max_a2)
// 	      {
// 		p = i; 
// 		q = j;
// 		max_a2 = a2;
// 	      }
// 	  }
//       }
//     std::cout<<"\nEntering givens rotation 1... "<<std::endl;
  
//     z = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
//     t = std::fabs(z) + std::sqrt(1.0 + z*z);
//     if (z < 0.0) t = - t;
//     t = 1.0 / t;
//     c = 1.0 / std::sqrt(1.0 + t*t);
//     s = c * t;
//     u = s / (1.0 + c);
//     std::cout<<"\nEntering givens rotation 2... "<<std::endl;
  
//     for (i=0; i<ndim; ++i)
//       for (j=0; j<ndim; ++j)
// 	b[i][j] = a[i][j];
  
//     b[p][p] = a[p][p] - t * a[p][q];
//     b[q][q] = a[q][q] + t * a[p][q];
//     b[p][q] = b[q][p] = 0.0;
//     for (j=0; j<ndim; ++j)
//       if ((j!=p) && (j!=q))
// 	{
// 	  b[p][j] = b[j][p] = a[p][j] - s * (a[q][j] + u*a[p][j]);
// 	  b[q][j] = b[j][q] = a[q][j] + s * (a[p][j] - u*a[q][j]);
// 	}
  
//     for (i=0; i<ndim; ++i)
//       for (j=0; j<ndim; ++j){
// 	a[i][j] = b[i][j];
// 	std::cout<<a[i][j]<<", ";

//       }
//     std::cout<<"\n Leaving givens rotation ... "<<std::endl;
  
//   }
//   /*
//     termination criterion
  
//     The iteration is terminated if F(k) <= eps^2 * F(0), 
//     where "eps" is the relative error tolerance, 
//     F(k) is the sum of all the non-diagonal elements squared of the matrix 
//     after the k-th iteration (Given's rotation).

//     (input) 
//     ndim    int :    the matrix size
//     a       double [][NDIM] : the matrix after each Givens rotation. 
//     a2      double : the sum of all the non-diagonal elements squared 
//     of the matrix A.
//     eps2    double : the square of the error tolerance "eps", i.e., eps^2.

//     (output) 
//     l_conv int * : l_conv = 1 if converged and 
//     l_conv = 1 if not yet converged. 
//   */
//   void piezoParamIdent::test_termination(Integer ndim, Matrix<Double> & a, Double a2, Double eps2, Integer *l_conv)
//   {
//     std::cout<<"\nEntering test_termin ... "<<std::endl;
//     Double a_nd2;
//     Integer i, j;
  
//     *l_conv = 0;
  
//     a_nd2 = 0.0;
//     for (i=0; i<ndim; ++i)
//       for (j=0; j<i-1; ++j)
// 	a_nd2 += a[i][j] * a[i][j];
//     a_nd2 *= 2.0;
  
//     if (a_nd2/a2 < eps2) *l_conv = 1;
// std::cout<<"\n leaving test_termin ... "<<std::endl;
//   }
//   /*
//     sorting of a 1-dimensional array d(1), d(2), ..., d(n)
//   */
//   void piezoParamIdent::sort_array(Integer ndim, Integer l_sort, Vector<Double> & d)
//   {
//     Double dv;
//     Integer k, i;

 
//     if (l_sort == 0)
//       {
// 	for (k=0; k<ndim-1; ++k)
// 	  for (i=k+1; i<ndim; ++i)
// 	    if (d[i] > d[k])
// 	      {
// 		dv = d[k];
// 		d[k] = d[i];
// 		d[i] = dv;
// 	      }
//       }
//     if (l_sort == 1)
//       {
// 	for (k=0; k<ndim-1; ++k)
// 	  for (i=k+1; i<ndim; ++i)
// 	    if (d[i] < d[k])
// 	      {
// 		dv = d[k];
// 		d[k] = d[i];
// 		d[i] = dv;
// 	      }
//       }
//   }
}// end namepsace
