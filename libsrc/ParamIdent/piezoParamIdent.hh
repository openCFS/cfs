#ifndef FILE_PIEZO_PARAM_IDENT
#define FILE_PIEZO_PARAM_IDENT

#include "Driver/singleDriver.hh"
#include "Driver/assemble.hh"


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


// forward class declaration
class SinglePDE;
class StdPDE;
class SingleDriver;
class PiezoPDE;


namespace CoupledField 
{

  

  //! Driver class for an inverse problem: The identification of material parameters in a piezoelectric body.
  class piezoParamIdent : public SingleDriver
  {

  public:

    //! constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time
    //! \param driverTag tag for current driver
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    piezoParamIdent(Domain * adomain,
                    Integer stepOffset = 0,
                    Double timeOffset = 0.0,
                    std::string driverTag = "anyTag",
                    Boolean isPartOfSequence = FALSE);

    //! Destructor
    ~piezoParamIdent();

    std::ifstream * allMeasuredData; // Contains selected measurements & further steering parameters
    std::ifstream * mess; // contains whole set of measurements
    std::ofstream * impedCurve;
    std::ofstream * piezoLog;
    std::ofstream * parLog;
    std::ofstream * parFinal;
    std::ofstream * mechDispl;
    std::ofstream * optimalFreqs;
    std::ofstream * confInterval;
    std::ofstream * rhosOut;
    std::ofstream * synMess;
    //  std::ofstream impedCurve("impedCurve.dat");
    //  std::ofstream impedCurve("impedCurve.dat");

    //! Starts parameter identification
    void SolveProblem();



  protected:
    //! Calculates the parameter to soution map F(p^k) at Newton iteration step k
    //! \param F_hat - contains calculated charge, each entry belongs to different frequency 
    void createF(MaterialData * ptMaterial, Vector<Complex> & F_hat, Boolean typeOut);

    // ! like create F but the frequencies are specified with Vector frequencies
    void createFVec(Complex & F_hat, Boolean typeOut,
                 Double frequency);

    // ! inverts a Matrix
    void invert(Matrix<Complex> & data);

    // ! inverts a Matrix with Lapacks ZGESV
    void invertWithLapack(Matrix<Complex> & data);
    
    // Descent Method for functional J(w)
    void descentMethod(Complex & functional);

    //  Descent Method for functional J(rho) with constraints
    void descentMethodRho(Complex & functional);

    // gradient of J(\rho)
    void createGradientRho(Vector<Complex> & grad, Double dOmega);


#ifdef USE_LAPACK
    //! Converts a fortran 77 matrix to C++ complex
    void F772CC( const F77complex16 &v, std::complex<double> &val ) {
      std::complex<double> aux(v.real,v.imag);
      val = aux;
    }

    void F772CC( const F77real8 &v, double &val ) {
      val = (double)v;
    }

    //! Converts cfs data to fortran 77 format
    void CC2F77( const std::complex<double> &v, F77complex16 &val ) {
      val.real = (F77real8)v.real();
      val.imag = (F77real8)v.imag();
    }
    void CC2F77( const double &v, F77real8 &val ) {
      val = (F77real8)v;
    }
#endif

    //! Calculates an approximation of the Jacobi Matrix of parameter to solution operator F 
    //! \param Jacobi Matrix - approximation of F'
  //   void createJacobiMatrix(MaterialData * ptMaterial, Vector<Complex> & F_hat, Vector<Double> & parameterIncrement,
//                             Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl);

    //! Creates special Jacobi Matrix, for optimal experiment design
    void createJacobian(Vector<Complex> & jacobi, Double omega);

    //! Calculates gradient of minimization problem \nabla J(w) in opt.Exp.design
    void createGradient(Vector<Complex> & grad, Double dOmega);

    //! Determines variance - covariance Matrix 
    void createCovA(Complex &J, Boolean writeOutCov);

    //! Determines J(omega) in case of flexible number of frequencies
    void createJRho(Complex &J, Boolean writeOutCov);

    //! Not in use in this version
//     void createJacobiMatrix2(Matrix<Complex> & JacobiMatrix);
//     void createJacobiMatrixC(Matrix<Complex> & JacobiMatrix);

    //! Calculates explicitely the Adjoint operator of F'
    void createAdjointJacobiMatrix(Matrix<Complex> & JacobiMatrix, Matrix<Complex> & adjJacobiMatrix);

    //! Method which reads Data from file measuredData.dat. The file contains measurements of amplitude, frequency, further 
    //according information concerning the piezhoelectric body (radius, thickness, ...)
    void readMeasuredData(Vector<Double> & freqs, Vector<Double> & real, Vector<Double> & imag ,Vector<Double> & parameter, 
                          Double & voltage, UInt & nrMeasuredData, Double & thickness, Double & radius, Double & delta);

    //! reads whole set of measured data
    void readInMeasurement(Vector<Double> & frequencies);

    //! updates the piezoMatrix in MaterialData parameter = 
    //! \f$(c_11, c_33, c_12, c_13, c_44, e_15, e_31, e_33, eps_11, eps_33)\f$
    //! \param parameter - new set of piezoelectric material parameters
    void updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial);

    void updateComplexMaterialData(Vector<Double> & parameterC, MaterialData * ptMaterial);

    //! overwrites values in paramter_new with paramter+step if whichParamterToUpdate ==1
    void setNewParameterSet(Vector<Double> & parameter,Vector<Double> &  parameter_new,Vector<Double> & scaling,Double & theta,
                            Vector<Double> & step, Vector<UInt> & whichParameterToUpdate);

    //! Calculates the impedance curve of piezo-simulation, writes results to file imped.dat
    void calcAbsImped(Complex & charge, Double & freq, UInt & fstep, Boolean typeOut);

    void calcImpedanceCurve();

    void calcMechDisplCurve();

    void updateRHS(Vector<Complex> & solElecPot, Vector<Complex> & mechDisplacement, Double omega);

    void updateRHS(Vector<Complex> & RHSsol);

    void updateRHS2(Vector<Complex> & RHSsol);

    //! types out nodal results of elecPot and mechanical displacement
    void typeOutSolutionOnScreen(Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl);

    void calcInitialResidual(Vector<Complex> & res, Vector<Complex> & y_hat, Vector<Complex> & PHI_p, UInt fstep,
                             Vector<Complex> & solElecPot, Double & meanValueMechDeformation);

    void measureMechDeformationInZ_Direction(Vector<Complex> & mechDisplacement, Double & Radius, 
                                             Double & meanValueMechDeformation, UInt dof);
   
    void calcNorm2Resid(Vector<Complex> &res, Double & anorm, UInt nrMeasuredData);

    Double calcEuclidianMatrixNorm(Matrix<Complex> & mat);
    
  
    //! see SFBReport F013 for details ;-)
  //   void NewtonCG();

//     //! The version from 10.12.04 ...
//     void NewtonCG1();

//     //! see SFBReport F013 for details ;-)
//     void NewtonCG2();

    //! see SFBReport F013 for details ;-)
    void NewtonCG3();

    //! The version schich is supposed to treat complex-valued material Parameter
    void NewtonCG4();

    //! Iterative Method to determine material parameter
    void NewtonLandweber();

    //! Iterative Method to determine complex valued material parameter
    void NewtonLandweberC();


    // ! nu - methods, semiiterative methods with optimal rate of convergence
    void nuMethods();

    // ! nu - methods, semiiterative method, complex version
    void nuMethodsC();

    // ! nu - methods, complex version - uses weighted norms
    void nuMethodsC2();

    //! methods which determines a set of parameters for an optimal experiment design
    void optimalExpDesign();

    // ! To be removed ... is now a kind of multilevel algo ...
    void optimalExpDesignDiffNumberFreqs();

    //! methods which determines a set of parameters for an optimal experiment design
    //! with flexible number of frequencies ...
    void optimalExpDesignVarNrFreqs();

    //! Try to solve problem with a least square approach
    void leastSquare();

    //! saves sysmat of forward problem, multiplication with \omega*\beta*j ...
    void createAndSetRHSforJacobian(UInt & fstep);

    //! calculates charges out of measurements of |Z|, phase and voltage for different frequencies
    void calc_measuredCharge(Vector<Double> freqs, Vector<Double> & absZ,
                             Vector<Double> & phi, Vector<Complex> & y_hat);
    
    //! calculates Euclidian vector norm
    Double a2norm(Vector<Double> &vec);

    //! calculates Euclidian vector norm
    Double a2norm(Vector<Complex> &vec);
    
    //! Calculates Euclidian norm of only real-parts of vec
    Double realA2norm(Vector<Complex> &vec);

    //! Calculates Euclidian norm of only real-parts of vec
    Double norm2Real(Vector<Complex> &vec);

    void norm(Vector<Complex> &  vec, Double & norm, Double & norm2,Vector<Complex> & q_meas);

    void maxAndEuclNorm(Vector<Complex> & vec, Double & maxNorm, Double & euclNorm);

    void logNorm(Vector<Complex> & vec, Double & logNorm);
  
    void maxAndWeightedResNorm(Vector<Complex> & vec, Double & maxNorm, Double & wNorm, Vector<Complex> & q_meas);

    //! Performs a forward simulation with exact data, adds to results alternating +- 10 Percent
    void calcSyntheticData(Vector<Complex> & y_hat);

    void createMaterialTensorMatrices(Vector<Double> & parameter, Matrix<Double> & couplingMatrix,
                                      Matrix<Double> & dielectricMatrix, UInt spaceDim);

    //! Tests, if JacobiMatrix is more or less approximated by F(p)-F(p+delta)/delta
    void testJacobiMatrix(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,
                          MaterialData * ptMaterial, Vector<Double> & parameterIncrement, 
                          Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl);

    void testJacobiMatrix2(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,
                           MaterialData * ptMaterial, Vector<Double> & parameterIncrement, 
                           Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl);


    void testJacobiMatrixC(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,
                           MaterialData * ptMaterial);
    // ! The following methods serve for the determination of eigenvalues ...

    //    void sort_array(UInt ndim, UInt l_sort, Vector<Double> & d);
    //     void test_termination(UInt ndim, Matrix<Double> & a, Double a2, Double eps2, UInt *l_conv);
    //     void givens_rotation(UInt ndim, Matrix<Double> & a);
    //     void jacobi(Matrix<Double>& a, Double eps, UInt l_sort, UInt l_print, Vector<Double> & d);

#ifdef USE_LAPACK
    //! Fortran f77 variables
    F77complex16 *lp_af77;
    F77real8 * lp_wf77;
    F77complex16 *lp_workf77;
    F77real8 *lp_rworkf77;
#endif

    Double fa_, fr_;
    Boolean CalcImpedanceCurve;
    Boolean CalcMechDisplCurve;
    Boolean directCoupling;
    UInt whichNewtonCG;
    UInt maxNumberInnerLoops;
    UInt maxNumberNewtonLoops;
    UInt mechDisplAtNode;
    UInt dofOfMechDispl;
    UInt whichNorm;
    Boolean considerMechDeformation;
    Vector<UInt> whichParameterToUpdate;
    Vector<UInt> whichParameterToUpdateC;
    Vector<UInt> whichParameterToUpdateRC;
    Vector<UInt> whichParToUpInd;
    Vector<UInt> whichParToUpIndC;
    Double sign;
    Matrix<Double> * piezoMatrix;
    UInt dofs;
    UInt numNodes;
    BaseSystem * ptAlgsys;
    Grid * ptGrid;
    NodeEQN * ptNodeEqn;
    Assemble * ptAssemble;
    StdVector<RegionIdType> subdoms;
    BaseNodeStoreSol * sol;
    Domain * ptDomain;
    Double startfreq;
    Double stopfreq;
    UInt nrfreq;
    Double finalnorm;
    UInt newtonCounter;
    Double inner_eta;


    // optimal experiment design:
    Double residuumParIdent;
    Double residuumParIdentOld;
    Vector<Double> projGradientFlags;
    Double normGradient;
    Double normGradientOld;
    Vector<Double> omegaDiffVec;


    Vector<Complex> solElecPot;
    Vector<Complex> solMechDispl;
    Vector<Complex> algSysSolVector;
    MaterialData * ptMaterial;

    Vector<Double> parameter;
    Vector<Double> parameterC;
    Vector<Double> parameter_new;
    Vector<Double> parameter_newC;
    Vector<Double> parameterIncrement;
    Vector<Double> parameterInitial;
    Vector<Double> omegas;
    Vector<Double> freqs;
    Vector<Double> real, imag;
    Vector<Complex> amplitude_phase;
    Vector<Complex> F_hat;
    Vector<Complex> y_hat;
    Vector<Complex> s_0;    
    Vector<Double> bas;
    Vector<Complex> basC;
    Vector<Complex> res_NE, res_NE_new, lin_res;
    Vector<Complex> res;
    Vector<Complex> bas_bar;
    Vector<Complex> s;
    Vector<Double> scaling;
    Vector<Double> scalingC;
    Double alpha_m, beta_m;
    Matrix<Complex> JacobiMatrix;
    Matrix<Complex> adjJacobiMatrix;
    Matrix<Complex> approxJacobiMatrix;
    Double voltage, thickness, radius, delta, meanValueMechDeformation, anorm, tau;
    UInt nrMeasuredData;
    UInt nrParameter;
    UInt actNrParameter;
    UInt actNrParameterC;
    Double relaxParameter;
    UInt nrSol, sizeSol;
    Double norm_res;
    Double eta;
    Double overall_res;
    Vector<Complex> overall_res0;
    Vector<Double> rhos; // rhos in optimalExpDesignVarNrFreqs

    Matrix<Complex> completeSolOf_F;
    Matrix<Complex> allElemsVec;

    //! pointer to SinglePDE

    //! it is the same as ptPDE, only that it is of the type
    //! SinglePDE *, since piezoParamIdebnt makes only sense
    //! with a PiezoPDE, which is of the latter type.
    SinglePDE * ptMyPDE_;

    //! Pointer to mechanic pde while direct-coupled
    SinglePDE * ptPDE1_; 

    //! Pointer to electrostatic pde while direct-coupled
    SinglePDE * ptPDE2_; 

    //    StdVector<DirectCoupledPDE*> ptDirectCoupledPde_;

    //! identification tag of PDE for algebraic system
    PdeIdType pdeId_;
  private:



  }; // end of class piezoParamIdent

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class piezoParamIdent
  //! 
  //! \purpose Determination of piezoelectric material parameters.
  //! Class implements a simulation based reconstruction principle.
  //! 
  //! \collab 
  //! 
  //! \implement Code in File piezoParamIdent provides setting for
  //! parameter identification. piezoParamIdent_methods is a container 
  //! for different subroutines, like the harmonic forward solution of 
  //! the piezoelectric PDEs, the calculation of the impedance curve ... It 
  //! further comprises different realizations (for real and comlex-valued
  //! parameters of the computation of
  //! the Jacobian of the so called 'parameter-to-solution' mapping.
  //!
  //! Files like nuMethods, newtonLandweber, newtonCG provide different
  //! regularization strategies which solve the nonlinear operator equation
  //! F(p)=y, where symbolically F denotes the forward solution (charges) 
  //! of the piezoelectric 
  //! equations for different frequencies using the parameter set p. The right 
  //! hand side contains a set of measured charges at the surface of the
  //! piezo.
  //! 
  //! \status In use, still in improvement phase
  //! 
  //! \unused 
  //! 
  //! \improve Clean code from not approved methods and ansatzes.
  //! Expand functionality to determination of parameters on different
  //! domains and of Lame parameters in non-piezoelectric clamping material
  //! 

#endif
}



#endif
