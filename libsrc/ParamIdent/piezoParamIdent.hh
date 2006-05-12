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
class PiezoCoupling;
class DirectCoupledPDE;


namespace CoupledField 
{

  

  //! Driver class for an inverse problem: 
  //! The identification of material parameters in a piezoelectric body.
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

    //! input file, here problem specific details will be set
    std::ifstream * allMeasuredData; 

    //! output file, simulated impedance curve
    std::ofstream * impedCurve;
    //! output file, norm of residue, 
    std::ofstream * piezoLog;
    //! output file, norm of residue and parameters
    std::ofstream * parLog;
    //! output file, finally calculated parameters
    std::ofstream * parFinal;
    //! output file, mechanical displacement at spec.node and dof
    std::ofstream * mechDispl;
    //! output file, set  of frequencies for optimal exp. design
    std::ofstream * optimalFreqs;
    //! output file, contains parameters and their confidence intervals
    std::ofstream * confInterval;
    //! output file, function rho(j) in optimal exp. design 
    std::ofstream * rhosOut;
    //! output file, writes synthetically created impedance curve
    std::ofstream * synMess;
    //! output file, writes synthetically created impedance curve
    std::ofstream * nrOfFreqs;

    //! Starts parameter identification
    void SolveProblem();

  protected:
    //! Calculates the parameter to soution map F(p^k) at Newton iteration step k
    //! \param F_hat - contains calculated charge, each entry belongs to different frequency 
    void createF(Vector<Complex> & F_hat, Boolean typeOut);

    //! like create F but the frequencies are specified with Vector frequencies
    void createFVec(Complex & F_hat, Boolean typeOut,
                    Double frequency);

    //! inverts a Matrix
    void invert(Matrix<Complex> & data);

#ifdef USE_LAPACK
    //! inverts a Matrix with Lapacks ZGESV
    void invertWithLapack(Matrix<Complex> & data);
#endif
    
    //! Descent Method for functional J(w)
    void descentMethod(Complex & functional);

    //!  Descent Method for functional J(rho) with constraints
    void descentMethodRho(Complex & functional);

    //! gradient of J(\rho)
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
    // Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl);

    //! Creates special Jacobi Matrix, for optimal experiment design
    void createJacobian(Vector<Complex> & jacobi, Double omega);

    //! Calculates gradient of minimization problem \nabla J(w) in opt.Exp.design
    void createGradient(Vector<Complex> & grad, Double dOmega);

    //! Determines variance - covariance Matrix 
    void createCovA(Complex &J, Boolean writeOutCov);

    //! writes Out Sum (cov(i)(i)) for frequencies selected by
    //! optimal experiment design with var. number of freqs 
    void writeOutConfInterval();

    //! Determines J(omega) in case of flexible number of frequencies
    void createJRho(Complex &J, Boolean writeOutCov);

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
    void updateMaterialData(Vector<Double> & parameter);
    
    void updateComplexMaterialData(Vector<Double> & parameterC);

    //! writes vector scaling_ which consists of reciprocals of all parameters
    void computeScaling();

    //! overwrites values in paramter_new with paramter+step if whichParamterToUpdate ==1
    void setNewParameterSet(Vector<Double> & parameter,Vector<Double> &  parameter_new,Vector<Double> & scaling,Double & theta,
                            Vector<Double> & step, Vector<UInt> & whichParameterToUpdate);

    //! Calculates the impedance curve of piezo-simulation, writes results to file imped.dat
    void calcAbsImped(Complex & charge, Double & freq, UInt & fstep, Boolean typeOut);

    //! Calculates and writes out impedance curve into file imped.dat
    //! Start - and stopfrequency will be specified in xml - file
    void calcImpedanceCurve();

    // calculates mech. Deformation at one selected node and dof and writes it into file mech.dat
    void calcMechDisplCurve();
    
    //! calculates l2 norm of residue   
    void calcNorm2Resid(Vector<Complex> &res, Double & anorm, UInt nrMeasuredData);

    //! square root of sum of all squared entires in matrix
    Double calcEuclidianMatrixNorm(Matrix<Complex> & mat);
    
    //! nu - methods, semiiterative methods with optimal rate of convergence
    void nuMethods();

     //! nu - methods, complex version - uses weighted norms
    void nuMethodsC2();

    //! methods which determines a set of parameters for an optimal experiment design
    void optimalExpDesign();

    //! To be removed ... is now a kind of multilevel algo ...
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

    //! function in which user specified norms (file:measuredData.dat) will be calculated
    void norm(Vector<Complex> &  vec, Double & norm, Double & norm2,Vector<Complex> & q_meas);

    void maxAndEuclNorm(Vector<Complex> & vec, Double & maxNorm, Double & euclNorm);

    //! calculates norm of logarithmic value of residue
    void logNorm(Vector<Complex> & vec, Double & logNorm);
  
    void maxAndWeightedResNorm(Vector<Complex> & vec, Double & maxNorm, Double & wNorm, Vector<Complex> & q_meas);

    //! Performs a forward simulation with exact data, adds to results alternating +- 10 Percent
    void calcSyntheticData(Vector<Complex> & y_hat);

    void createMaterialTensorMatrices(Vector<Double> & parameter, Matrix<Double> & couplingMatrix,
                                      Matrix<Double> & dielectricMatrix, UInt spaceDim);

    //! Tests, if JacobiMatrix is more or less approximated by F(p)-F(p+delta)/delta
    void testJacobiMatrix(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,
                          Vector<Double> & parameterIncrement, 
                          Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl);

    //! approximates Jacobian by a second order FDM
    void testJacobiMatrix2(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,
                           Vector<Double> & parameterIncrement, 
                           Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl);

    //! approximates Jacobian by a second order FDM for real and complex values
    void testJacobiMatrixC(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter);
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
    Boolean CalcImpedanceCurve_;
    Boolean CalcMechDisplCurve_;
    Boolean directCoupling_;
    UInt whichNewtonCG_;
    UInt maxNumberInnerLoops_;
    UInt maxNumberNewtonLoops_;
    UInt mechDisplAtNode_;
    UInt dofOfMechDispl_;

    UInt whichNorm_;
    Boolean considerMechDeformation_;
    Vector<UInt> whichParameterToUpdate_;
    Vector<UInt> whichParameterToUpdateC_;
    Vector<UInt> whichParToUpInd_;
    Vector<UInt> whichParToUpIndC_;
    Double sign_;

    BaseSystem * ptAlgsys_;
    Assemble * ptAssemble_;
    Assemble * ptAssemble2_;

    StdVector<RegionIdType> subdomsMech_;
    StdVector<RegionIdType> subdomsElec_;
    Domain * ptDomain_;
    Double startfreq_;
    Double stopfreq_;
    UInt nrfreq_;
    UInt newtonCounter_;
    Double inner_eta_;

    // optimal experiment design:
    Double residuumParIdent_;
    Double residuumParIdentOld_;
    Vector<Double> projGradientFlags_;
    Double normGradient_;
    Double normGradientOld_;
    Vector<Double> omegaDiffVec_;
    Matrix<Complex> globalCov_;
    Vector<UInt> globalIndexSet_;
    Vector<Complex> globalJacobi_;
    Vector<Complex> globalJacobiH_;

    Vector<Complex> solElecPot_;
    Vector<Complex> solMechDispl_;
    //   Vector<Complex> algSysSolVector_;

    std::map<RegionIdType, BaseMaterial*> ptMaterialMech_;
    std::map<RegionIdType, BaseMaterial*> ptMaterialElec_;
    std::map<RegionIdType, BaseMaterial*> ptMaterialPiezo_;

    //    BaseMaterial * ptMaterial_;
    BaseMaterial * ptMaterial1_;
    BaseMaterial * ptMaterial2_;

    Vector<Double> parameter_;
    Vector<Double> parameterC_;

    Vector<Double> parameter_new_;
    Vector<Double> parameter_newC_;
    Vector<Double> parameterIncrement_;
    Vector<Double> parameterInitial_;
    Vector<Double> omegas_;


    Vector<Double> freqs_;
    Vector<Double> real_, imag_;
    Vector<Complex> amplitude_phase;
    Vector<Complex> F_hat_;
    Vector<Complex> y_hat_;

    Vector<Double> bas;
    Vector<Complex> basC;
    Vector<Complex> res_NE, res_NE_new, lin_res;
    Vector<Complex> res;
    Vector<Complex> bas_bar;
    Vector<Complex> s;

    Vector<Double> scaling_;
    Vector<Double> scalingC_;

    Matrix<Complex> JacobiMatrix_;
    Matrix<Complex> adjJacobiMatrix_;
    Matrix<Complex> approxJacobiMatrix_;

    Double voltage_, thickness_, radius_, delta_;
    Double anorm_, tau_;

    UInt nrMeasuredData;
    UInt nrParameter_;
    UInt actNrParameter;
    UInt actNrParameterC;
    Double relaxParameter;
    UInt nrSol, sizeSol;
    Double norm_res;
    Double eta;
    Double overall_res;
    Vector<Complex> overall_res0;
    Vector<Double> rhos; // rhos in optimalExpDesignVarNrFreqs

    Double resonanceFrequency_;
    Double antiResonanceFrequency_;

    Matrix<Complex> completeSolOf_F;
    Matrix<Complex> allElemsVec;

    //! pointer to SinglePDE

    //! it is the same as ptPDE, only that it is of the type
    //! SinglePDE *, since piezoParamIdebnt makes only sense
    //! with a PiezoPDE, which is of the latter type.
    SinglePDE * ptMyPDE_;

    //! Pointer to mechanic pde while direct-coulped
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
