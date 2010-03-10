// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PIEZO_PARAM_IDENT
#define FILE_PIEZO_PARAM_IDENT

#include <def_use_lapack.hh>

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

namespace CoupledField
{
// forward class declaration
class SinglePDE;
class StdPDE;
class SingleDriver;
class PiezoPDE;
class PiezoCoupling;
class DirectCoupledPDE;


  //! Driver class for an inverse problem:
  //! The identification of material parameters in a piezoelectric body.
  class piezoParamIdent : public SingleDriver
  {

  public:

    //! constructor
    //! \param driverTag tag for current driver
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    piezoParamIdent( UInt sequenceStep,
                     bool isPartOfSequence = false);

    //! Destructor
    ~piezoParamIdent();

    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) {
      WARN("Tom Lahmer: Does there exist a meaningful value for this?");
      return 0;
    }

    //! Initialization method
    void Init();

    //! Starts parameter identification
    void SolveProblem(bool write_results = true, InfoNode* given_analysis_id = NULL, const bool reAssembleMatrices = true);

  protected:
    //! Calculates the parameter to soution map F(p^k) at Newton iteration step k
    //! \param F_hat - contains calculated charge, each entry belongs to different frequency
    void createF(Vector<Complex> & F_hat, bool typeOut);

    //! like create F but the frequencies are specified with vector frequencies
    void createFVec(Complex & F_hat, bool typeOut,
                    Double frequency);

    //! Write all Tensors(cE, sE, ...) in File
    void writeTensorsInFile();

    //! inverts a Matrix
    void invert(Matrix<Complex> & data);

#ifdef USE_LAPACK
    //! inverts a Matrix with Lapacks ZGESV
    void invertWithLapack(Matrix<Complex> & data);

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

    //! Creates sensitivity vector, needed for computation of confidence intervals
    void computeVecOfSingleDeriv(Vector<Complex> & jacobi, Double omega);

    //! Calculates gradient of minimization problem \nabla J(w) in opt.Exp.design
  //  void createGradient(Vector<Complex> & grad, Double dOmega);

    //! Calculates explicitely the Adjoint operator of F'
    void createAdjointJacobiMatrix();

    //! Method which reads data from file measuredData.dat
    //! which includes frequencies where to fit and initial guesses
    void readMeasuredData();

    //! reads  measured data
    void readInMeasurement(Vector<Double> & frequencies);

    //! updates the piezoMatrix in MaterialData parameter =
    //! \f$(c_11, c_33, c_12, c_13, c_44, e_15, e_31, e_33, eps_11, eps_33)\f$
    //! \param parameter - new set of piezoelectric material parameters
    void updateMaterialData(Vector<Double> & parameter);

    //! updates the imaginary parts
    void updateComplexMaterialData(Vector<Double> & parameterC);

    //! computes vector scaling_ which consists of reciprocals of all parameters
    void computeScaling();

    //! overwrites values in paramter_new with paramter+step if whichParamterToUpdate ==1
    void setNewParameterSet(Vector<Double> & parameter,Vector<Double> &  parameter_new,Vector<Double> & scaling,Double & theta,
                            Vector<Double> & step, Vector<UInt> & whichParameterToUpdate);

    //! Calculates and writes out impedance curve into file imped.dat
    //! Start - and stopfrequency will be specified in xml - file
    void calcImpedanceCurve();

    //! square root of sum of all squared entires in matrix
    Double calcEuclidianMatrixNorm(Matrix<Complex> & mat);

    //! nu - methods, semiiterative methods with optimal rate of convergence
    void nuMethods();

     //! nu - methods, complex version - uses weighted norms
    void nuMethodsC2();

    //! Solve problem with a least square error minimization
    void leastSquare();

    //! nonlinear Landweber's iteration p^{k+1} = p^k - \omega F'*(p^k)(F(p^k)-y)
    void nonlinLandweber();

    //! saves sysmat of forward problem, multiplication with \omega*\beta*j ...
    void createAndSetRHSforJacobian(UInt & fstep);

    //! calculates charges out of measurements of |Z|, phase and voltage for different frequencies
    void evaluateMeasuredData(Vector<Double> freqs, Vector<Double> & absZ,
                             Vector<Double> & phi, Vector<Complex> & y_hat);

    //! calculates Euclidian vector norm
    Double a2norm(Vector<Double> &vec);

    //! calculates Euclidian vector norm
    Double a2norm(Vector<Complex> &vec);

        //! function in which user specified norms (file:measuredData.dat) will be calculated
    void norm(Vector<Complex> &  vec, Double & norm, Double & norm2,Vector<Complex> & q_meas);


    void maxAndWeightedResNorm(Vector<Complex> & vec, Double & maxNorm, Double & wNorm, Vector<Complex> & q_meas);


    //! approximates Jacobian by a second order FDM
    void computeJacobiMatrix();

    //! approximates Jacobian by a second order FDM for real and complex values
    void computeJacobiMatrixC();

    //! This method implements the different sampling variants for the
    //! frequency domain. We currently support linear, logarithmic and
    //! reverse logarithmic sampling.
    //! \param freqIndex index of the frequency that shall be computed, i.e.
    //!                  an integral value from [1:numFreq_]
    Double ComputeNextFrequency( UInt freqIndex ) const;


#ifdef USE_LAPACK
    //! Fortran f77 variables
    F77complex16 *lp_af77;
    F77real8 * lp_wf77;
    F77complex16 *lp_workf77;
    F77real8 *lp_rworkf77;
#endif


  private:


    //! input file, here problem specific details will be set
      std::ifstream * allMeasuredData_;
      // ! input of impedance measurements
      std::ifstream inMess_;
      // ! input of mechanical measurements
      std::ifstream inMechMess_;
      //! output file, simulated impedance curve
      std::ofstream * impedCurve_;
      //! output file, norm of residue and parameters
      std::ofstream * parLog_;
      //! output file, finally calculated parameters
      std::ofstream * parFinal_;
      //! output file, mechanical displacement at spec.node and dof
      std::ofstream * mechDispl_;
      //! output file, contains parameters and their confidence intervals
      std::ofstream * confInterval_;
      //! output file, function rho(j) in optimal exp. design
      //std::ofstream * rhosOut;
      //! output file, writes synthetically created impedance curve
      std::ofstream * synMess_;
      //! output file, writes synthetically created impedance curve
      //std::ofstream * nrOfFreqs;
      // output file, contains all Tensors
      std::ofstream * allTensors_;

    //! Parameter node for "paramIden"-element from xml-parameter file
        ParamNode * myParam_;

        bool CalcImpedanceCurve_;
        bool CalcMechDisplCurve_;
        bool directCoupling_;
        bool isPartOfSequence_;
        bool imagMaterialParam_;
        bool writeResults_;

        UInt maxNumberInnerLoops_;
        UInt maxNumberNewtonLoops_;
        UInt mechDisplAtNode_;
        UInt dofOfMechDispl_;
        UInt computeImpedanceCurveAfterStep_;
        UInt globalIterationNr_;

        std::string whichMethod_;

        // Parameter vectors
        // contains all real parts
        Vector<Double> parameter_;
        // contains all imaginary parts
        Vector<Double> parameterC_;

        // contains new iterates
        Vector<Double> parameter_new_;
        Vector<Double> parameter_newC_;

        // parameter increments
        Vector<Double> parameterIncrement_;

        // stores initial guess
        Vector<Double> parameterInitial_;

        Vector<UInt> whichParameterToUpdate_;
        Vector<UInt> whichParameterToUpdateC_;
        Vector<UInt> whichParToUpInd_;
        Vector<UInt> whichParToUpIndC_;

        UInt nrParameter_;
        UInt actNrParameter_;
        UInt actNrParameterC_;

        StdVector<RegionIdType> subdomsMech_;
        StdVector<RegionIdType> subdomsElec_;

        // nu_ is the steering parameter in inner iteration of Newtons-nu Method
        Double nu_;

        // parameter update in Newtons Method, ie.e
        // solution of the linearized problem
        Vector<Complex> s_;

        Double startfreq_;
        Double stopfreq_;
        Double startFreq_;
        Double stopFreq_;
        UInt nrfreq_;
        FreqSamplingType samplingType_;
        UInt newtonCounter_;
        Double stopRes_;

        // confidence interval computation:
        Double residuumParIdent_;
        Matrix<Complex> globalCov_;
        Vector<Complex> globalJacobi_;
        Vector<Complex> globalJacobiH_;

        //Vector<Complex> solElecPot_;
        Vector<Complex> solMechDispl_;
        //   Vector<Complex> algSysSolVector_;

        // Pointers to other classes of cfs
        BaseSystem * ptAlgsys_;
        Assemble * ptAssemble_;
        Assemble * ptAssemble2_;

        std::map<RegionIdType, BaseMaterial*> ptMaterialMech_;
        std::map<RegionIdType, BaseMaterial*> ptMaterialElec_;
        std::map<RegionIdType, BaseMaterial*> ptMaterialPiezo_;

        BaseMaterial * ptMaterial1_;
        BaseMaterial * ptMaterial2_;

        Vector<Double> freqs_;

        // Vectors containing measured data
        Vector<Double> real_, imag_;
        Vector<Double> realMech_, imagMech_;

        UInt numMechMeasurements_;
        UInt nrMeasuredData_;

        // contains right hand side, i.e. we solve F_hat_(p)=y_hat
        Vector<Complex> y_hat_;
        // solution of forward problem
        Vector<Complex> F_hat_;


        Vector<Double> scaling_;
        Vector<Double> scalingC_;

        // Numerical linearizations of forward soution operator
        Matrix<Complex> JacobiMatrix_;
        Matrix<Complex> adjJacobiMatrix_;
        Matrix<Complex> approxJacobiMatrix_;

        Double voltage_, delta_, anorm_;

        std::string whichNormCriteria_;
        std::string whichNormCriteriaMech_;

        Double resonanceFrequency_;
        Double antiResonanceFrequency_;

        // Pointer to Single PDE
        SinglePDE * ptMyPDE_;

        //! Pointer to mechanic pde while direct-coulped
        SinglePDE * ptPDE1_;

        //! Pointer to electrostatic pde while direct-coupled
        SinglePDE * ptPDE2_;

        //! identification tag of PDE for algebraic system
        PdeIdType pdeId_;

        //! pointer to piezo coupling object
        PiezoCoupling * piezoCpl_;

        //! result object containing charges
        shared_ptr<Result<Complex> > charges_;

        //! neighboring region for charge computation
        RegionIdType chargeNeighborRegion_;


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
