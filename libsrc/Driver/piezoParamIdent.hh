#ifndef FILE_PIEZO_PARAM_IDENT
#define FILE_PIEZO_PARAM_IDENT

#include "singleDriver.hh"


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

    std::ifstream * allMeasuredData;
    std::ofstream * impedCurve;
    std::ofstream * piezoLog;
    //  std::ofstream impedCurve("impedCurve.dat");
    //  std::ofstream impedCurve("impedCurve.dat");

    //! Starts parameter identification
    void SolveProblem();

  protected:
    //! Calculates the parameter to soution map F(p^k) at Newton iteration step k
    //! \param F_hat - contains calculated charge, each entry belongs to different frequency 
    void createF(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat);

    //! Calculates an approximation of the Jacobi Matrix of parameter to solution operator F 
    //! \param Jacobi Matrix - approximation of F'
    void createJacobiMatrix(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat, Vector<Double> & parameterIncrement, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl);

    void createJacobiMatrix2(Matrix<Complex> & JacobiMatrix);

    //! Calculates explicitely the Adjoint operator of F'
    void createAdjointJacobiMatrix(Matrix<Complex> & JacobiMatrix, Matrix<Complex> & adjJacobiMatrix);

    //! Method which reads Data from file measuredData.dat. The file contains measurements of amplitude, frequency, further according information concerning the piezhoelectric body (radius, thickness, ...)
    void readMeasuredData(Vector<Double> & freqs, Vector<Double> & real, Vector<Double> & imag ,Vector<Double> & parameter, Double & voltage, Integer & nrMeasuredData, Double & thickness, Double & radius, Double & delta);

    //! updates the piezoMatrix in MaterialData parameter = \f$(c_11, c_33, c_12, c_13, c_44, e_15, e_31, e_33, eps_11, eps_33)$\f
    //! \param parameter - new set of piezoelectric material parameters
    void updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial);

    //! Calculates the impedance curve of piezo-simulation, writes results to file imped.dat
    void calcAbsImped(Complex & charge, Double & freq, Integer & fstep);

    void calcImpedanceCurve();

    void updateRHS(Vector<Complex> & solElecPot, Vector<Complex> & mechDisplacement, Double omega);

    void updateRHS(Vector<Complex> & RHSsol);

    void updateRHS2(Vector<Complex> & RHSsol);

    //! types out nodal results of elecPot and mechanical displacement
    void typeOutSolutionOnScreen(Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl);

    void calcInitialResidual(Vector<Complex> & res, Vector<Complex> & y_hat, Vector<Complex> & PHI_p, Integer fstep, Vector<Complex> & solElecPot, Double & meanValueMechDeformation);

    void measureMechDeformationInZ_Direction(Vector<Complex> & mechDisplacement, Double & Radius, Double & meanValueMechDeformation, int dof);
   
    void calcNorm2Resid(Vector<Complex> &res, Double & anorm, Integer nrMeasuredData);

    Double calcEuclidianMatrixNorm(Matrix<Complex> & mat);
    
    //! CG method, approximates F'^(-1)(F-y_hat)
    void CG();

    //! see SFBReport F013 for details ;-)
    void NewtonCG();

    //! see SFBReport F013 for details ;-)
    void NewtonCG2();

    //! see SFBReport F013 for details ;-)
    void NewtonCG3();

    //! Iterative Method to determine material parameter
    void NewtonLandweber();

    //! The classical regularisation strategy for ill-posed systems of equations
    void tichonov();

    //! Implementation of a linesearchstrategy, Idea by Eisenstat, Walker
    void backtracking(Double & eta, Double & theta, Vector<Complex> & s, Double & norm_res, Double & norm_res_new);

    //! saves sysmat of forward problem, multiplication with \omega*\beta*j ...
    void createAndSetRHSforJacobian(Integer & fstep);

    //! calculates charges out of measurements of |Z|, phase and voltage for different frequencies
    void calc_measuredCharge(Vector<Double> freqs, Vector<Double> & absZ, Vector<Double> & phi, Vector<Complex> & y_hat);
    
    //! calculates Euclidian vector norm
    Double a2norm(Vector<Double> &vec);

    //! calculates Euclidian vector norm
    Double a2norm(Vector<Complex> &vec);
    
    //! Calculates Euclidian norm of only real-parts of vec
    Double realA2norm(Vector<Complex> &vec);

    //! Calculates Euclidian norm of only real-parts of vec
    Double norm2Real(Vector<Complex> &vec);

    //! Performs a forward simulation with exact data, adds to results alternating +- 10 Percent
    void calcSyntheticData(Vector<Complex> & y_hat);

    void createMaterialTensorMatrices(Vector<Double> & parameter, Matrix<Double> & couplingMatrix, Matrix<Double> & dielectricMatrix, Integer spaceDim);

    //! Tests, if JacobiMatrix is more or less approximated by F(p)-F(p+delta)/delta
    void testJacobiMatrix(Vector<Complex> & F_hat, Matrix<Complex> & JacobiMatrix, Vector<Double> & parameter,BCs * ptBCs,MaterialData * ptMaterial, Vector<Double> & parameterIncrement, Vector<Complex>& solElecPot,Vector<Complex> &solMechDispl);


    Boolean considerMechDeformation;
    Matrix<Double> * piezoMatrix;
    Integer dofs;
    Integer numNodes;
    BaseSystem * ptAlgsys;
    BCs * ptBCs;
    Grid * ptGrid;
    NodeEQN * ptNodeEqn;
    Assemble * ptAssemble;
    StdVector<std::string> subdoms;
    BaseNodeStoreSol * sol;
    Domain * ptDomain;

    Vector<Complex> solElecPot;
    Vector<Complex> solMechDispl;
    Vector<Complex> algSysSolVector;
    MaterialData * ptMaterial;

    Vector<Double> parameter;
    Vector<Double> parameter_new;
    Vector<Double> parameterIncrement;
    Vector<Double> omegas;
    Vector<Double> freqs;
    Vector<Double> real, imag;
    Vector<Complex> amplitude_phase;
    Vector<Complex> F_hat;
    Vector<Complex> y_hat;
    Vector<Complex> s_0;    
    Vector<Double> bas;
    Vector<Complex> res_NE, res_NE_new, lin_res;
    Vector<Complex> res;
    Vector<Complex> bas_bar;
    Vector<Complex> s;
    Vector<Double> scaling;
    Double alpha_m, beta_m;
    Matrix<Complex> JacobiMatrix;
    Matrix<Complex> adjJacobiMatrix;
    Double voltage, thickness, radius, delta, meanValueMechDeformation, anorm, tau;
    Integer nrMeasuredData;
    Integer nrParameter;
    Integer nrSol, sizeSol;
    Double norm_res;
    Double eta;
    Double overall_res;
    Vector<Complex> overall_res0;

    Matrix<Complex> completeSolOf_F;
    Matrix<Complex> allElemsVec;
 

  private:



  }; // end of class piezoParamIdent

}



#endif
