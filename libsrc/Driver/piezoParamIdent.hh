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

    std::fstream allMeasuredData;

    void SolveProblem();

protected:
  //! Calculates the parameter to soution map F(p^k) at Newton iteration step k
    void createF(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_hat);

  //! Calculates an approximation of the Jacobi Matrix of parameter to solution operator F 
    void createJacobiMatrix(MaterialData * ptMaterial, BCs * ptBCs, Vector<Complex> & F_ha, Vector<Double> & parameterIncrement, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl);

  //! Calculates explicitely the Adjoint operator of F'
    void createAdjointJacobiMatrix(Vector<Double> & parameterIncrement,Vector<Double> &  parameter, Matrix<Complex> & JacobiMatrix, Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl, Vector<Double> & freqs, Matrix<Complex> & adjJacobiMatrix);

  //! Method which reads Data from file measuredData.dat. The file contains measurements of amplitude, frequency, further according information concerning the piezhoelectric body (radius, thickness, ...)
     void readMeasuredData(Vector<Double> & freqs, Vector<Double> & real, Vector<Double> & imag ,Vector<Double> & parameter, Double & voltage, Integer & nrMeasuredData, Double & thickness, Double & radius, Double & delta);

    void updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial);

    void updateRHS(Vector<Complex> & solElecPot, Vector<Complex> & mechDisplacement, Double omega);

    void updateRHS(Vector<Complex> & RHSsol);

    void typeOutSolutionOnScreen(Vector<Complex> & solElecPot,Vector<Complex> & solMechDispl);

    void calcInitialResidual(Vector<Complex> & res, Vector<Complex> & y_hat, Vector<Complex> & PHI_p, Integer fstep, Vector<Complex> & solElecPot, Double & meanValueMechDeformation);

    void measureMechDeformationInZ_Direction(Vector<Complex> & mechDisplacement, Double & Radius, Double meanValueMechDeformation, int dof);

    void calcNorm2Resid(Vector<Complex> &res, Double & anorm, Integer nrMeasuredData);

		void createMaterialTensorMatrices(Vector<Double> & parameter, Matrix<Double> & couplingMatrix, Matrix<Double> & dielectricMatrix, Integer spaceDim);

    Matrix<Double> * piezoMatrix;
    Integer dofs;
    Integer numNodes;
    BaseSystem * ptAlgsys;
    BCs * ptBCs;
    Vector<Complex> solElecPot;
    Vector<Complex> solMechDispl;
    Vector<Complex> algSysSolVector;
    MaterialData * ptMaterial;

    Vector<Double> parameter;
    Vector<Double> parameterIncrement;
    Vector<Double> omegas;
    Vector<Double> freqs;
    Vector<Double> real, imag;
    Vector<Complex> amplitude_phase;
    Vector<Complex> F_hat;
    Matrix<Complex> JacobiMatrix;
    Matrix<Complex> adjJacobiMatrix;
    Double voltage, thickness, radius, delta, meanValueMechDeformation, anorm, tau;
    Integer nrMeasuredData;
    Integer nrParameter;
    Integer nrSol, sizeSol;


    //void getParamsFromPiezoParamIdent(Vector<Double> &parameter);
    //	void updateParams(Vector<params>, MaterialData * material)

    // void newtonCG(Vector<Double> &params, Double &omegas);

  private:

    //BasePDE * actPDE;
    //MaterialData * ptMaterial;

  }; // ned of class piezoParamIdent

}



#endif
