#ifndef FILE_MECHANIC2DPDE_2001
#define FILE_MECHANIC2DPDE_2001

#include "basepde.hh"

 
namespace CoupledField
{

  //! Class for mechanic equation 2D
  /*! 
    This class is derived from class BasePDE. It is used for solving acoustic equation on one time step.  We set rules for assembling global system matrix according to weak form of PDE, define right hand side and set boundary conditions. Then we cause one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

  class Mech2dPDE: public BasePDE
  {
  public:
    
    //!  Constructor. here we read integration parameters
    /*!
      \param  aptalgsys pointer to class Algebraic system
      \param aGrid pointer to grid
      \param aMatFile pointer to class Material. material data.
      \param aInFile pointer to class FileType. input data.
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    Mech2dPDE(AbstractAlgebraicSys * aptalgsys, Grid * aGrid , Material * aMatFile, TimeFunc * aTimeFunc ,FileType * aInFile, WriteResults * aOutFile );

    //!  Deconstructor
    virtual ~Mech2dPDE(){;}


    //! specify type of system matrix for AlgebraicSystem
    /*!
      \param matrixtype out: 0..NOCLASS, 1..RSPARSE, 2..CSPARSE, 3..RBLOCK, 4.. CBLOCK = 0,
      5..RFULL, 6..CFULL, 7..MIXED
      \param matrixsystype out:define need we memory for different types of element-matrix or not                     
      \param graphtype out: type of graph
      \param numdofpernode out: number of dof per node
      \param numdirichlets out:number of nodes for dirichlets conditions
      \param numconstraints out:number of nodes for constraints conditions
    */
    virtual void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

    //! set information for algebraic system about PDE. set matrix factors
    virtual void SetMatrixFactors();

    //! specify type of system matrix for AlgebraicSystem
    /*!
      \param level (input) level of Grid
    */
    virtual void SetupMatrices(const Integer level, BCs * ptBCs=NULL);

    //! set boundary condition
    /*!
      \param ptBCs pointer to boundary condition
      \param level level of grid
      \param update indicator: do we update boundary condition in algebraic system ot set new
      \param atime time of calculation
    */
    virtual void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

    /*!
      \param ptBCs pointer to class with data about boundary condition
      \param level level of grid
    */
    virtual void SolveStepStatic(BCs * ptBCs ,const Integer level);

    //! solve one step for transient problem 
    /*!
      \param ptBCs pointer to class with data about boundary condition
      \param kstep number of calculating step
      \param steptime time of calculation
      \param level level of grid
      \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
    */
    virtual void SolveStepTrans(BCs * ptBCs ,const Integer kstep, const Double steptime, const Integer level, const Boolean updatesysmat);

    //! write results in file
    virtual void WriteResultsInFile();


    //! return pointer to vector with solution
    virtual const Vector<Double> & getS() const { return sol_;};

  protected:

    //!
    Integer dofspernode_;

    //!
    Grid * ptgrid_;

    //! Calculation parameters for Newmark method
    virtual void CalcParameters(const Double dt);

    //!  calculation of coefficient.
    /*!
      \param coeffmass coefficient before mass matrix
      \param coeffstiff coefficient before stiffness matrix
    */
    void CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff);

    //! coefficients from Newmark method
    Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;

    //! Integration parameters
    Double alpha_,gamma_, beta_;

    //! store solution, 1st derivative , 2nd derivative solution
    Vector<Double> sol_, sol_der1_, sol_der2_, sol_old_, sol_der1_old_, sol_der2_old_;  

    //! Last time on which we have calculated solution
    Double lasttimecalc_;

    //! Number of last timestep on which we have calculated our solution
    Integer laststepcalc_;

    //! size of solution and etc.
    Integer size_;

    //! function for RHS
    Integer arg_rhs_;
    pfn1var ptRHSFnc_;
    std::vector<Double> directionFnc_;

    //! Indicator: is there RHS function
    Boolean SetRHSFnc;

  };

} // end of namespace
#endif
