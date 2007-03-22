// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MAGEDGEPDE_2002
#define FILE_MAGEDGEPDE_2002

#include "SinglePDE.hh" 

namespace CoupledField
{

  //! Class for electrostatic equation in 3D (no adaptivity)
  class MagEdgePDE: public SinglePDE
  {
  public:

    //! Constructor. here we read integration parameters
    /*!
      \param 
      \param aGrid pointer to grid
    */
    MagEdgePDE(Grid * aptgrid );

    //! Deconstructor
    virtual ~MagEdgePDE();

    //! define discrete PDE
    virtual void DiscreteParamsPDE();

    //! set information for algebraic system about PDE. set matrix factors.
    void SetMatrixFactors();


    //! set up the edge graph
    void SetupEdgeGraph();

    //! define algebraic system 
    /*!
      \param AS_sysid id of PDE in algebraic system
    */
    virtual void SetAlgSys(const Integer AS_sysid);

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();

    //! compute the number of edge) dirichlets
    void EvalNumDirichlet();

    //! get the edge numbers and theis signs (from AlgSys)
    /*!
      \param nodes node numbers of the element
      \param epos integer array conatining the edge numbers
      \param esign integer array containing the edge signs
    */
    void GetEdgeNumber(Integer *nodes, std::vector<Integer>& epos, 
                       std::vector<Integer>& esign, BaseFE * ptElem);

    //! setup element matrices for AlgebraicSystem for assembling procedure
    void SetupMatrices( );
  
    //! set boundary condition
    /*!
      \param level level of grid
      \param update indicator: do we update boundary condition in algebraic system ot set new
      \param atimestep time step of claculation
    */
    virtual void SetBCs(const Integer level, const Integer update, const Double atimestep);

    //! solve one step for static problem 
    /*!
      \param level level of grid
    */
    virtual void SolveStepStatic(const Integer level, const Double aTime=0);

    virtual void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, 
                                const bool updatesysmat)
    { 
      Error("Makes no sense for Electrostatics to perform transient step",__FILE__,__LINE__);
    }


    //! write results in file
    virtual void WriteResultsInFile();

    //! return pointer to vector with solution
    virtual const Array<Double>& getS() const { ;}

    /// calc element quantities
    void PostProcess(const Integer level);

    /// solve one harmonic step
    virtual void SolveStepHarmonic(const Integer level);
  
    /// parameters necessary to describe coils
    struct coilDefStruct
    {
      Integer iDir;
      Double  current;
      Double  coilArea;
      Double currentPhase;
      std::vector<Double> coilMidPt;
    };
  

  protected:
    /// setup source term
    void SetupRHS(const Integer level);

    /// reads all data in the config-file belonging to coils
    void MagEdgePDE::ReadCoils();
  
    /// correct the sign in the element matrix due to orientation of edges
    /*!
      \param elemmat Element matrix
      \param esign   Vector of edge orientations
    */
    void CorrectEdgeDir(Matrix<Double> & elemmat, std::vector<Integer>& esign);
  

    //   Array<Double> E_;   //!< store Electric Field of each element
    //   Array<Double> Force_;  //!< stores Electric pressure force of each element
    //   Integer size_;        //!< size of solution (number of equations)
    //   Integer nElements_;
    UInt dofsperedge_;  //!< number of unknowns per edge
    Vector<Double> solRe_;      //!< store real part of solution,
    Vector<Double> solIm_;      //!< store imaginary part of solution
    Integer size_;              //!< size of solution (number of edges)
    Integer NumElems_;          //!< number of elements belonging to PDE
    Integer * EdgeDir_;         //!< stores the Dirichlet-edges
    Integer numEdgedir_;        //!< number of Dirichlet edges
    Vector<Double>* bFieldRe_;  //!< real part of vector of magnetic field induction
    Vector<Double>* bFieldIm_;  //!< imaginary part of vector of magnetic field induction
    Vector<Double>* magVecPotRe_; //!< real part of vector of magnetic magnetic vector potential
    Vector<Double>* magVecPotIm_; //!< imaginary part of vector of magnetic magnetic vector potential
    std::vector <std::string> coilDomain_;  //!< name of all subdomains containing coils
    // Double freq_;            //!< excitation frequency for harmonic analysis

  

    std::vector<struct coilDefStruct> coilDef_; //!< vector of paramters describing coils
  
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgePDE
  //! 
  //! \purpose 
  //! This class is derived from class SinglePDE. It is used for solving 
  //! the magnetic field equation in 3D
  //!
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status Not used
  //! 
  //! \unused 
  //! 
  //! \improve
  //! The concept of edge elements has to be re-thougt. We need a 
  //! generalization of the Unknown/dof concept, because at the moment
  //! each integrator only gets a vector with edge nodes from grid and 
  //! calculates the according matrix. \n
  //! Therefore the following is needed:
  //! -# Introduce a strict separtaion ofreference element (BaseFE) and the
  //! according degrees of freedom (nodes, edges, faces, ...)
  //! -# Introduce an additional equation class, derived from BaseEQN, which
  //! e.g. is repsonsible for giving equation numbers to edges, faces etc.
  //! 

#endif

} // end of namespace
#endif
