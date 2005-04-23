#ifndef FILE_ELECPDE_NEW
#define FILE_ELECPDE_NEW

#include "scalarnodeEQN.hh"
#include "SinglePDE.hh" 


namespace CoupledField
{

  //! Class for electrostatic equation in 3D (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving
    electrostatic equation in 3D. 
  */

  class ElecPDE : public SinglePDE {

  public:

    //! Constructor. here we read integration parameters
    /*!
      \param 
      \param aGrid pointer to grid
      \param aBCs pointer to Boundary condition object
      \param aGrid pointer to class Grid
      \param aInFile pointer to class FileType. input data.
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    ElecPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
	    FileType *aptFileType, WriteResults *aptOut);

    //! Deconstructor
    virtual ~ElecPDE(){};

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators(const Integer level);

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //! Init the time stepping: nothing to do
    virtual void InitTimeStepping() {;};

    //! nothing to do
    virtual void SetTimeStep(const Double dt) {;};

    //! return size of solution
    virtual Integer getSize() const 
    { return numPDENodes_*dofspernode_;}


    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! do PostProcessing step
    virtual void PostProcess(const Integer level);

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(const Integer kstep = 0,
				    const Double asteptime = 0.0,
				    Integer stepOffset = 0,
				    Double timeOffset = 0.0);
    
    //! computes the electric energy for each subdomain
    void CalcEnergy();

    //! callculates nodal forces
    void CalcNodeForce(Vector<Double> & force, 
		       StdVector<Integer> & nodes, 
		       StdVector<Elem*> & elems,
		       StdVector<StdVector<ShortInt> > &isBoundaryNode,
		       StdVector<StdVector<Integer> > &elemNodeToCouplingNode);

    //!
    void CalcInterfaceForces(Integer actCoupling);

    //! GET SOLUTION AT ALL NODES OF AN ELEMENT
    //void GetSolOfElement( Vector<Double>& elpot, StdVector<Integer>& connect_PDE);


    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! calculate coupling terms
    void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    Boolean HasOutput(SolutionType output);
  
    //! turn the piezo coupling on

    //! triggers the correct assembly of the electrostatic block in a 
    //! piezo-coupled simulation, because the coupled electrostatic block
    //! is negative compared to the normal one
    void SetPiezoCoupling();
    

  protected:

    /// calculated the electric field at the integration points of the couple element
    void CalcEfieldAtCoupleElemIP(Elem * actVolElem,
				  Elem * actCoupleElem,
				  Vector<Double>& coordAtIP, 
				  StdVector<Integer>& boundNodesOfVolElem,
				  Vector<Double>& tempE);
  


    //  Boolean nonLinGeo_;  //! switch for geometric update 
  
    // ---- Electric Force variables ---
    ElemStoreSol<Double> Force_;        //!< stores Electric force of each element
    StdVector<StdVector<Elem*> > F_Interface_; //!<vector of vectors conaining Elements with acting force
    StdVector<StdVector<StdVector<ShortInt> > > isBoundaryNode_; //!< vector containing flag array for element boundary nodes
    StdVector<StdVector<StdVector<Integer> > > elemNodeToCouplingNode_; //!< assigns each coupling element node the according Coupling Node number
    //  StdVector<StdVector<Integer> > numBoundaryNodes_;               //!< contains number of surface nodes per element
    
    // *****************
    //  POSTPROCESSING
    // *****************

    //! callculate electrid field intensity
    void CalcElectricField();

    //! calculate electric charges
    void CalcCharges();

    //! contains the subdomains, on which the electric field is computed
    StdVector<std::string> calcEfield_; 
    
    //! contains the subdomains, on which the electric energy is computed
    StdVector<std::string> calcEnergy_;  

    //! contains the subdomains, on which the electric charges  are computed
    StdVector<std::string> calcCharges_;

    //! contains the (Volume) subdomains next to the surface
    //! elements where the charges are computed
    StdVector<std::string> chargeNeighborRegion_;

    //! conatins electric field
    ElemStoreSol<Double>  E_;  

    //! contains electric charges
    ElemStoreSol<Double>  charges_;

    // for check: own solver
    Boolean SolverCFS_; //<! parameter indicator: TRUE, if you want to use Solver CFS. reading from config-file
    Matrix<Double> sysmat_;
    Vector<Double> vecrhs_;

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the electrostatics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the ElecPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Keyword</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>potential</td>
    //!     <td>nodeResults</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>efield</td>
    //!     <td>elemResults</td>
    //!     <td>calcEfield_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>energy</td>
    //!     <td>elemResults</td>
    //!     <td>calcEnergy_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();

    //! flag for piezo-coupling
    Boolean isPiezoCoupled_;

  };

} // end of namespace

#endif
