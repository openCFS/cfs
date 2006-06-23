#ifndef FILE_MAGNETICPDE
#define FILE_MAGNETICPDE

#include "SinglePDE.hh" 
#include "Forms/magforceop.hh"

namespace CoupledField
{

  // Forward declaration of classes
  class Coil;

  //! Class for defining the magnetic field in 2D
  class MagPDE : public SinglePDE
  {
  public:

    MagPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut);

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagPDE();


    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! read special boundary conditions (coils, magnets)
    virtual void ReadSpecialBCs();

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators();

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();


    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! do PostProcessing step
    virtual void PostProcess();

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(const UInt kstep = 0,
                                    const Double asteptime = 0.0,
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0);
    
    //! write history results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteHistoryInFile(const UInt kstep,
                            const Double asteptime,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);

    //! computes the electric energy for each subdomain
    void CalcEnergy();


    //! calculates nodal forces
    void CalcNodeForce(ElemStoreSol<Double> & force, 
                       StdVector<UInt> & nodes, 
                       StdVector<Elem*> & elems,
                       StdVector<StdVector<ShortInt> > & isBoundaryNode,
                       StdVector<StdVector<UInt> > &elemNodeToCouplingNode)
    {Error("CalcNodeForce not implemented",__FILE__,__LINE__);}



    // ======================================================
    // COUPLING SECTION
    // ======================================================


    //! initalize PDE coupling
    virtual void InitCoupling(PDECoupling * Coupling);
  
    //! calculate coupling terms
    virtual void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output);

    //! computation of Lorentz force
    void CalcNodeForceLorentz(Vector<Double> & force, 
                              StdVector<StdVector<UInt> > & 
                              elemNodeToCouplingNode,
                              UInt actCoupling, 
                              UInt numCouplingNodes);


  protected:

    //! Query parameter object for information on coils
    void ReadCoils();
    
    //! Query parameter object for information on permanent magnets
    void ReadMagnets();

    //! Init the time stepping
    virtual void InitTimeStepping();

    //!
    void ComputeUI(Vector<Double>& uiSD);
  
    void WriteUI2File(Vector<Double>& uiSD);

    //! contains first derivative of magnetic vector potential
    NodeStoreSol<Double> solDeriv1_;

    //!
    Vector<Double> RhsLinVal_;

    ElemStoreSol<Double> B_;  //!< conatins magnetic field
    ElemStoreSol<Double> Jeddy_;  //!< conatins eddy currents field
  
    // ---- Electric Force variables ---
    //!  stores Magnetic force (VWP) of each element
    NodeStoreSol<Double> Force_;        

    //! vector of vectors conaining Elements with acting force
    StdVector<StdVector<Elem*> > F_Interface_; //!<

    //! vector containing flag array for element boundary nodes
    StdVector<StdVector<StdVector<ShortInt> > > isBoundaryNode_; 

    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; 

    //! contains number of surface nodes per element
    StdVector<StdVector<UInt> > numBoundaryNodes_;               

    // =======================================================================
    //   COILS
    // =======================================================================
    
    //@{ \name Attributes related to coils
    
    //! Names of coils resp. their subdomains
    StdVector<RegionIdType> coilRegionId_;  
    
    //! Parameters of the individual coils;
    StdVector<Coil*> coilDef_;
    
    //@}

    // =======================================================================
    //   PERMANENT MAGNETS
    // =======================================================================
    
    //@{ \name Attributes related to permanent magnets
    
    //! Subdomains containing permanent magnets
    StdVector <RegionIdType> magnetsDomain_;
    
    //! x-component of direction of magnetisation for each magnet
    
    //! x-component of direction of magnetisation for each magnet
    //! \todo As suggested by Fred Hofer, the direction of magnetisation of a
    //! permanent magnet must now be specified in the XML parameter file and
    //! no longer in the material data file. While magneticPDE already reads
    //! these data, they are not yet used in the simulation.
    StdVector<Double> magnetsOriX_;
    
    //! y-component of direction of magnetisation for each magnet
    StdVector<Double> magnetsOriY_;

    //! z-component of direction of magnetisation for each magnet
    StdVector<Double> magnetsOriZ_;

    //@}

    //postprocessing
    //! contains the subdomains, on which the magnetic field is computed
    StdVector<RegionIdType> calcBfield_;  
    
    //! contains the subdomains, on which the magnetic energy is computed
    StdVector<RegionIdType> calcEnergy_;  

    //! contains the subdomains, on which the eddy currents are computed
    StdVector<RegionIdType> calcEddy_;  

    //! contains the name of the nodes, where the force is computed via VWP
    StdVector<std::string> calcForceVWP_;  

    //! nodes, for wich the force has to be computed
    StdVector<UInt> ForceNodes_;           

    //! force operator (for coupling as well as postprocessing)
    MagForceOp* ForceOpVWP_;

    //! file for informational output of coils
    std::ofstream * UIfile_; 

    //! name of file for saving current/voltage values
    std::string UIfilename_;
   
  private:

    //! List of regions with non-linearity
    StdVector<std::string> nonLinType_;

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the electrostatics PDE 
    //! are given in the following table. Here 'Keyword' and 'Result Type' 
    //! refer to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the MagPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Physical quantity</b></td>
    //!     <td><b>Keyword in parameter file</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>Magnetic flux density \f$\vec{B}\f$</td>
    //!     <td>bfield</td>
    //!     <td>Element result</td>
    //!     <td>calcBfield_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>Magnetic energy \f$W\f$</td>
    //!     <td>energy</td>
    //!     <td>Element result</td>
    //!     <td>calcEnergy_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>Eddy current density \f$I_E\f$</td>
    //!     <td>eddy</td>
    //!     <td>Element result</td>
    //!     <td>calcEddy_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>depends on formulation</td>
    //!     <td>none -- always true</td>
    //!     <td>Nodal results</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td align="center">---</td>
    //!     <td>none -- always false</td>
    //!     <td>Nodal results</td>
    //!     <td>deriv1_</td>
    //!   </tr>
    //!   <tr>
    //!     <td align="center">---</td>
    //!     <td>none -- always false</td>
    //!     <td>Nodal results</td>
    //!     <td>deriv2_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagPDE
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace
#endif

