#ifndef FILE_ELECPDE_NEW
#define FILE_ELECPDE_NEW

#include "SinglePDE.hh" 
#include "Forms/elecforceop.hh"

namespace CoupledField
{

  //! Class for electrostatic equation (no adaptivity)
  class ElecPDE : public SinglePDE {

  public:

    //! Helper struct for defining a general impedance
    struct Impedance {

      //! Constructor
      Impedance() : 
        node1(0), node2(0), resistance(0.0),
        inductance(0.0), capacitance(0.0) {};
      
      //@{
      //! Node numbers the impedance is connected with
      UInt node1, node2;
      //@}

      //@{
      //! Defining quantities of the impedance
      Double resistance, inductance, capacitance;
      //@}
    };

    //! Constructor. here we read integration parameters
    /*!
      \param 
      \param aGrid pointer to grid
      \param aGrid pointer to class Grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    ElecPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut);

    //! Destructor
    virtual ~ElecPDE(){};

    //! Define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( );

    //! Define the SolveStep-Driver
    virtual void DefineSolveStep();

    //! Init the time stepping: nothing to do
    virtual void InitTimeStepping() {;};

    //! Nothing to do
    virtual void SetTimeStep(const Double dt) {;};

    //! Read special boundary conditions
    void ReadSpecialBCs();


    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! Do PostProcessing step
    virtual void PostProcess( );

    //! Write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(const UInt kstep = 0,
                                    const Double asteptime = 0.0,
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0);
    
    //! Write history results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteHistoryInFile(const UInt kstep,
                            const Double asteptime,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);

   

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! Initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! Calculate coupling terms
    void CalcOutputCoupling();

    //! Returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);
  
    //! Turn the piezo coupling on

    //! Triggers the correct assembly of the electrostatic block in a 
    //! piezo-coupled simulation, because the coupled electrostatic block
    //! is negative compared to the normal one
    void SetPiezoCoupling();


  protected:

    //! SubType of electrostatic section
    std::string subType_;
    
    // *****************
    //  POSTPROCESSING
    // *****************

    //! Calculate electrid field intensity
    template <class TYPE>
    void CalcElectricField();

    //! calculates the polarization vector
    void CalcPolarizationField();

    //! Calculate electric charges
    template <class TYPE>
    void CalcCharges();

    //! Computes the electric energy for each subdomain
    template <class TYPE>
    void CalcEnergy();

    //! Contains the subdomains, on which the electric field is computed
    StdVector<RegionIdType> calcEfield_; 

    //! Contains the subdomains, on which the electric field is computed
    StdVector<RegionIdType> calcPolarization_; 
    
    //! Contains the subdomains, on which the electric energy is computed
    StdVector<RegionIdType> calcEnergy_;  

    //! Contains the subdomains, on which the electric charges  are computed
    StdVector<RegionIdType> calcCharges_;

    //! Contains the (Volume) subdomains next to the surface
    //! elements where the charges are computed
    StdVector<RegionIdType> chargeNeighborRegion_;

    //! Conatins electric field
    BaseElemStoreSol * E_;  

    //! names for elements, for which polarization is saved
    StdVector<std::string> saveElemElecIntensityHist_;

    //! Conatins electric field
    BaseElemStoreSol * P_;  

    //! names for elements, for which polarization is saved
    StdVector<std::string> saveElemPolarizationHist_;

    //! Contains electric charges
    BaseElemStoreSol * charges_;

    //! Electric impedances
    StdVector<Impedance> impedances_;

  private:

    //! List of regions with non-linearity
    StdVector<std::string> nonLinType_;

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
    bool isPiezoCoupled_;

    //! force operator (for coupling as well as postprocessing)
    ElecForceOp* ForceOp_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class ElecPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! electrostatic equation in 3D. 
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
