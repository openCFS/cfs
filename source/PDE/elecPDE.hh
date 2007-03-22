// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ELECPDE_NEW
#define FILE_ELECPDE_NEW

#include "SinglePDE.hh" 
#include "Forms/elecforceop.hh"

namespace CoupledField
{

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  
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

    //! Constructor
    /*!
      \param 
      \param aGrid pointer to grid
      \param aGrid pointer to class Grid
    */
    ElecPDE( Grid* aptgrid, ParamNode* paramNode );

    //! Destructor
    virtual ~ElecPDE(){};

    //! Initialize NonLinearities
    void InitNonLin();

    //! Define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping() {;};

    //! Nothing to do
    void SetTimeStep(const Double dt) {;};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> result );

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

    //! Define availabe result types
    void DefineAvailResults();

    //! Calculate electrid field intensity
    template <class TYPE>
    void CalcElectricField( shared_ptr<BaseResult> vals );

    //! Calculates the polarization vector
    void CalcPolarizationField( shared_ptr<BaseResult> vals );

    //! Calculate electric charges
    template <class TYPE>
    void CalcCharges( shared_ptr<BaseResult> vals );

    //! Computes the electric energy for each subdomain
    template <class TYPE>
    void CalcEnergy( shared_ptr<BaseResult> vals );

    //! Contains the (Volume) subdomains next to the surface
    //! elements where the charges are computed
    std::map<shared_ptr<EntityList>,RegionIdType> chargeNeighborRegion_;

    //! Electric impedances
    StdVector<Impedance> impedances_;

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
    void InitStoreResults();

    //! flag for piezo-coupling
    bool isPiezoCoupled_;

    //! force operator (for coupling as well as postprocessing)
    ElecForceOp* ForceOp_;

    //! vector containing regionIds of non-conforming interfaces
    StdVector<RegionIdType> ncIFaces_;
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
