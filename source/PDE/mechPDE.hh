// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWBASEMECHPDE
#define FILE_NEWBASEMECHPDE

#include <map>

#include "SinglePDE.hh"
 
namespace CoupledField
{

  // forward class declarations
  class MechVolForceInt;
  class BaseForm;

  //! Class for mechanic equation (no adaptivity)
  class MechPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    MechPDE( Grid *aGrid, ParamNode* paramNode );

    //!  Deconstructor
    virtual ~MechPDE();

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! Initialize NonLinearities
    void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! Read special results definition
    void ReadSpecialResults();

    /// returns a stiffness integrator appropriate to the actual problem (e.g.3D)
    BaseForm * GetStiffIntegrator(BaseMaterial* actSDMat,
                                  RegionIdType regionId,
                                  bool reducedInt=false);
  
    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! calculate coupling terms
    void CalcOutputCoupling();
  
    //! returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);

    /// setup source term
    void SetupRHS( );
  

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> results );

  protected:

    // ======================================================
    // POSTPROCESSING METHODS
    // ======================================================
  
    //computes mechanical deformation energy
    template <class TYPE>
    void CalcEnergy( shared_ptr<BaseResult> vals );

    //computes mechanical stresses
    template <class TYPE>
    void CalcStresses(  shared_ptr<BaseResult> vals );

    //computes mechanical strains
    template <class TYPE>
    void CalcStrains(  shared_ptr<BaseResult> vals );

    //! compute volume above a deformed surface
    template <class TYPE>
    void ComputeVolDefSurf( shared_ptr<BaseResult> vals );
    
    //! computes averaged volume of an deformed element
    template<class TYPE>
    TYPE ComputeVolElem(BaseFE * ptSurfEl, Matrix<Double>& SurfCoord, 
                        Vector<TYPE> disp);
    

    //! Nodestoresol for RHS
    BaseNodeStoreSol * rhs_;
    
    // ======================================================
    // SCRIPTING SECTION
    // ======================================================
    //@{
    //! \name Scripting Methods
    
    //! Register scriptable functions
    void RegisterFunctions();
    //@}
    
    //! Obtain information on desired output quantities from parameter file
  
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the mechanics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the MechPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Keyword</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>displacement</td>
    //!     <td>nodeResults</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>velocity</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>acceleration</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv2_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();
  
    //! Init the time stepping
    void InitTimeStepping();

  private:

    // calc rhs coupling to acoustic pde
    // void CalcAcousticCouplingRHS(StdVector<Elem*> * couplingElems,
    // Vector<Double>& forceOnElem);
  
    /// calc rhs coupling to acoustic pde
    void CalcAcousticCouplingRHS( StdVector<Elem*> * couplingElems, 
                                  StdVector<BaseMaterial*> & materials, 
                                  StdVector<UInt>& couplingNodes,
                                  Vector<Double> & forceOnElem,
                                  UInt couplingdof );
                               
  

    /// does a line search and returns the optimal residual norm
    Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                      Double& etaLineSearch, bool trans=false);


    /// Write nonlin iteration norms to the cla-file
    void WriteClaNlNorms( const UInt iterationCounter,
                          const Double residualL2Norm,
                          const Double extForcesL2Norm, const Double residualErr, 
                          const Double solIncrL2Norm, const Double actSolL2Norm, 
                          const Double incrementalErr );
  

    //! read in the domains with prestressing
    void ReadPreStressing();

    //! read in volume sources
    void ReadRegionLoads();

    //! read in surface stress boundary conditions
    void ReadSurfStress();

    //! read pressure loads
    void ReadPressureLoads();

    //! read in softening types
    void ReadSoftening();

    //! define available result types
    void DefineAvailResults();

    //! Class defining data needed for region loads
    class RegionLoad {

    public:

      //! Constructor
      RegionLoad( UInt dim, bool isaxi );

      //! Print region definition to info-file
      void Print( bool onlyHeader, std::string pdeName );
      
      //! Returns the RHS-integrator
      MechVolForceInt *  GetIntegrator();
      
      // ----------------------------
      //   Data members
      // ----------------------------

      //@{
      // \name Data members
      
      //! Name of region
      std::string name;

      //! Value of load
      StdVector<std::string>  value;

      //! Phase value
      std::string phase;

      //! Name of reference coordinate system
      std::string refCoord;

      //! Type of load (total/unit)
      std::string type;

      //! Volume of region
      Double volume;

      //! Flag for axisymmetry
      bool isAxi_;
      //@}

    };

    struct SurfStress {
      
      //! Name of surface elements
      RegionIdType surface;

      //! Name of neighbouring region
      RegionIdType region;
      
      //! Vector of load
      Vector<Double> stress;

      //! Phase value
      std::string phase;
    };
      
    
    //! List of region loads
    std::map<RegionIdType, RegionLoad> regionLoads_;

    //! List of surface stresses
    std::map<RegionIdType, SurfStress> surfStresses_;

    //! Read information about springs and define related integrators
    void DefineSprings();

    /// stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    /// returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);
  
    /// flag for reduced Integration for each subdomain
    StdVector<std::string> reducedIntegration_;

    /// returns the solution matrix belonging to all nodes of the actual element
    void GetSolOfElement( Matrix<Double>& elDisp, StdVector<UInt>& connect_PDE);

    //! Number of dimension for stresses
    UInt stressDim_;
    
    //! Flag indicating the use of fractional damping
    //bool fracDamping_;

    //! surface of pressure loads
    StdVector<RegionIdType> pressSurf_;  

    //! values of the pressure loads
    StdVector<std::string>  pressVals_; 

    //! phase of the pressure loads  
    StdVector<std::string>  pressPhase_; 

    //! list of prestressing types
    std::map<RegionIdType,std::string> preStressList_;

    //! prestress-values: numSubdoms x 3 
    std::map< RegionIdType, Vector<Double> > preStressVal_; 

     //@{ \name Attributes related to post-processing

     //! Contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;
  
    //! Contains mechanic acceleration
    NodeStoreSol<Double> solDeriv2_;

    //! Contains the directions for which the deformed volume is computed
    std::map<shared_ptr<EntityList>,std::string> volAboveDefSurfDir_;
    
    //! Stores softening for each region
    std::map<RegionIdType, std::string> regionSoftening_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MechPDE
  //! 
  //! \purpose 
  //! This class defines the mechanical field PDE and the according 
  //! postprocessing methods.
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

