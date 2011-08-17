// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SINGLEPDE
#define FILE_SINGLEPDE

#include "PDE/StdPDE.hh"

#include <list>
#include <map>

#include "DataInOut/Scripting/scriptable.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/resultInfo.hh"
#include "Utils/elemstoresol.hh"
#include "Domain/bcs.hh"
#include "Utils/result.hh"

namespace CoupledField
{
  // forward class declaration
	class VolForceInt;
	class SpaceErrorEstimator;
  class BasePairCoupling;
  class DirectCoupledPDE;
  class Assemble;
  class BaseForm;
  class PDEMemento;

  
  //! Base class for all kinds of single field problems.

  class SinglePDE : public StdPDE, public Scriptable
  {
  
  public:

    // friend declaration
    friend class BasePairCoupling;
    friend class DirectCoupledPDE;

    
    //! typedefs for result handling
    typedef std::set<shared_ptr<ResultInfo> > ResultSet;
    typedef StdVector<shared_ptr<BaseResult> > ResultList;
    typedef std::map<shared_ptr<ResultInfo> , ResultList > ResultMap;

    //class RegionLoad;

    bool boolComplexMaterialData_;
    
    /** Initialize PDEs 
     * @param base a coupled PDE case we do not choose our own base */ 
    virtual void Init( UInt sequenceStep, InfoNode* base = NULL ); 

    // ---------------------- ***** --------------------------------

    //! destructor
    virtual ~SinglePDE();
  
    //! MpCCI gets the geometry
    virtual void PreparePDE4Computation() {;};
  
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! define algebraic system 
    virtual void DefineAlgSys();
 
  
    // ======================================================
    // ADATPTIVITY SECTION
    // ======================================================
#ifdef ADAPTGRID
    //@{
    //! \name Methods used for performing adaptivity
  
    //! test error of calculation. return true, if it is more then tolerance
    virtual bool TestError(const Integer level);
  
    //! refine mesh
    virtual void RefineMesh(const Integer level=0);

    // Does this method belong to postproc section?
    //! write information about relative error of calculation
    void WriteErrorInfo(WriteResults * ptmeshes);
  
    //@}
#endif

    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling (only done once)
    virtual void InitCoupling(PDECoupling * Coupling) = 0;

    //! reset coupling counters and data (done after each timestep)
    virtual void ResetCoupling();

    //! Fill in input coupling terms
    virtual void CalcInputCoupling();
  
  
    //! calculate coupling terms
    virtual void CalcOutputCoupling() = 0;

  
    // ======================================================
    // GET /SET  METHODS
    // ======================================================

    //! get algsys identification tag of PDE
    FeFctIdType GetPDEId()
    { return pdeId_; }

    //! return subtype
    virtual std::string GetSubType() {return subType_;}



    //! Set Direct coupling information
    virtual void SetDirectCoupling();

    //! return number of restraints
    UInt GetNumRestraints( );
 
    //! set boundary condition OBSOLETE
    void SetBCs();

    //! Transforms a given BoundaryCondition value according to Timestepping (i.e. TransientSim)
    virtual void TransformBC(Double& transVal, Double initValue, Integer eqnNumber);

    //! set special PDE dependent boundary conditions
    virtual void SetSpecialBCs(){ return; }

    //! Method for modifying an inhomogeneous boundary condition
    void SetIDBC( const std::string &name,
                  const std::string &dofString, 
                  const std::string &value, 
                  const std::string &phase );

    /** Helper method for ReadBCs() which reads the loads. It is also called
     * by the SIMP mechanism design optimization to read the outputs
     * @param loadNodes the potential empty array from the xml file
     * @param either the loads_ of StdPDE or for optimization */
    void ReadLoads(StdVector<ParamNode*> loadNodes, LoadList& out_list);
    
    //! write general defines (BCs, loads, etc.) to info-file
    void WriteGeneralPDEdefines();

    //! get the encapsulated state of the PDE
  
    //! returns the current state of the PDE (solution, derivative,
    //! coupling-objects) in an encapsulated object. This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the solution, the derivative and perhaps coupling 
    //! values like geometry update have to be passed. 
    //! The PDEMemento object encapsulates this information. 
    //! Later on the information can be given back to the PDE
    //! with the method SetMemento();
    //! \param memento (output) Object where the current state gets saved
    void GetMemento(shared_ptr<PDEMemento>& memento);
  
    //! set the encapsulated state of the PDE
  
    //! set the current state of this PDE (solution, derivative,
    //! coupling-objects) from an encapsulated object. This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the solution, the derivative and perhaps coupling 
    //! values like geometry update have to be passed. 
    //! The PDEMemento object encapsulates this information. 
    //! With this method the previous stored information can be set
    //! to the current PDE.
    //! \param memento (input) Previously saved state of the PDE
    //! \param useAsDirichlet (input) Usage type of values (start-value / 
             //!                      dirichlet value )
    void SetMemento( shared_ptr<PDEMemento>&  memento, 
                     bool useAsDirichlet );
                   

    //! write the PDE state (pdememento) to a restart file "simname_pdename.restart"
    void WriteRestart( );

    //! read the PDE state (pdememento)from a restart file: "simname_pdename.restart"
    void ReadRestart(UInt &startStep );

    //! Map containing the result types and the results
    ResultMap& GetResults() { return resultLists_; }
    
    /**<p>This is part of ReadStoreResults(). If candiate is defined in the xml file
     * it is added to resultLists_.</p>
     * <p>This method is to be called by ReadStoreResults() for every element in
     * availResults_. Additionally an Otimization instance calls when there a
     * result element defines one of the solution types optResult_1/2/3 in more detail  
     * @param candidate normally an element of the (mathematical) set availResults_
     * @return true if in xml and added */
    bool CheckStoreResult(shared_ptr<ResultInfo> canditate);

    //@{
    //! store the new solution returned by the algebraic system
    //! \param ptSol pointer to solution array
    //! \param size legnth of solution array
    void SaveSolution( const Double * ptSol, UInt size );
    void SaveSolution( const Complex * ptSol, UInt size );
    
    void SavePrevSolution( const Double * ptSol, UInt size );

    //@}

//    /** do the actual reading of loads, this is also called from optimization 
//     * @param bcNode paramnode that has "regionLoad" nodes as children 
//     * @param pressSurf StdVector containing the RegionLoads */
//    void ReadRegionLoadsFromXML(ParamNode* bcNode, std::map<RegionIdType, RegionLoad>& regionLoads_);

  protected:

  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
    */
    SinglePDE( Grid *aptgrid, ParamNode *paramNode );

    //! private copy constructor
    SinglePDE & operator= (const StdPDE & myPDE) {
      EXCEPTION( "Not implemented" );
      
      // For compiler
      return *this;
      ;}

   
    // ======================================================
    // INITIALIZATION METHODS
    // ======================================================

    //! define all computable results
    virtual void DefineAvailResults() {};  

    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();

    //! Read special store results information
    virtual void ReadSpecialResults(){};

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( )=0;

    //! Trigger calculation of results
    virtual void CalcResults( shared_ptr<BaseResult> result ) {};
    
    //! Calculate special results, not handled by resulthandler
    virtual void CalcSpecialResults( ) {};

    //! read damping information
    virtual void ReadDampingInformation( ){
    };
    
    //! read material data
    virtual void ReadMaterialData();

    //! read from config-file info about BCs
    void ReadBCs();

    // overloaded version of ReadBCs for special
    // boundary conditions in derived classes
    virtual void ReadSpecialBCs(){}
    
    //! read in volume sources
    void ReadRegionLoads();
    
    //! write results in file
    void WriteResultsInFile( const UInt kstep, 
                             const Double actTimeFreq );
    
    //! incorporate information of memento object, if set
    virtual void IncorporateMemento();

    //! Initialize NonLinearities
    virtual void InitNonLin(){};

    //! Define FeFunctions for this PDE according to the
    //! definition in the XML file
    virtual void DefineFeFunctions();

    //! Create FeSpaces according to formulation
    virtual std::map<SolutionType, shared_ptr<FeSpace> > CreateFeSpaces(std::string formulation) = 0;

    //@{
    //! Save load part of RHS to private variable
    void SaveRHS( const Double * ptSol, UInt size );
    void SaveRHS( const Complex * ptSol, UInt size );
    
  public:
//    
//    //! Class defining data needed for region loads
//    class RegionLoad {
//
//    public:
//
//      //! Constructor
//      RegionLoad( UInt dim, bool isaxi );
//
//      //! Print region definition to info-file
//      void Print( bool onlyHeader, std::string pdeName );
//      
//      //! Returns the RHS-integrator
//      VolForceInt *  GetIntegrator();
//      
//      // ----------------------------
//      //   Data members
//      // ----------------------------
//
//      //@{
//      // \name Data members
//      
//      //! Name of region
//      std::string name;
//
//      //! Value of load
//      StdVector<std::string>  value;
//
//      //! Phase value
//      std::string phase;
//
//      //! Name of reference coordinate system
//      std::string refCoord;
//
//      //! Type of load (total/unit)
//      std::string type;
//
//      //! Volume of region
//      Double volume;
//
//      //! Flag for axisymmetry
//      bool isAxi_;
//      //@}

//    };
    
  protected:
    
    //! List of region loads
     //std::map<RegionIdType, RegionLoad> regionLoads_;
     
    //@}

    // ======================================================
    // SCRIPTING SECTION
    // ======================================================
    //@{
    //! \name Scripting Methods
    
    //! Register scriptable functions
    virtual void RegisterFunctions();

    //! Wrapper for setting an inhomogeneous boundary condition
    void Wrap_IDBC();

    //! Wrapper for getting the nodal result
    void Wrap_GetValue();

    //! check if subdomain is a coupled piezo subdomain with hystersis
    bool IsRegionPiezoHyst( std::string regionName );
    
    //! check if PDE is a coupled piezo subdomain with hystersis
    bool BelongsPDE2PiezoHyst();    
    //@}
    
    // ======================================================
    // DATA SECTION
    // ======================================================

    // reads in the PML data
    void ReadDataPML(std::string& typePML, Matrix<Double>& inner, 
		     Double& dampPML, ParamNode * actNode);

    //! computes the PML layer dimensions
    void GetPMLLayerData(Matrix<Double>& inner, Matrix<Double>& outer,
			 RegionIdType regionId );
  
    // -----------------------------------------------------------------------
    // Storing information
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to storing information
    
    //! Copy the solution from a nodeStoresolution object to a solution object
    template<class TYPE>
    void ExtractResult( shared_ptr<BaseResult> res,
                        BaseNodeStoreSol * ptStoreSol );

    //! Copy the time derivative of the solution to a solution objet

    //! This method fills the given result objects
    void ExtractDerivResult( shared_ptr<BaseResult> res, UInt deriv );

    //! Copy the linear rhs from the internal vector to a solution object
    template<class TYPE>
    void ExtractRhsResult( shared_ptr<BaseResult> res, 
                           shared_ptr<ResultInfo> eqnResultInfo );

    //! Set containing the types of possible results
    ResultSet availResults_;

    //! Map containing the result types and the results
    ResultMap resultLists_;

    //! map with neighboring volume regions for surface results
    std::map<shared_ptr<BaseResult>,RegionIdType> surfNeighborRegions_;
    //@}
    
    // -----------------------------------------------------------------------
    // Adaptivity
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to adaptivity

    //! object with methods for error estimation
    SpaceErrorEstimator * ptError_;

    //! array where  we store the number of refinement for the element
    StdVector<Double> markingElems_;

    Vector<Double> errorMap_;  //!< array with error map
    Double tolSpaceErr_;       //!< tolerance
    //@}
    

    // -----------------------------------------------------------------------
    // Miscellaneous paramters
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellanous paramters

    //! Identifier of the PDE which is used in the algebraic system
    FeFctIdType  pdeId_;
  
    //! flag for direct coupling
    bool isDirectCoupled_;

    //! flag indicating if Init() was already called
    bool isInitialized_;

    //! flag indicating if penalty method is used
    bool usePenalty_;

    //! maximum order of partial derivatives w.r.t. time
    UInt maxTimeDerivOrder_;

    //! PDEMemento
    shared_ptr<PDEMemento> memento_;

    //! usage type of memento
    bool mementoAsDirichlet_;

    //! Handle for MathParser object
    MathParser::HandleType mHandle_;

    //! map for storing bilinear forms needed for postprocessing
    std::map< RegionIdType, std::map< std::string, BaseForm* > > pdeBilinearForms_;


    //@}
  private:
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SinglePDE
  //! 
  //! \purpose 
  //! This class serves as base class for all single field problems, 
  //! like electrostatic,  acoustic, mechanic and others.
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
