#ifndef FILE_SINGLEPDE
#define FILE_SINGLEPDE

#include "PDE/StdPDE.hh"

#include <list>
#include <map>

#include "Utils/mathParser/mathParser.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Domain/BCs.hh"
#include "Domain/Results/BaseResults.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"

namespace CoupledField
{
  // forward class declaration
  class SpaceErrorEstimator;
  class BasePairCoupling;
  class DirectCoupledPDE;
  class Assemble;
  class BaseForm;
  class BiLinearForm;
  class BaseBDBInt;
  class PDEMemento;
  class ResultFunctor;
  class BaseFieldFunctor;
  class LinearFormContext;
  class CoefFunction;
  class CoefFunctionMulti;
  
  //! Base class for all kinds of single field problems.

  class SinglePDE : public StdPDE {
  
  public:

    // friend declaration
    friend class BasePairCoupling;
    friend class DirectCoupledPDE;

    
    typedef StdVector<shared_ptr<BaseResult> > ResultList;
    typedef std::map<shared_ptr<ResultInfo> , ResultList > ResultMap;

    //! Initialize PDEs 
    //! @param base pointer to InfoNode of this PDE */
    virtual void Init( UInt sequenceStep, PtrParamNode base = PtrParamNode() ); 

    //! Finalize initialization of PDEs. Either called in SinglePDE::Init() if
    //! directCoupled_ == false, otherwise called from DirectCoupledPDE::Init()
    void FinalizeInit( );

    // ---------------------- ***** --------------------------------

    //! Destructor
    virtual ~SinglePDE();
  
  
    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! Initalize PDE coupling (only done once)
    virtual void InitCoupling(PDECoupling * Coupling) = 0;

    //! Reset coupling counters and data (done after each timestep)
    virtual void ResetCoupling();

    //! Fill in input coupling terms
    virtual void CalcInputCoupling();
  
  
    //! Calculate coupling terms
    virtual void CalcOutputCoupling() = 0;

  
    // ======================================================
    // GET /SET  METHODS
    // ======================================================

    /** return sub type. The string is stored internally any we need to convert. :(
     * @return if StdPDE::subType_ is not set we return NO_TENSOR  */
    SubTensorType GetSubTensorType() const;

    //! Set Direct coupling information
    virtual void SetDirectCoupling();

    //! set boundary condition
    void SetBCs();

    //! set special PDE dependent boundary conditions
    virtual void SetSpecialBCs(){ return; }

    //! Update PDE due to updated step in multistep solution strategy
    virtual void UpdateToSolStrategy();
    
    /** Write general defines (BCs, loads, etc.) to info.xml.
     * Note, that only the current state is (over) written! */
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

    //! Obtain coefficient function of given type
    PtrCoefFct GetCoefFct( SolutionType solType );
    
    //! Read general external field information from given xml node
    //! The node has to contain either a values tag, a number of comp tags or
    //! a grid node
    //! \param[in] name The name of the entityList the Field should be applied to
    //! \param[in] valueNode The xml node of the user parameters
    //! \param[in] compNames Names of the components (vector, tensor)
    //! \param[in] type Type of CoefFunction to be read in (scalar, vector, tensor)
    //! \param[in] isComplex Indicates if we need to account for Complex results
    //! \param[out] coef The generated coefficient function
    //! \param[out] definedDofs Set containing all defined dofs in case of a
    //!             vector-valued quantity.
    void ReadUserFieldValues( const std::string name,
                              PtrParamNode valueNode,
                              const StdVector<std::string>& compNames,
                              ResultInfo::EntryType type,
                              bool isComplex,
                              PtrCoefFct & coef,
                              std::set<UInt>& definedDofs );

  protected:

    //! Constructor
    /*!
      \param aptgrid pointer to grid
    */
    SinglePDE( Grid *aptgrid, PtrParamNode );

    //! private copy constructor
    SinglePDE & operator= (const StdPDE & myPDE) {
      EXCEPTION( "Not implemented" );

      // For compiler
      return *this;
      }

    // ======================================================
    // INITIALIZATION METHODS
    // ======================================================

    //! Define primary results
    
    //! Initialize the primary results, i.e. the results corresponding to the
    //! primary variables and their unknowns.
    virtual void DefinePrimaryResults() { };
    
    //! Define post-processing results
    
    //! This method defines the post-processing results, which are computed
    //! mostly using the bilienarform of the main problem.
    //! Thus, this method gets called after the integrators are defined
    virtual void DefinePostProcResults() {};

    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();


    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( )=0;
    
    //! define surface integrators needed for this pde
    virtual void DefineSurfaceIntegrators( )=0;

    //! Define all RHS linearforms for load / excitation 
    virtual void DefineRhsLoadIntegrators( ) {}


    //! Read single RHS excitation
    
    //! This method reads an xml element for a general RHS excitation and
    //! returns the entityList and CoefFunction. 
    //! \param elemName Name of ParamNode within <bcsAndLoads> to be read 
    //! \param compNames Names of the components (vector, tensor)
    //! \param type Type of CoefFunction to be read in (scalar, vector, tensor)
    //! \param isComplex Denotes  if a complex valued coef-function is to be 
    //!                  generated
    //! \param entities Vector of entityLists of the boundary condition
    //! \param coef Vector of coefficients function for the values
    void ReadRhsExcitation( const std::string& elemName, 
                            const StdVector<std::string>& compNames,
                            ResultInfo::EntryType type,
                            bool isComplex,
                            StdVector<shared_ptr<EntityList> >& entities, 
                            StdVector<PtrCoefFct>& coef );



    //! Read results information for interpolation of continuous fields
    virtual void ReadFieldResults();
    
    // =======================================================================
    //   INTERPOLATION OF FIELD VARIABLES
    // =======================================================================
    
    //! Helper struct for interpolating field variables at arbitrary points
    struct FieldAtPoints {

      //! Physical Quantity
      shared_ptr<ResultInfo> resultInfo;
     
      //! Filename where points get written to
      std::string fileName;

      //! Format output file as CSV (comma separated values)
      bool csv;

      //! Delimiter for CSV fields
      char delim;

      //! Vector with elements 
      StdVector<const Elem*> elems;
      
      //! Vector with local points
      StdVector<LocPoint> locPoints;

      //! Flux values for each point
      SingleVector* field;
    };
    
    //! List of fields to be interpolated
    StdVector<FieldAtPoints> fields_;
    
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
    
    //! write results in file
    void WriteResultsInFile( const UInt kstep, 
                             const Double actTimeFreq );
    
    //! incorporate information of memento object, if set
    virtual void IncorporateMemento();

    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! Define the time FeFunctions for this PDE according to the
    //! definition in the XML file
    virtual void DefineFeFunctions();

    //! Create FeSpaces according to formulation
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode ) = 0;

    //@{
    
  public:

    //! Class defining data needed for defining Rayleigh damping
    struct RaylDampingData {
      
      //! Damping parameters used for MASS and STIFFNESS integrator
      std::string alpha, beta;
      
      //! Ratio for calculation of deltaF
      Double ratioDeltaF;
      
      //! Target frequency, for which alpha and beta should get computed
      Double freq;
      
      //! Use damping adjustment to achieve constant tanDelta
      bool adjustDamping;
    };
    
    
  protected:
    
     
    //@}

    // ======================================================
    // DATA SECTION
    // ======================================================

    // reads in the PML data
    void ReadDataPML( std::string& typePML, Matrix<Double>& inner, 
		                  Double& dampPML, 
		                  std::string& coordSysId,
		                  PtrParamNode actNode);

    //! computes the PML layer dimensions
    void GetPMLLayerData( Matrix<Double>& inner, 
                          Matrix<Double>& outer,
			                    RegionIdType regionId,
			                    std::string& coordSysId );
  
    // -----------------------------------------------------------------------
    // Storing information
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to storing information
    
    //! Define a field result
    void DefineFieldResult( PtrCoefFct coef, shared_ptr<ResultInfo> res );
    
    //! Define result based on the time derivative of the main results
    void DefineTimeDerivResult( SolutionType derivSolType,
                                UInt timeDerivOrder,
                                SolutionType primSolType );
    
    //! Map containing the result types and the results
    ResultMap resultLists_;

    //@}
    
    // -----------------------------------------------------------------------
    // Miscellaneous parameters
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous parameters

    //! flag for direct coupling
    bool isDirectCoupled_;

    //! flag indicating if Init() was already called
    bool isInitialized_;

    //! maximum order of partial derivatives w.r.t. time
    UInt maxTimeDerivOrder_;

    //! PDEMemento
    shared_ptr<PDEMemento> memento_;

    //! usage type of memento
    bool mementoAsDirichlet_;

    //! Handle for MathParser object
    MathParser::HandleType mHandle_;

    //! Map for storing the primary BDB integrators of the problem
    
    //! This map stores the primary BDB integrators, which can be used for 
    //! calculating spatial derivatives, fluxes and energy.
    std::map<RegionIdType, BaseBDBInt*> bdbInts_;
    

    //! Map for storing the primary mass integrator of the problem
    
    //! This map stores the primary MASS integrators, which can be used for 
    //! calculating spatial derivatives, fluxes and energy.
    std::map<RegionIdType, BaseBDBInt*> massInts_;

    //! Map for storing "material" parameters as coefficient functions
    
    //! This map holds the coefficient functions for the different material
    //! parameters. As material parameters are typically defined per region,
    //! we store them in a CoefFunctionMulti object
    std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs_;
    
    
    //! Map storing functors for calculating general results
    std::map<SolutionType, shared_ptr<ResultFunctor> > resultFunctors_;
    
    //! Store field coefficient functions
    std::map<SolutionType, PtrCoefFct > fieldCoefs_;
    
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
