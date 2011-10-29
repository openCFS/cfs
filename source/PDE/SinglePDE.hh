
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
#include "DataInOut/SimInOut/hdf5/simOutputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/simInputHDF5.hh"

namespace CoupledField
{
  // forward class declaration
  class VolForceInt;
  class VolumeSrcInt;
  class SpaceErrorEstimator;
  class BasePairCoupling;
  class DirectCoupledPDE;
  class Assemble;
  class BaseForm;
  class Integrator;
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

    class RegionLoad;

    /** Initialize PDEs 
     * @param base pointer to InfoNode of this PDE */
    virtual void Init( UInt sequenceStep, PtrParamNode base = PtrParamNode() ); 

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

    /** return sub type. The string is stored internally any we need to convert. :(
     * @return if StdPDE::subType_ is not set we return NO_TENSOR  */
    SubTensorType GetSubTensorType() const;

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
    void ReadLoads(ParamNodeList loadNodes, LoadList& out_list);
    
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

    //@{
    //! store the new solution returned by the algebraic system
    //! \param ptSol pointer to solution array
    //! \param size legnth of solution array
    void SaveSolution( SBM_Vector& );
    
    void SavePrevSolution( SBM_Vector& );

    /** Return the native solution type, MECH_DISPLACEMENT, ... */
    virtual SolutionType GetNativeSolutionType() const { EXCEPTION("not implemented"); }

    /** Return the number of dofs of the native solution type, MECH_DISPLACEMENT, ...
     * @return 1 or the number of dimensions (displacement) */
    virtual UInt GetNativeDOF() const { EXCEPTION("not implemented"); }


    //@}

    /** do the actual reading of loads, this is also called from optimization 
     * @param bcNode paramnode that has "regionLoad" nodes as children 
     * @param pressSurf StdVector containing the RegionLoads */
    void ReadRegionLoadsFromXML(PtrParamNode bcNode, std::map<RegionIdType, RegionLoad>& regionLoads_);
//    
//    /** reread the results, this is needed in transient optimization, so that the results are re registered for the new multisequencestep 
//     */
//    void ReReadResults(){
//      ReadStoreResults();
//      ReadSpecialResults();
//    }

//    /** Calculate the element gradient matrix for all nodal positions of the matrix where only the
//     * dof component of the gradient is stored. This matrix just needs to be multiplied with the element
//     * solution vector.
//     * @param elem a volume element
//     * @param wanted which nodes of the element shall be considered. elem->connect for all!
//     *        Only the nodenumbers are of interest, not the ordering
//     * @param dof which dof of the gradient to use (1 based)
//     * @param fctType do we use non LAGRANGE?
//     * @param q_mat will be resized to (elem nodes)^2 with rows for all wanted nodes */
//    template <class TYPE>
//    void CalcElemGradMatrix(const Elem* elem, const StdVector<UInt>& wanted, UInt dof,
//                               shared_ptr<AnsatzFct> fctType, Matrix<TYPE>& q_mat);
//
//    /** Calculate the gradient of all nodes of a single element based on the elements ansatz functions.
//     * Strictly the gradient is only defined in the interior of the element, so for multiple elements
//     * you need to average and know what you are doing!
//     * @param dof mechDisplacement has x=1,y=2(,z=3) solution per node. Select the node. 1 for scalar results (potential)
//     * @param grad_out will resized to the number of nodes of the element, will get the gradients of size dim.
//     * @param elem_data the element solution of NULL if the current OLAS solution vector shall be used!*/
//    template <class TYPE>
//    void CalcGradNodeSolution(const Elem* elem, shared_ptr<AnsatzFct> fctType, UInt dof,
//                                 StdVector<Vector<TYPE> >& grad_out, Vector<TYPE>* elem_data = NULL);
//
//
//    /** Calculate the nodal gradients for several elements, does averaging.
//     *  @param list needs to be a ELEM_LIST
//     *  @param dof mechDisplacement has x=1,y=2(,z=3) solution per node. Select the node. 1 for scalar results (potential)
//     *  @param nodal_grad by the *global* node numbers of the elements the averaged gradients are written.
//     *  @param counter is internally needed but as parameter can be reused. Is resized to global number of nodes */
//    template <class TYPE>
//    void CalcGradNodeSolution(shared_ptr<EntityList> list, UInt dof, StdVector<Vector<TYPE> >& nodal_grad, StdVector<UInt>& counter);

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
      ;}

   
    // ======================================================
    // INITIALIZATION METHODS
    // ======================================================

    /** define all computable results. */
    virtual void DefineAvailResults() { };

    /** The gradSolution is common for (almost) all PDEs.
     * @param dof for mech displacement, 1 = grad of x displacement, ... */
    template<class TYPE>
    void CalcGradSolution(shared_ptr<BaseResult> res, UInt dof = 1);

    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();


    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( )=0;

    /** Trigger calculation of results. */
    virtual void CalcResults( shared_ptr<BaseResult> result ) { };
    
    //! Read results information for interpolation of continuous fields
    virtual void ReadFieldResults();
    
    //! Calculate field variables at arbitrary points
    virtual void CalcField( SolutionType solType, StdVector<const Elem*>& elems,
                            StdVector<LocPoint>& points, SingleVector& values ) {}

    // =======================================================================
    //   INTERPOLATION OF FIELD VARIABLES
    // =======================================================================
    
    //! Helper struct for interpolating field variables at arbitrary points
    struct FieldAtPoints {

      //! Physical Quantity
      shared_ptr<ResultInfo> resultInfo;
     
      //! Filename where points get written to
      std::string fileName;

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
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode ) = 0;

    //@{
    //! Save load part of RHS to private variable
    void SaveRHS( SBM_Vector& );
    
  public:

    //    //! Class defining data needed for region loads
    class RegionLoad {

    public:

      //! Constructor
      RegionLoad( UInt dim, bool isaxi );

      void ToInfo(PtrParamNode in) const;
      
      //! Returns the RHS-integrator
      void GetIntegrator();

      //! Returns the RHS-integrator for scalar source
      void GetSrcScalarIntegrator();
      
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
    
    //! List of region loads
    std::map<RegionIdType, RegionLoad> regionLoads_;
     
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

    //! check if subdomain is a coupled piezo subdomain with micro-piezo-model
    bool IsRegionMicroPiezo( std::string regionName );
    
    //! check if PDE is a coupled piezo subdomain with  micro-piezo-model
    bool BelongsPDE2MicroPiezo();    
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
    
    //! Copy the solution from a vector to a BaseNodeStoreSol object
    template<class TYPE>
    void ExtractResult( shared_ptr<BaseResult> toBaseStore, \
                                   const Vector<TYPE>& fromVec );
    //! Copy the solution from a nodeStoresolution object to a solution object
    template<class TYPE>
    void ExtractResult( shared_ptr<BaseResult> toBaseResult,
                        BaseNodeStoreSol* ptStoreSol );

    //! Copy the time derivative of the solution to a solution objet

    //! This method fills the given result objects
    void ExtractDerivResult( shared_ptr<BaseResult> res, DERIVType derivType );

    //! Copy the linear rhs from the internal vector to a solution object
    template<class TYPE>
    void ExtractRhsResult( shared_ptr<BaseResult> res, 
                           shared_ptr<ResultInfo> eqnResultInfo );
    //! Copy the solution from a solution object to a nodeStoresolution object
    //! So the oposite way if ExtractResult.
    template<class TYPE>
    void InsertResult( Vector<TYPE> &resVec,
                       shared_ptr<BaseResult> res );

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
    /** write out results for restart
     * @param outFile the file to write to
     * @param outResults the results to write
     */
    inline void writeOutTimeStep(shared_ptr<SimOutput>& outFile, \
        StdVector<shared_ptr<BaseResult> >& outResults);
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
