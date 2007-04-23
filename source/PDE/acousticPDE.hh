// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include "SinglePDE.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODESolve/ODESolver_RKF45.hh"
#include "Forms/bubbleDampInt.hh"
#include "Forms/bubbleStiffInt.hh"
#include "Utils/mathParser/mathParser.hh"
 
namespace CoupledField {

  //! Class for acoustic equation (no adaptivity)


  class AcousticPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    AcousticPDE(Grid* aGrid, ParamNode* paramNode );

    //! Destructor
    virtual ~AcousticPDE(){};


   //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> result );


#ifdef ADAPTGRID
    //! test error of computation
    virtual bool TestError(const UInt level);
#endif


    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);
  
    //! calculate coupling terms
    void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output);
  
    //! calculate the vector of coupling forces to the mechanical PDE
    void CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                              StdVector<UInt> & couplingNodes,
                              Vector<Double>& elemCouplingSols,
                              UInt couplingdof );

    //! calculate the vector of coupling surface nodes to the nrbcPDE  
    void CalcNRBCCouplingRHS( StdVector<Elem*> * couplingElems, 
                              StdVector<UInt> & couplingNodes,
                              Vector<Double>& elemCouplingSols,
                              UInt couplingdof );
    
    //! calculate the heat source term for heatConduction PDE
    template <class TYPE>
    void CalcHeatCouplingRHS( Vector<Double> & energy, 
                              StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                              UInt actCoupling, UInt numCouplingNodes );

    //! 
    void SetMechanicCoupling() {
      isMechCoupled_ = true;
    }
    //! 
    void SetNrbcCoupling() {
      isNrbcCoupled_ = true;
    }

    //! returns formulation (pressure/scalar potential)
    SolutionType GetFormulation() {
      return formulation_;
    }

  protected:

    // ========================
    // initialization
    // ========================    

    //! Init the time stepping
    void InitTimeStepping();

    //! Define availabe result types
    void DefineAvailResults();

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! Read special boundary conditions (here: bubble information)
    void ReadSpecialBCs();
    
    //! Initialize NonLinearities
    void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();
  
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //     //! define the algbraic system
    //     void DefineAlgSys();

    // ========================
    // Postprocessing
    // ========================

    
    //! calculate Force acting on specified surface elements
    template <class TYPE> 
    void CalcForce( shared_ptr<BaseResult> vals );
    
    //! calculate pressure from acoustic potential
    template <class TYPE>
    void CalcElemPressure( shared_ptr<BaseResult> vals );
    
    //! Calculate acoustic power
    template <class TYPE>
    void CalcAcouPower( shared_ptr<BaseResult> vals );
    
    //! Calclulate acoustic intensity
    template <class TYPE>
    void CalcAcouIntensity( shared_ptr<BaseResult> vals );

   //! calculate element mean pressure and derivative for bubble PDE
    void CalcBubblePressure( StdVector<Elem*>& couplingElems,
                             Vector<Double>& elemCouplingSols,
                             SolutionType solType );


    // ========================
    // set solution information
    // ========================
    //! Returns the resultInfo related to the specified solutionType

    bool isMechCoupled_; //!< indicator for mechanic coupling
    
    SolutionType formulation_; //!< variable in which PDE is formulated


    //! indicator for mechanic coupling
    bool isNrbcCoupled_;    

    StdVector<RegionIdType> absBCs_; //!< subdomains, which form absorbing BCs

    bool absorbingBCs_; //!< switch for absorbing BCs     

    //! indicator for bubble coupling
    bool isBubbleCoupled_;   
    
    //bool fracDamping_; //!< switch indicating use of fractional damping
    
    // ========================
    // time stepping
    // ========================
    // solving of nonlinear acoustics
    Vector<Double> RhsLinVal_;


    // ========================
    // coupling
    // ========================
    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; 


    // ========================
    // Postprocessing results
    // ========================

    //! right hand side vector
    NodeStoreSol<Double> rhs_; 

    //! right hand side vector if we have a multi-dim RHS (e.g. vel)
    NodeStoreSol<Double> rhs2_;
    
    //!  nodal acoustic sources ( combustion noise )
    NodeStoreSol<Double> rhsNodalSrc_;

    //! variable speed of sound( combustion noise )
    NodeStoreSol<Double> speedOfSound_; 

    //! Attribute describing model for bubble dynamics
    BubbleDynType bubbleDynType_;

    //! bubbledensity
    Double bubbleDensity_;

    bool plotRHS_; //!< Flag for saving of rhs for output
    bool plotRHSVel_; //!< Flag for saving of rhs as a vector field
    
    bool saveNodalSourcesRHS_;  //!< Flag for saving nodal acoustic sources

    // variable speed of sound
    bool variableSpeedOfSoundCN_;

    //! vector containing regionIds of non-conforming interfaces
    StdVector<RegionIdType> ncIFaces_;

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the acoustics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the AcousticPDE class that is set, if the
    //! keyword is specified.\n
    //! \todo Specification of ReadStoreResults for AcousticPDE!!!
    void ReadStoreResults();

    //! read flow data
    void ReadFlowData();

    //! information for damoing layer
    void ReadDataDampLayer(std::string& dampingTypePML, 
			   Vector<Double>& mPoint, 
			   Double& dampFactor, 
			   Double& dampFactorMax, 
			   Double& startRadius, 
			   Double& endRadius, 
			   ParamNode * actNode );

    

    //! 

    std::map<RegionIdType, BubbleDampInt*> bubbleDampIntMap_; 
    std::map<RegionIdType, BubbleStiffInt*> bubbleStiffIntMap_; 
  
    //! map storing for each region the related flowData node
    std::map<RegionIdType, ParamNode*> regionFlowNodes_;

    //! Handle for MathParser object
    MathParser::HandleType mHandle_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class AcousticPDE
  //! 
  //! \purpose 
  //! This class is derived from class SinglePDE.
  //! It is used for solving acoustic equation on one time step.  
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
