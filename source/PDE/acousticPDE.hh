// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include <map>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/bcs.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/nodestoresol.hh"
 
namespace CoupledField {

  //! Class for acoustic equation (no adaptivity)


class BaseResult;
class EntityList;
class Grid;
class PDECoupling;
struct Elem;

  class AcousticPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    AcousticPDE(Grid* aGrid, PtrParamNode paramNode );

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

    //! calculate the heat source term for heatConduction PDE
    template <class TYPE>
    void CalcHeatCouplingRHS( Vector<Double> & energy, 
                              StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                              UInt actCoupling, UInt numCouplingNodes );

    //! 
    void SetMechanicCoupling() {
      isMechCoupled_ = true;
    }

    //! returns formulation (pressure/scalar potential)
    SolutionType GetFormulation() {
      return formulation_;
    }

    /** Helper for CalcAcouSurfIntensity() and eventually also for other methods to fight the copy and past stuff :(.
     * Used by optimization for acoustic near field optimization.
     * Evaluates for the surface center points the gradient by the help of the neighbor volume element.
     * @param apply_normal if true then the gradient is multiplied by -1 if vol2 is used -> see CalcAcouPower().
     *                     It seems that his corrects actSurfElem->normalSign
     * @param grad_out for every surface element within val the gradient is evaluated. Eventually multiplied by -1 (see apply_normal) */
    template <class TYPE>
    void CalcGradSurfaceElement( shared_ptr<BaseResult> vals, SolutionType solType, bool apply_normal, StdVector<Vector<TYPE> >& grad_out );

    /** @see virtual SinglePDE::GetNativeSolutionType()
     *  @return depending on the formulation: ACOU_POTENTIAL or ACOU_PRESSURE. */
    SolutionType GetNativeSolutionType() const { return formulation_; }

    /** @see virtual SinglePDE::GetNativeDOF() */
    virtual UInt GetNativeDOF() const { return 1; }

  protected:

    // ========================
    // initialization
    // ========================    

    //! Init the time stepping
    void InitTimeStepping();

    //! Define available result types
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

    //     //! define the algebraic system
    //     void DefineAlgSys();

    // ========================
    // Postprocessing
    // ========================

    /** Calculate acoustic power. */
    template <class TYPE>
    void CalcAcouPower( shared_ptr<BaseResult> vals);

    /** Calculate acoustic energy. */
    template <class TYPE>
    void CalcAcouEnergy( shared_ptr<BaseResult> vals);

    //! computes particle velocity out of pressure    
    void  CalcVelFromPressure( shared_ptr<BaseResult> vals );

    //! calculate Force acting on specified surface elements
    template <class TYPE> 
    void CalcForce( shared_ptr<BaseResult> vals );
    
    //! calculate pressure from acoustic potential
    template <class TYPE>
    void CalcElemPressure( shared_ptr<BaseResult> vals );
    

    //! Calculate Kaltenbacher's intensity projected on a surface
    template <class TYPE>
    void CalcAcouSurfIntensity( shared_ptr<BaseResult> vals );
    

    //! Calculate acoustic intensity
    template <class TYPE>
    void CalcAcouIntensity( shared_ptr<BaseResult> vals );

    //! Do some jobs before computation  starts
    virtual void PreparePDE4Computation();


    // ========================
    // set solution information
    // ========================
    //! Returns the resultInfo related to the specified solutionType

    bool isMechCoupled_; //!< indicator for mechanic coupling
    
    SolutionType formulation_; //!< variable in which PDE is formulated

    //! surface elements with absorbing boundary conditions
    StdVector<shared_ptr<EntityList> > absBCs_; 

    //! Stores Rayleigh damping definition for each region
    std::map<RegionIdType, RaylDampingData> regionRaylDamping_;

    //bool fracDamping_; //!< switch indicating use of fractional damping
    
    IdBcList impedanceBCs_;
    
    
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

    // stores acoustic particle velocity; needed for energy computation
    std::map<RegionIdType, Vector<Double> > acouParticleVelocity_;

    //! right hand side vector
    NodeStoreSol<Double> rhs_; 

    //! right hand side vector if we have a multi-dim RHS (e.g. vel)
    NodeStoreSol<Double> rhs2_;
    
    //!  nodal acoustic sources ( combustion noise )
    NodeStoreSol<Double> rhsNodalSrc_;

    //! variable speed of sound( combustion noise )
    NodeStoreSol<Double> speedOfSound_; 

    bool isAPML;  //flag for almost PML formulation
    bool plotRHS_; // Flag for saving of rhs for output
    bool plotRHSVel_; // Flag for saving of rhs as a vector field
    bool justInterpolate_; // Should only the RHS interpolation be performed?
    
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

    //! information for damping layer
    void ReadDataDampLayer(std::string& dampingTypePML, 
			   Vector<Double>& mPoint, 
			   Double& dampFactor, 
			   Double& dampFactorMax, 
			   Double& startRadius, 
			   Double& endRadius, 
			   PtrParamNode actNode );

    //! 

    //! map storing for each region the related flowData node
    std::map<RegionIdType, PtrParamNode> regionFlowNodes_;

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
  //! It is used for solving acoustic equation in transient and harmonic
  //! case, linear and nonlinear.  
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
