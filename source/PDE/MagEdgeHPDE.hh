/** DOCUMENTATION HEADER
 * This PDE is used to calculate the magnetostatic field via
 * the penalized h-formulation that originates from the mixed ha-formulation [1]:
 * 
 * (I):  (curlh,a') = (j,a')              for all a',
 * (II): (mu h,h') - (a,curlh') = 0       for all h'.
 * 
 * Eliminate the mag. vector potential a yields the penalized h-formulation
 * 
 * (I): (mu h,h') + (rho_art curlh,curlh') = (rho_art j, curlh') for all h'.
 * 
 * This is the equation that is solved using this class.
 * The equation can involve nonlinear materials, which is realized by the
 * energy-based hysteresis model [2], which is able to perform anhystertic and
 * hysteretic materials.
 * 
 * REFERENCES:
 * [1]  C. Gaier and H. Haas, "Edge element analysis of magnetostatic and transient eddy current fields using
 *      H/spl I.oarr/-formulations," in IEEE Transactions on Magnetics, vol. 31, no. 3, pp. 1460-1463, May 1995, 
 *      doi: 10.1109/20.376304.
 * [2] L. Prigozhin, V. Sokolovsky, J. W. Barrett and S. E. Zirka, "On the Energy-Based Variational Model for 
 *     Vector Magnetic Hysteresis," in IEEE Transactions on Magnetics, vol. 52, no. 12, pp. 1-11, Dec. 2016, 
 *     Art no. 7301211, doi: 10.1109/TMAG.2016.2599143.
 */
#ifndef FILE_MAGNETICEDGEHPDE
#define FILE_MAGNETICEDGEHPDE

#include <map>
#include "MagBasePDE.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>


namespace CoupledField
{

  // Forward declaration of classes
  class CurlCurlEdgeInt;
  class nLinCurlCurlEdgeInt;

  //! Class for 3D magnetics H-formulation using edge elements
  class MagEdgeHPDE : public MagBasePDE 
  {
  public:

    //! Constructor
    MagEdgeHPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeHPDE();

//    //! Get mehtod for specific coils. Needed e.g. by the SinglePDE for
//    //! specifying coil results.
//    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return MAG_FIELD_INTENSITY; }

    //stores the flux for hystersis and nonlinear models
    std::map<RegionIdType, shared_ptr<CoefFunctionMulti> > nlFluxCoefm_;
    std::map<RegionIdType, shared_ptr<CoefFunctionMulti> > nlScalCoefm_;

    //! Initialize NonLinearities
    virtual void InitNonLin();
  protected:
    
    LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "" );

    // define all integrators: bilinear forms and linearforms (left hand side and right hand side)
    void DefineIntegrators();
    //! define all (bilinearform) integrators needed for this pde
    void DefineStandardIntegrators();
    //! Define all RHS linearforms for load / excitation
    void DefineCoilIntegrators();
    void DefineRhsLoadIntegrators(); 

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators(){}; // not implemented

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){}; // not implemented

    LinearForm* GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation ){return NULL;}; // not implemented
    BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity ){return NULL;}; // not implemented
    void InitMagnetization(){}; // not implemented

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Define available primary result types
    void DefinePrimaryResults();

    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
     
    //! Define available postprocessing results
    void DefinePostProcResults();

    //! Init the time stepping:
    void InitTimeStepping();

    //! map containing the flux density remanence field, e.g., from a previous sequence step or from .xml input file
    std::map<RegionIdType, PtrCoefFct> Brmap_;

    //! Coefficient function, containing the overall permeability
    shared_ptr<CoefFunctionMulti> perm_;
    shared_ptr<CoefFunctionMulti> resistivity_;

    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );

    //! Map containing the remanence (B excitation on RHS)
    std::map<RegionIdType,PtrCoefFct> bRHSRegions_;
    

    //! Set containing all regions with regularized conductivity
    
    //! This set contains all regions, which have no physical conductivity,
    //! but only a non-zero one due to regularization. This must be considered
    //! when calculating eddy currents, which would be "unphysical" in these
    //! regions.
    std::set<RegionIdType> regularizedRegions_;

    std::map<RegionIdType, PtrCoefFct> coilCurrentDens_;

  private:

    //! Use gradient fields in shape functions (Edge elements of second kind)
    bool useGradFields_;

    double dt_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgeHPDE
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
  //! \improve Encapsulate the defintiion of magnets within an own struct
  //! 

#endif

} // end of namespace
#endif