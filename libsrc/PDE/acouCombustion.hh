#ifndef FILE_ACOUCOMBUSTION_2006
#define FILE_ACOUCOMBUSTION_2006

#include "basePDE.hh"
#include "acousticPDE.hh"

namespace CoupledField
{
  typedef enum {NOSRC,  REYNOLDSTRESS, MOMENTUM, CHEMICAL, DENSITY2DERIV, 
                HEATRELEASE, GASCONSTDERIV2, SHEARTERM, REYNOLDSTRESS_MOMENTUM_CHEMICAL, 
                HEATRELEASE_SHEARTERM, GASCONSTDERIV2_SHEARTERM,
                HEATRELEASE_GASCONSTDERIV2, HEATRELEASE_GASCONSTDERIV2_SHEARTERM} SrcType;

  //! Class for acoustic equation
  class AcouCombustionNoise: public AcousticPDE
  {
  public:

    //!  
    AcouCombustionNoise(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut);

    //!
    virtual ~AcouCombustionNoise();

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //! put nodal values to RHS
    void ComputeRHS(const Double atime);

    //! Reads at every time the flowdatafile from CFD code
    void ReadFlowData(const char * aname);

  private:

    //get the correct source term
    void GetSrcTerm(Double& val, UInt idx);

    //CFD Data
    Double** dataCFD_;

    //CFD Data
    Integer* nodesCFD_;

    //! type of source term
    SrcType srcType_;

    // true, if harmonic analysis
    bool isHarmonic_;

    //!< name of subdomain to be coupled with CFD
    StdVector<RegionIdType> couplSubDomId_; 
    
    //! number of nodes in coupling region
    UInt numNodesInCoupledRegion_;

    //! number of source terms in CFD-file
    UInt numDataInCFD_;

    //! CFD-file name 
    std::string filenameCFD_;

    //! index for variable speed of sound
    UInt varSpeedOfSoundIdx_;

    //! index for source term: Reynold stress
    UInt srcReynoldStressIdx_;

    //! index for source term: momentum fluctuation
    UInt srcMomentumIdx_;

    //! index for source term: chemical
    UInt srcChemicalIdx_;

    //! index for source term: chemical
    UInt srcDensityDeriv2Idx_;

    //! index for source term: heat release
    UInt srcHeatReleaseIdx_;

    //! index for source term: gas const. 2nd time derivative
    UInt srcGasConstDeriv2Idx_;

    //! index for source term: shear term
    UInt srcShearTermIdx_;

    //! time step counter
    UInt timeStep_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class AcouCombustionNoise
  //! 
  //! \purpose 
  //! This class is derived from class AcousticPDE. It is used for solving 
  //! the acoustic analogy inhomogeneous acoustic equation for combustion
  //! noise on one time step. We define the right hand side by setting the
  //! sound sources using the results from the external CFD computation, set 
  //! rules for assembling global system matrix according to weak form of
  //! this PDE, define right hand side and set first order absorbing boundary
  //! conditions. Then we call one of methods of LinSystem for solving linear
  //! system. On the last step we calculate first and second derivatives of
  //! the solution.
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
