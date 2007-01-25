#ifndef FILE_ACOUFLOWNOISE_2003
#define FILE_ACOUFLOWNOISE_2003

#include "basePDE.hh"
#include "acousticPDE.hh"
#include <MpCCIcpl/MpCCIexch.hh>
//#include <MpCCIcpl/MpCCIexch.hh>

namespace CoupledField
{

  //! Class for acoustic equation

  class AcouFlowNoise: public AcousticPDE
  {
  public:

    //!  
    AcouFlowNoise(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut);

    //!
    virtual ~AcouFlowNoise();

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //!
    void ComputeRHS(const Double atime);

    //! Writes out fine RHS sources computed using VortexAnalytical() in
    //! transient files
    void WriteOutVortexFineSources(UInt timestep);

    //! Computes RHS using VortexAnalytical() to get source values
    void ComputeRHSwithVortexSource(const Double atime);

    //! Reads at every time the flowdatafile from the Fluid's Computation
    void ReadFlowData(const char * aname, const UInt timestep,
                      Matrix<Double> &nodedata );

    void VortexAnalytical(Double & press, Vector<Double>& dTij_di, const Double x,
                                       const Double y, const Double t, 
                        const UInt outType);
    

    Double JJ2(const Double x);
    Double YY2(const Double x);
    

  private:


    UInt arg_rhs_; //!< function for RHS


    //Flow Data
    Matrix<Double> flowdata_;

    //!type of FlowData
    bool nodalSrc_;
    //!flag when using vortex source
    bool vortexSrc_;
    //!flag to verify if analysis harmonic
    bool isHarmonic_;
     //!< name of subdomain to be coupled with MpCCI or where vortex is applied
     StdVector<RegionIdType> couplSubDomId_; 
    //!< mapping of linear nodes only (3D MpCCI allows only linear elements)
     StdVector<UInt> mapSD_onlyLinNodes_; 
    //<! mapping containing all nodes of coupled region
     StdVector<UInt> mapSD_allNodes_;
    
    bool  writeGridFile_; //!<flags to write grid with coupled vals in file
    bool  writeSrcFileperTS_; //!<flags to write coarse srcs in time step files
    bool pressFormul_;//!<Flag for recognizing sources from pressure formulation
    
    //!Objects for topology files
    std::ofstream * outelemfile_;
    std::ofstream * outnodefile_;
    //!Object to file for storing src in time step files (NrFiles=NrTimeSteps)
    std::ofstream * outsrcfile_;

    UInt vortexFlag_;//! To choose output for vortex sound computation


    //!MpCCI
#ifdef MpCCI
    MpCCIexch * ptMpCCIexch_;
    Integer MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
#endif

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class AcouFlowNoise
  //! 
  //! \purpose 
  //! This class is derived from class AcousticPDE. It is used for solving 
  //! the acoustic analogy inhomogeneous acoustic equation for flow induced
  //! noise on one time step. We define the right hand side by setting the
  //! sound sources using the results from the external fluid computation, set 
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
