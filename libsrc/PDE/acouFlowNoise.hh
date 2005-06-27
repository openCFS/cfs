#ifndef FILE_ACOUFLOWNOISE_2003
#define FILE_ACOUFLOWNOISE_2003

#include "basePDE.hh"
#include "acousticPDE.hh"
//#include <MpCCIcpl/MpCCIexch-patchQuad3D.hh>
#include <MpCCIcpl/MpCCIexch.hh>

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

    //! Reads at every time the flowdatafile from the Fluid's Computation
    void ReadFlowData(const char * aname, const UInt timestep,
                      Matrix<Double> &nodedata );



  private:


    UInt arg_rhs_; //!< function for RHS


    //Flow Data
    Matrix<Double> flowdata_;
    //!< name of subdomain to be coupled with MpCCI
    StdVector<RegionIdType> couplSubDomId_; 

    //!type of FlowData
    Boolean nodalSrc_;


    //  BaseNodeStoreSol * rhs_;        //!< For eventual saving of rhs as sol
    
    //!MpCCI
#ifdef MpCCI

    StdVector<UInt> mapSD_;
    StdVector<Integer> mapSD_allNodes_;
    MpCCIexch * ptMpCCIexch_;
    Integer MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
    Integer meshId_;
    Integer partId_;
    Integer nNodeIds_;
    Integer *nodeIds_;
    Integer GlobalDim_;
    Integer nElemIds_;
    Integer *elemIds_;
    Integer nElemTypes_;
    Integer *nNodesPerElem_;
    Integer *elemTypes_;
    Integer MpCCIprocess_;

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
