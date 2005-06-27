#ifndef FILE_MPCCIEXCH_2003
#define FILE_MPCCIEXCH_2003

#include <Matrix/matrix.hh>
#include <Domain/grid.hh>
#include <General/environment.hh>
#include "PDE/nodeEQN.hh"

namespace CoupledField
{

  //! Class for mesh coupling
  /*! 
    This class handles the subroutines calls concerning MpCCI for coupling the fluid and acoustic computations.
  */

  class MpCCIexch
  {
  public:

    //!AcouFlowNoise constructor
    MpCCIexch(Grid * aptgrid, UInt nNodesSD);

    //!Fluid Structure constructor
    MpCCIexch(Grid * aptgrid, StdVector<RegionIdType> subdoms);

    //!destructor
    virtual ~MpCCIexch();

    //! Reorganizing grid info for MpCCi and hand over to MpCCI.
    //! This method is called during a acoustic-fluid coupled computation
    void PutExchangeGrid2MpCCI(StdVector<RegionIdType> subdoms);

    //! Define a partition for coupling via MpCCI
    void DefMpcciPartition(UInt meshId, UInt partId);

    //! Define Nodes belonging to a MpCCI partition and makes them known in MpCCI
    void DefMpcciNodes(UInt  meshId, UInt partId, UInt nrNodes, UInt* Nodes, NodeEQN & eqnData);

    //! Define Elements belonging to a MpCCI partition and makes them known in MpCCI
    void DefMpcciElements(UInt meshId, UInt partId, NodeEQN & eqnData);

    //! Define the communicator and call CloseSetup
    void FinishMpcciSetup(std::string couplingType);

    /////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////

    //! Performs the coupled computation phase
    void CouplCompPhase(Matrix<Double> & flowdata, UInt timestep);

    //! Receive values of all partitions/subdomains
    void RecvAllPartitions(std::string couplingType);

    //! Get the nodal value of one partition/subdomain
    void GetNodalValOfOnePartition(UInt partId, Vector<Double> & forceData, UInt nrNodesSD, UInt* nodeIds, std::string couplingType);
  
    ////////////////////////////////////////////////////////////////

    //! put values of one partitions/subdomain into the sending queue of MpCCI
    void PutPartition(UInt partId, const Vector<Double>  & displData, UInt nrNodesSD, UInt* nodeIds, Boolean conv);

    //! Sends values of all partitions/subdomains to MpCCI
    void SendAllPartitions();

  private:

    // pointers to objects
    Grid * ptgrid_;           //!< pointer to Grid
    StdVector<RegionIdType> subdoms_;  //!< subdomain-levels belongig to PDE
    ShortInt Dim_;         //!< space dimension of pde  

    //!MpCCI
    UInt MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
    // Integer MpCCI_; //<! if TRUE: coupling via MpCCI to low simulator
    UInt meshId_;
    UInt partId_;
    UInt nNodeIds_;
    Integer *nodeIds_;
    UInt GlobalDim_;
    UInt nElemIds_;
    Integer *elemIds_;
    UInt nElemTypes_;
    Integer *nNodesPerElem_;
    Integer *elemTypes_;
    Integer MpCCIprocess_;
    int ** TOPOLOGYDATA_;
    double ** NODEDATA_;
  };

} // end of namespace
#endif
