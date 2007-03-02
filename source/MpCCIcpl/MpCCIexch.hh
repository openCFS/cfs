#ifndef FILE_MPCCIEXCH_2003
#define FILE_MPCCIEXCH_2003

#include <def_use_mpcci.hh>

#include <Matrix/matrix.hh>
#include <Domain/grid.hh>
#include <General/environment.hh>

namespace CoupledField
{

  //! forward class declaration
  class ParamNode;

  //! Class for mesh coupling
  /*! 
    This class handles the subroutines calls concerning MpCCI for coupling the fluid and acoustic computations.
  */

class MpCCIexch
{
public:

  //!FlowNoise Constructor
  MpCCIexch(Grid * aptgrid, StdVector<UInt> & mapSD);

  //!Fluid Structure constructor
  MpCCIexch(Grid * aptgrid, StdVector<RegionIdType> subdoms);

  //!
  virtual ~MpCCIexch();

  //! Reorganizing grid info for MpCCi and hand over to MpCCI
  void PutExchangeGrid2MpCCI(StdVector<RegionIdType> subdoms);

  //! Define a partition for coupling via MpCCI
  void DefMpcciPartition(UInt meshId, UInt partId);
  
  //! Define Nodes belonging to a MpCCI partition and makes them known in MpCCI
  void DefMpcciNodes(UInt  meshId, UInt partId, UInt nrNodes, UInt* Nodes, shared_ptr<EqnMap> eqnMap);
  
  //! Define Elements belonging to a MpCCI partition and makes them known in MpCCI
  void DefMpcciElements(UInt meshId, UInt partId, shared_ptr<EqnMap> eqnMap);
  
  //! Define the communicator and call CloseSetup
  void FinishMpcciSetup(std::string couplingType);
  
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////

  //! Performs the coupled computation phase
  void CouplCompPhase(Matrix<Double> & flowdata, Double acttime);
  
  //! Receive values of all partitions/subdomains
  void RecvAllPartitions(std::string couplingType);
  
  //! Get the nodal value of one partition/subdomain
  void GetNodalValOfOnePartition(UInt partId, Vector<Double> & forceData, UInt nrNodesSD, UInt* nodeIds, std::string couplingType);
  
  ////////////////////////////////////////////////////////////////
  
  //! put values of one partitions/subdomain into the sending queue of MpCCI
  void PutPartition(UInt partId, const Vector<Double>  & displData, UInt nrNodesSD, UInt* nodeIds, bool conv);
  
  //! Sends values of all partitions/subdomains to MpCCI
  void SendAllPartitions();

private:

  // pointers to objects
  ParamNode * myParam_; //!< pointer to "MpCCI-flownoise"-element
  Grid * ptgrid_;           //!< pointer to Grid
  StdVector<RegionIdType> subdoms_;  //!< subdomain-levels belongig to PDE
  ShortInt Dim_;         //!< space dimension of pde  
  StdVector<UInt> mapSD_;//!< local vector containing coupled domain
  StdVector<RegionIdType> couplSubDomId_; 
  bool  writeGridFile_; //!<flags to write grid with coupled vals in file
  bool  writeSrcFileperTS_; //!<flags to write coarse srcs in time step files
  bool  writeSrcFileperNode_; //!<flags to write coarse srcs in nodal files
  //!Objects for topology files
  std::ofstream * outelemfile_;
  std::ofstream * outnodefile_;
  UInt numSteps_; //!< to retrieve number of time steps from xml file
  Double firstdt_;//!< to retrieve time step size from xml file
  Double ** nodalSrcMat_; //!< to store srcs in a nodal times numTimeSteps matrix
  //!Object to file for storing src in time step files (NrFiles=NrTimeSteps)
  std::ofstream * outsrcfile_;

  //!MpCCI
  UInt MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
  // Integer MpCCI_; //<! if true: coupling via MpCCI to low simulator
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
