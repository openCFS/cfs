#ifndef FILE_MPCCIEXCH_2003
#define FILE_MPCCIEXCH_2003

#include <Matrix/matrix.hh>
#include <Domain/grid.hh>
#include <General/environment.hh>


namespace CoupledField
{

  //! Class for mesh coupling
  /*! 
    This class handles the subroutines calls concerning MpCCI for coupling the fluid and acoustic computations.
  */

class MpCCIexch
{
public:

  //!
  MpCCIexch(Grid * aptgrid, Integer nNodesSD);

  //!
  virtual ~MpCCIexch();

  //! Reorganizing grid info for MpCCi and hand over to MpCCI
  void PutExchangeGrid2MpCCI(std::vector<std::string> subdoms);

  //! Performs the coupled computation phase
  void CouplCompPhase(Matrix<Double> & flowdata, Integer timestep);

private:

  // pointers to objects
  Grid * ptgrid_;           //!< pointer to Grid
  std::vector<std::string> subdoms_;  //!< subdomain-levels belongig to PDE
  Integer actlevel_; //! actual level
  ShortInt Dim_;         //!< space dimension of pde  

  //!MpCCI
  Integer MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
  // Integer MpCCI_; //<! if TRUE: coupling via MpCCI to low simulator
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
};

} // end of namespace
#endif
