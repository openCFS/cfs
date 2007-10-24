#ifndef FILE_MPCCIEXCH_2006
#define FILE_MPCCIEXCH_2006

#include "filereader.hh"

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
    MpCCIexch(int argc, char *argv[], FileReader * ptFileReader);


  //!
  virtual ~MpCCIexch();

  //! Reorganizing grid info for MpCCi and hand over to MpCCI
  void PutExchangeGrid2MpCCI();

  //! Performs the coupled computation phase
  void Couple();

private:

    void CalculateAcouSrcs(const int partitionIdx,
                           const std::vector<double>& velocity,
                           // const std::vector<double>& pressure,
                           std::vector<double>& acouSrc);
    

  FileReader *  ptFileReader_;
  std::map<int, std::vector<Realtype> > NodalCoords_;
  std::map<int, std::vector<int> > Topology_;
  
  std::map<int, std::vector<int> > numNodesPerElem_;
  std::map<int, std::vector<int> > elemTypes_;

  std::map<int, std::map<int, int> > nodesInPartition_;
  
#if 0
  //!MpCCI
  int meshId_;
  int GlobalDim_;
  int MpCCIprocess_;
#endif

};

} // end of namespace
#endif
