#ifndef FILE_MPCCIEXCH_2006
#define FILE_MPCCIEXCH_2006

#include <General/environment.hh>

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
                           std::vector<double>& flowdata);
    
    void ClearMeshTempFiles(std::string& coordCfgFileName,
                            std::string& coordDatFileName,
                            std::string& connCfgFileName,
                            std::string& connDatFileName,
                            std::string& typesCfgFileName,
                            std::string& typesDatFileName,
                            std::string& importMeshCmd);

    void ClearDatasetTempFiles(std::string& stepNumsFileName,
                               std::string& stepValuesFileName,
                               std::string& resultScriptFileName,
                               std::string& resultDatFileName,
                               std::string& writeResultCmd,
                               std::string& writeAuxCmd);

    
    FileReader *  ptFileReader_;
    std::map<UInt, std::vector<Double> > NodalCoords_;
    std::map<UInt, std::vector<UInt> > Topology_;

    std::map<UInt, std::vector<UInt> > numNodesPerElem_;
    std::map<UInt, std::vector<UInt> > elemTypes_;

    std::map<UInt, std::map<UInt, UInt> > nodesInPartition_;
  
  };

} // end of namespace
#endif
