#ifndef FILE_FILEREADER_OPENFOAM_2007
#define FILE_FILEREADER_OPENFOAM_2007

#include <cplreaderdefs.hh>
#include "../filereader.hh"

class vtkOpenFOAMReader;

namespace CoupledField
{

  class FileReader_OPENFOAM : public FileReader
  {
  public:

    //! Constructor
    FileReader_OPENFOAM(const std::string& name,
                        const UInt dim,
                        const UInt numFiles);
    
    //! Deconstructor
    virtual ~FileReader_OPENFOAM();

    virtual void Init();

    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD,
                                 const UInt partitionIdx);


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                              std::vector<UInt> & numNodesPerElem,
                              std::vector<UInt> & elemTypes,
                              const UInt partitionIdx);

    //! get nodal values from the corresponding fluid datafile
    virtual void ReadNodalValues(std::vector<double> & flowdata,
                                 const UInt partitionIdx,
                                 const UInt timeStepIdx);

    virtual double GetTimeStep(UInt t);
    virtual std::string GetPartitionName(const UInt partitionIdx);

  protected:
      std::vector<Integer> dataColumns_;
      vtkOpenFOAMReader* reader_;
  };

      
} // end of namespace
#endif
