#ifndef FILE_FILEREADER_FASTEST_2006
#define FILE_FILEREADER_FASTEST_2006

#include <cplreaderdefs.hh>
#include "../filereader.hh"

namespace CoupledField
{

    class FileReader_FASTEST : public FileReader
    {
    public:

        //! Constructor
        FileReader_FASTEST(const std::string& name,
                           const UInt dim,
                           const UInt numFiles);
    
        //! Deconstructor
        virtual ~FileReader_FASTEST();

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

    protected:

      std::vector<Integer> dataColumns_;
    };

      
} // end of namespace
#endif
