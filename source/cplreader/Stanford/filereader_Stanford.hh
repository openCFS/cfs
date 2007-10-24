#ifndef FILE_FILEREADER_STANFORD_2006
#define FILE_FILEREADER_STANFORD_2006

#include "filereader.hh"

namespace CoupledField
{

    class FileReader_Stanford : public FileReader
    {
    public:

        //! Constructor
        FileReader_Stanford(const std::string& name,
                           const int dim,
                           const int numFiles);
    
        //! Deconstructor
        virtual ~FileReader_Stanford();

        virtual void Init();


        //! get node coordinates from the corresponding file
        virtual void ReadNodalCoords(std::vector<Realtype> & NODECOORD,
                             const int partitionIdx);


        //! get topology information from the corresponding topology file
        virtual void ReadTopology(std::vector<int> & TOPOLOGYDATA,
                                  std::vector<int> & numNodesPerElem,
                                  std::vector<int> & elemTypes,
                                  const int partitionIdx);

        //! get nodal values from the corresponding fluid datafile
        virtual void ReadNodalValues(std::vector<double> & flowdata,
                             const int partitionIdx,
                             const int timeStepIdx);

    protected:

    };

      
} // end of namespace
#endif
