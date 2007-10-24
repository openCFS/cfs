#ifndef FILE_FILEREADER_2006
#define FILE_FILEREADER_2006

#include <fstream>
#include <vector>

#define REALTYPE CCI_DOUBLE
typedef double Realtype;

//! Base class for reading topology and nodal data from results of fluid computations 
/*!
  This class contains methods for the reading information 
  from a fluid topology file and from the fluid results files. 
  It reads and couples flow results from a set of files (or set of
  partition files coming from a blocked structured fluid computation).
*/
 
namespace CoupledField
{

    class FileReader
    {
    public:

        //! Constructor
        FileReader(const std::string& name, const int dim, const int numFiles);
    
        //! Deconstructor
        virtual ~FileReader();

        virtual void Init() = 0;

        //! get node coordinates from the corresponding file
        virtual void ReadNodalCoords(std::vector<Realtype> & NODECOORD,
                             const int partitionIdx) = 0;


        //! get topology information from the corresponding topology file
        virtual void ReadTopology(std::vector<int> & TOPOLOGYDATA,
                                  std::vector<int> & numNodesPerElem,
                                  std::vector<int> & elemTypes,
                                  const int partitionIdx) = 0;

        //! get nodal values from the corresponding fluid datafile
        virtual void ReadNodalValues(std::vector<double> & flowdata,
                             const int partitionIdx,
                             const int timeStepIdx) = 0;

        //! Return number of partitions
        virtual int GetNumPartitions() { return numPartitions_;}

        //! Return number number of timestep files
        virtual int GetNumFiles() { return numFiles_;}

        //! Return maximum number of nodes
        virtual int GetNumNodes(const int partitionIdx) { 
            return MpCCInodes_[partitionIdx];
        }

        //! Return maximum number of nodes
        virtual int GetNumElems(const int partitionIdx) { 
            return MpCCIelems_[partitionIdx];
        }

        //! return dimension of grid
        virtual int GetDim() { return dim_;}

        //! return size of element
        virtual int GetElemSize(const int partitionIdx) {
            return elsize_[partitionIdx];
        }

        //! return size of element
        virtual int GetNumResults() { 
            return numResults_;
        }

        //! return size of element
        virtual double GetTimeStep();

    protected:


        //! ifstream which are associated with the fluid files
        std::ifstream infile;
        std::ifstream flowdatafile;
        std::string name_;
        std::string basename_;
        std::vector<int> elsize_;
        //<! current number of partitions
        int numPartitions_; 
        //<! dimension of grid
        int dim_;
        //<! maximum number of nodes
        std::vector<int> MpCCInodes_;
        //<! maximum number of elements
        std::vector<int> MpCCIelems_;
        int numResults_;
        int numFiles_;
    };

      
} // end of namespace
#endif
