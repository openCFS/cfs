#ifndef FILE_FILEREADER_2006
#define FILE_FILEREADER_2006

#include <fstream>
#include <vector>

#include "cplreaderdefs.hh"

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
        FileReader(const std::string& name, const UInt dim, const UInt numFiles);
    
        //! Deconstructor
        virtual ~FileReader();

        virtual void Init() = 0;

        //! get node coordinates from the corresponding file
        virtual void ReadNodalCoords(std::vector<Double> & NODECOORD,
                             const UInt partitionIdx) = 0;


        //! get topology information from the corresponding topology file
        virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                  std::vector<UInt> & numNodesPerElem,
                                  std::vector<UInt> & elemTypes,
                                  const UInt partitionIdx) = 0;

        //! get nodal values from the corresponding fluid datafile
        virtual void ReadNodalValues(std::vector<Double> & flowdata,
                             const UInt partitionIdx,
                             const UInt timeStepIdx) = 0;

        //! Return number of partitions
        virtual UInt GetNumPartitions() { return numPartitions_;}

        //! Return number number of timestep files
        virtual UInt GetNumFiles() { return numFiles_;}

        //! Return maximum number of nodes
        virtual UInt GetNumNodes(const UInt partitionIdx) { 
            return MpCCInodes_[partitionIdx];
        }

        //! Return maximum number of nodes
        virtual UInt GetNumElems(const UInt partitionIdx) { 
            return MpCCIelems_[partitionIdx];
        }

        //! return dimension of grid
        virtual UInt GetDim() { return dim_;}

        //! return size of element
        virtual UInt GetElemSize(const UInt partitionIdx) {
            return elsize_[partitionIdx];
        }

        //! return number of result quantities
        virtual UInt GetNumResults() { 
            return numResults_;
        }

        //! return time step
        virtual double GetTimeStep();

    protected:


        //! ifstream which are associated with the fluid files
        std::ifstream infile;
        std::ifstream flowdatafile;
        std::string name_;
        std::string basename_;
        std::vector<UInt> elsize_;
        //<! current number of partitions
        UInt numPartitions_; 
        //<! dimension of grid
        UInt dim_;
        //<! maximum number of nodes
        std::vector<UInt> MpCCInodes_;
        //<! maximum number of elements
        std::vector<UInt> MpCCIelems_;
        UInt numResults_;
        UInt numFiles_;
    };

      
} // end of namespace
#endif
