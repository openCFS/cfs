#ifndef FILE_FILEREADER_CFX_2006
#define FILE_FILEREADER_CFX_2006

#include <map>

#include "filereader.hh"

namespace CoupledField
{

    class FileReader_CFX : public FileReader
    {
    public:

        //! Constructor
        FileReader_CFX(const std::string& name,
                           const UInt dim,
                           const UInt numFiles);
    
        //! Deconstructor
        virtual ~FileReader_CFX();

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

        virtual double GetTimeStep();

    protected:

        void GetInfosFromCommand();

        void ParseCommand(std::vector<char>& cmdstr,
                          int& pos,
                          std::string& cmd,
                          std::string& attrib,
                          std::string indent,
                          std::ofstream& outFile);
        
        void ParseOption(std::vector<char>& cmdstr,
                         int& pos,
                         std::string& option,
                         std::string& value,
                         std::string indent,
                         std::ofstream& outFile);

        void IOErrorToString(int ioerr, std::string& errStr);

    private:

        double timeStep;
        std::string timeStepStr;
        std::string solTimeUnit;
        std::string defFile;
        std::map<std::string, std::string> exprMap;
        std::vector< std::string > transientFNs_;
        std::vector< int > timeStepNumbers_;

        char fn[80];
        int nerr;
        int whatfile;
        char what[80];
        char where[80];
        int when;
        int dattyp;
        int length;
        int nsize;
        int iopt;
        int ioptar;
        float rarr[80];
        int iarr[80];
        int nvx;
        char carr[80];
        int larr[80];
        double darr[80];
        char sarr[80];

        const static int BIGMEM;
        static std::vector<int> intvec;
        static std::vector<double> doublevec;
        static std::vector<float> floatvec;
        static std::vector<char> charvec;


    };

      
} // end of namespace
#endif
