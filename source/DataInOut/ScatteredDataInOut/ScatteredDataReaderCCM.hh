#ifndef SCATTERED_DATA_READER_CCM_H
#define SCATTERED_DATA_READER_CCM_H

#include <vector>
#include <map>
#include <set>

#include <ccmio.h>

#include "ScatteredDataReader.hh"

namespace CoupledField 
{
  
  class ScatteredDataReaderCCM : public ScatteredDataReader
  {
  public:
    ScatteredDataReaderCCM(const std::string& fileName,
                           bool verbose = false) :
      ScatteredDataReader(fileName, verbose) 
    {};
    virtual ~ScatteredDataReaderCCM();
  
  
    enum DataType { kScalar, kVector, kVertex, kCell, kInternalFace, kBoundaryFace,
                    kBoundaryData, kBoundaryFaceData, kCellType };
  
    void ReadInput();
    void WriteCellCenters();

    void SetComponentShortNames(const std::vector<std::string>& componentShortNames) 
    {
      componentShortNames_ = componentShortNames;
    }

    virtual void Read(std::vector< std::vector<double> >& scatteredData);

  private:

    void ReadMesh( CCMIOError &err, CCMIOID &vertices, CCMIOID &topology );
    void ReadVertices( CCMIOError &err, CCMIOID &vertices,
                       char const *counter, int offset = 0 );
    void ReadSolverInfo( CCMIOError &err, CCMIOID &solution);
    void ReadPost( CCMIOError &err, CCMIOID &solution );
    void ReadScalar( CCMIOError &err, CCMIOID field, std::vector<int> &mapData,
                     std::vector<float> &data, bool readingVector = false );
    void ReadSets( CCMIOError &err, CCMIOID &problem );
    void CheckError( CCMIOError const &err, char const *str );
    void PrintData( int n, int id, DataType type, void *data,
                    void *data2 = NULL );
  
  private:
    static int const kNValues; // Number of values of each element to print
    static char const kUnitsName[];
    static int const kVertOffset;
    static int const kCellInc;

    typedef std::map< int, std::vector< double > > VertexMap;
    VertexMap vertices_;

    typedef std::map< int, std::set< int > > Cell2VertexMap;
    Cell2VertexMap cell2Verts_;

    typedef std::map< std::string, float > Name2FloatMap;
    typedef std::map< int, Name2FloatMap > Cell2DataMap;

    Cell2DataMap dMap_;

    std::vector<std::string> componentShortNames_;
  };
  
}

#endif  
