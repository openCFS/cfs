#ifndef SCATTERED_DATA_READER_CCM_H
#define SCATTERED_DATA_READER_CCM_H

#include <vector>
#include <map>
#include <set>

#include <ccmio.h>

#include "ScatteredDataReader.hh"

namespace CoupledField 
{
  /**
   * Class for reading cell-centered data from STAR-CCM+ .ccm files.
   *
   * This class reads the vertex coordinates and stores them in a map with
   * their global id as key. The vertex ids belonging to the polyhedral cells
   * are read from the internal faces data sets and cell centers are computed
   * from them. Data arrays are read and provided using their short names.
   * The member functions in this class are mostly based on the readexample.cpp
   * provided in the documentation of the CCMIO 2.6.1 library.
   */  
  class ScatteredDataReaderCCM : public ScatteredDataReader
  {
  public:
    ScatteredDataReaderCCM(const std::string& fileName,
                           bool verbose = false) :
      ScatteredDataReader(fileName, verbose) 
    {};
    virtual ~ScatteredDataReaderCCM();
  
  
    void ReadInput();
    void WriteCellCenters();

    void SetComponentShortNames(const std::vector<std::string>& componentShortNames) 
    {
      componentShortNames_ = componentShortNames;
    }

    virtual void Read(std::vector< std::vector<double> >& scatteredData);

  private:

    enum DataType { kScalar, kVector, kVertex, kCell, kInternalFace, kBoundaryFace,
                    kBoundaryData, kBoundaryFaceData, kCellType };

    void ReadMesh( CCMIOError &err, CCMIOID &vertices, CCMIOID &topology );
    void ReadVertices( CCMIOError &err, CCMIOID &vertices,
                       char const *counter, int offset = 0 );
    void ReadSolverInfo( CCMIOError &err, CCMIOID &solution);
    void ReadPost( CCMIOError &err, CCMIOID &solution );
    void ReadScalar( CCMIOError &err, CCMIOID& field, std::vector<int> &mapData,
                     std::vector<double> &data, bool readingVector = false,
                     const char* shortName = NULL);
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

    typedef std::map< int, std::set< int > > Entity2VertexMap;
    Entity2VertexMap cell2Verts_;
    Entity2VertexMap face2Verts_;

    typedef std::map< std::string, double > Name2FloatMap;
    typedef std::map< int, Name2FloatMap > Entity2DataMap;

    Entity2DataMap dMap_;
    Entity2DataMap f2dMap_;

    std::vector<std::string> componentShortNames_;
  };
  
}

#endif
