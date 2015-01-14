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
    ScatteredDataReaderCCM(PtrParamNode& scatteredDataNode,
                           bool verbose = false);
    virtual ~ScatteredDataReaderCCM();
  
  
    void ReadInput();
    void Dump();

  protected:
    virtual void ReadData();

    void ParseParamNode();

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
    //! File name of input CCM file.
    std::string fileName_;

    //! Id of reader.
    std::string id_;

    //! Shall we dump the point cloud data to a CSV file?
    bool dump_;

    //! Do we need to read data from cell centers?
    bool cellCenters_;

    //! Do we need to read data from face centers?
    bool faceCenters_;

    //! Number of values of each element to print when in verbose mode
    static int const kNValues;
    static char const kUnitsName[];
    static int const kVertOffset;
    static int const kCellInc;

    //! Data type for associating a node number with coordinates
    typedef std::map< int, std::vector< double > > VertexMap;
    //! Vertices of face and cell nodes
    VertexMap vertices_;

    //! Data type for uniquely associating vertex ids to cell or face ids
    typedef std::map< int, std::set< int > > Entity2VertexMap;
    //! Vertices belonging to cells
    Entity2VertexMap cellVertices_;
    //! Vertices belonging to faces
    Entity2VertexMap faceVertices_;

    typedef std::map< std::string, double > Name2FloatMap;
    typedef std::map< int, Name2FloatMap > Entity2DataMap;

    //! Map from global cell number to data defined on it.
    Entity2DataMap cell2Data_;

    //! Map from global face number to data defined on it.
    Entity2DataMap face2Data_;

    //! Map from quantity id to map of quantity dof indices to short names.
    std::map< std::string, std::map<UInt, std::string> > qidDof2ShortName_;

    //! Set of short names of quantity dofs which are to be read.
    std::set< std::string > componentShortNames_;
  };  
}

#endif
