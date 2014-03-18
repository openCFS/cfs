#ifndef SCATTERED_DATA_READER_H
#define SCATTERED_DATA_READER_H

#include <string>
#include <vector>

#include "General/defs.hh"

namespace CoupledField 
{  
  class ScatteredDataReader;
  typedef boost::shared_ptr<ScatteredDataReader> ScatteredDataReaderPtr;

  /**
   * Base class for all scattered data readers
   */
  class ScatteredDataReader
  {
  public:
    ScatteredDataReader(const std::string& fileName, bool verbose = false) :
      fileName_(fileName),
      verbose_(verbose)
    {};
    virtual ~ScatteredDataReader() {};
  
    virtual void Read(std::vector< std::vector<double> >& scatteredData) = 0;
  protected:
    //! Name of scattered data file.
    std::string fileName_;

    //! Should verbose messages be printed?
    bool verbose_;
  };
  
}

#endif
