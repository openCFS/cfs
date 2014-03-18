#ifndef SCATTERED_DATA_READER_CSV_H
#define SCATTERED_DATA_READER_CSV_H

#include "ScatteredDataReader.hh"

namespace CoupledField 
{
  /**
   * Class for reading scattered data from comma separated files.
   */  
  class ScatteredDataReaderCSV : public ScatteredDataReader
  {
  public:
    ScatteredDataReaderCSV(const std::string& fileName,
                           bool verbose = false) :
      ScatteredDataReader(fileName, verbose)
    {};
    virtual ~ScatteredDataReaderCSV();
  
  
    void SetNumSkipLines(UInt skipLines) 
    {
      skipLines_ = skipLines;
    }

    virtual void Read(std::vector< std::vector<double> >& scatteredData);

  private:

    //! Number of heading lines in .csv to skip.
    UInt skipLines_;
  };
  
}

#endif  
