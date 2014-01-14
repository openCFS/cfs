#ifndef SCATTERED_DATA_READER_CSV_H
#define SCATTERED_DATA_READER_CSV_H

#include "ScatteredDataReader.hh"

namespace CoupledField 
{
  
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

    UInt skipLines_;
  };
  
}

#endif  
