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
    ScatteredDataReaderCSV(PtrParamNode& scatteredDataNode,
                           bool verbose = false);
    virtual ~ScatteredDataReaderCSV();
  
  
    void SetNumSkipLines(UInt skipLines) 
    {
      skipLines_ = skipLines;
    }

  protected:
    virtual void ReadData();

    void ParseParamNode();

  private:
    //! File name of input CSV file.
    std::string fileName_;

    //! Number of heading lines in .csv to skip.
    UInt skipLines_;

    std::map<UInt, UInt> dof2CoordColumn_;
    std::map< std::string, std::map<UInt, UInt> > qidDof2Column_;
  };
  
}

#endif  
