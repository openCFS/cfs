#ifndef SCATTERED_DATA_READER_CSVT_H
#define SCATTERED_DATA_READER_CSVT_H

#include <boost/signals2.hpp>
#include "ScatteredDataReader.hh"

namespace CoupledField 
{
  class MathParser;
  /**
   * Class for reading scattered data from comma separated files.
   */  
  class ScatteredDataReaderCSVT : public ScatteredDataReader
  {
  public:
    ScatteredDataReaderCSVT(PtrParamNode& scatteredDataNode,
                           bool verbose = false);
    virtual ~ScatteredDataReaderCSVT();
  
  
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
    std::string fileNameGeom_; //ADDED

    //! Number of heading lines in .csv to skip.
    UInt skipLines_;

    std::map<UInt, UInt> dof2CoordColumn_;
    std::map< std::string, std::map<UInt, UInt> > qidDof2Column_;

    UInt time2Col_;
    UInt file2Col_;

    StdVector< Double > valuetime_;
    StdVector< std::string > valueFileName_;

    // current time step
    Integer step_;

    //! Pointer to math parser instance
    MathParser* mp_ = nullptr;

    //! Handle for expression determines current time/freq value
    unsigned int mHandleStep_ = 0;

    //! Handle callback time/freq value changed
    boost::signals2::connection connReal_;
  };
  
}

#endif  
