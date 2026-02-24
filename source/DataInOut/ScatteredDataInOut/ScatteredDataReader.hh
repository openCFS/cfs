#ifndef SCATTERED_DATA_READER_H
#define SCATTERED_DATA_READER_H

#include <string>
#include <vector>
#include <map>
#include <set>

#include "DataInOut/ParamHandling/ParamNode.hh"

#include "General/defs.hh"

namespace CoupledField 
{  
  class ScatteredDataReader;
  typedef std::shared_ptr<ScatteredDataReader> ScatteredDataReaderPtr;

  /**
   * Base class for all scattered data readers
   */
  class ScatteredDataReader
  {
  public:

      enum TransientMode{
        TF, STAT
      };

    ScatteredDataReader(PtrParamNode& scatteredDataNode, bool verbose = false) :
      myParamNode_(scatteredDataNode),
      verbose_(verbose),
      readerMode_(STAT)
    {};
    virtual ~ScatteredDataReader() {};


    TransientMode GetMode(){
      return readerMode_;
    }
  
    //! The CreateReaders function generates the different readers.
    static void CreateReaders(PtrParamNode& scatteredDataNode); 

    //! CoefFunctionScatteredData   instances  can   register  their   desired
    //! quantities using RegisterQuantity.
    static void RegisterQuantity(const std::string& quantity);

    //! The static Read function dispatches calls to the different readers.
    static void Read(bool updateMode=false);

    //! The  GetQuantity   function  retrieves  the  desired   data  from  the
    //! corresponding reader.
    static void GetQuantity(const std::string& quantity,
                            std::vector< std::vector<double> >& coordinates,
                            std::vector< std::vector<double> >& scatteredData);

    static void GetQuantity(const std::string& quantity,
                            std::vector< std::vector<double> >& coordinates,
                            std::vector< std::vector<Complex> >& scatteredData);


  protected:

    //! Pure  virtual  function  for  actually  reading  scattered  data  from
    //! files. Must be implemented by each derived reader class.
    virtual void ReadData() = 0;

    //! ParamNode containing information for this reader instance.
    PtrParamNode myParamNode_;

    //! Should verbose messages be printed?
    bool verbose_;
    
    //! A map from reader ids to the actual readers.
    static std::map<std::string, std::shared_ptr<ScatteredDataReader> > readers_;

    //! A map from quantity ids to reader ids.
    static std::map<std::string, std::string > quantities2Readers_;

    //! What quantities does this reader handle?
    std::set<std::string> registeredQuantities_;

    //! Array for point coordinates.
    std::vector< std::vector<double> > coordinates_;

    //! Map from  quantity id to a  vector containing the actual  data for the
    //! quantity.
    std::map< std::string, std::vector< std::vector<double> > > scatteredDataPerQuantity_;
    std::map< std::string, std::vector< std::vector<double> > > scatteredDataPerQuantityImag_;

    //! store mode of reader
    TransientMode readerMode_;


  };
  
}

#endif
