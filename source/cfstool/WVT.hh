// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef WVT_HEADER_2013
#define WVT_HEADER_2013

#include <string>

#include "General/Environment.hh"
#include "DataInOut/SimInput.hh"

namespace CFSTool
{

  /** This is the WVT class for weight vector theory-based post-processing.
   *
   *  This is the WVT class for the evaluation of post-processing results
   *  specific to coriolis mass flow meters (CMFs) according to the
   *  weight vector theory developed in \cite Hemp1994 and \cite Hemp1995. */
  class WVT 
  {
  public:

    WVT( const PtrParamNode& param,
         const PtrParamNode& info );

    ~WVT() {};
    
    typedef enum {PRIMARY_MODE, SECONDARY_MODE, MEAN_FLOW} InputFileType;

    static Enum<InputFileType> inputFileType;
    
  private:

    const PtrParamNode& param_;
    const PtrParamNode& info_;

    //! File name of input file for primary mode. Usually, this is the lateral
    //! mode.
    std::string priModeFile_;

    //! File name of input file for secondary mode. Usually, this is the coriolis
    //! mode.
    std::string secModeFile_;

    //! File name of input file for mean flow.
    std::string meanFlowFile_;

    //! File name of output file, if output for post-processing weight vectors
    //! should be written.
    std::string outFile_;
        
    bool writeOutputFile_;
    
    typedef std::map< InputFileType, shared_ptr<SimInput> > InputsType;
    InputsType inputs_;
  };

}

#endif
