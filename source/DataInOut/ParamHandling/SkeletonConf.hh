// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SKELETONCONF
#define FILE_SKELETONCONF

/**************************************************************************/
/* File:   SkeletonCONF.hh                                                */
/* Author: Manfred Kaltenbacher                                           */
/* Date:   14. Nov. 2003                                                  */
/*                                                                        */
/* Writes a skeleton of the config file.                                  */
/**************************************************************************/

#include "DataInOut/simInput.hh"

namespace CoupledField
{

  //! Class for writing  a skeleton of the config file
  class SkeletonConf {

  public:

    //! Constructor
    SkeletonConf(SimInput* simInput);

    //! destructor
    virtual ~SkeletonConf();

    //! writes the conf-file
    void WriteConf();

    //! writes general information 
    void WriteGeneral();

    //! writes number and list of subdomains 
    void WriteSubdomains();

    //! writes list_surfaces, list_edges and list_nodes
    void WriteLists();

    //! writes skeleton for PDE specific parameters
    void WritePDE();

    //! writes skeleton for PDE coupling list
    void WriteCouplingList();

    //! write analysis types
    void WriteAnalysisTypes();

  private:

    //! Helper methods

    std::string Indent( Integer indentChange = 0) {
      level_ = level_ + indentChange;
      std::string help;
      if( level_ < 0 ) 
        EXCEPTION( "Level must not be < 0 ");
      
      for(Integer i=0; i<level_; i++) help+= "  ";
      return help;
    }

    std::string Quote( const std::string& in, Integer indentChange = 0) {
      std::string ret = Indent(indentChange);
      ret += "<!-- " + in + " -->";
      return ret;
    }

    //! Current indentation level
    Integer level_;

    //! file pointer fo xml-file
    std::ofstream * out_; 

    //! file pointer of mesh-file
    SimInput* simInput_;   

  };
} // end namespace CoupledField
 
#endif
