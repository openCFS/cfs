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

#include <DataInOut/simInput.hh>

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

    std::ofstream * skelfile_; //!< file pointer fo conf-file
    SimInput* simInput_;      //!< file pointer of mesh-file
  };
} // end namespace CoupledField
 
#endif
