#ifndef FILE_SKELETONCONF
#define FILE_SKELETONCONF

/**************************************************************************/
/* File:   SkeletonCONF.hh                                                */
/* Author: Manfred Kaltenbacher                                           */
/* Date:   14. Nov. 2003                                                  */
/*                                                                        */
/* Writes a skeleton of the config file.                                  */
/**************************************************************************/

#include "DataInOut/filetype.hh"

namespace CoupledField
{

  //! Class for writing  a skeleton of the config file
  
  class SkeletonConf
  {
  public:
    //! constructor
    /*!
      \param name name of input file without extension
    */
    SkeletonConf(const Char * name);
    
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

  private:

    std::ofstream * skelfile_; //!< file pointer fo conf-file
    FileType * meshfile_;      //!< file pointer of mesh-file
    Char * name_;              //!< stores name of input file (without extension)
  };
} // end namespace CoupledField
 
#endif
