#ifndef DEFINEFILES_2001
#define DEFINEFILES_2001

#include "DataInOut/filetype.hh"
#include "DataInOut/writeresults.hh"

namespace CoupledField
{ 

  //! Define trace, debug, info files

  //! In this class we define auxiliary files, such as the {\it trace }-file
  //! (file,where we list all methods and classes, that are used during
  //! running of program; this file should help developer to trace a mistake),
  //! {\it debug}-file (file, where intermediate results are stored),
  //! {\it info} - file (in this file we print specific information about
  //! methods and types of data, which were used in code)
  class DefineInOutFiles
  {

  public:

    //! constructor
    //! \param name This is the basename for all auxilliary file.
    DefineInOutFiles(const Char * name);

    //! deconstructor
    ~DefineInOutFiles();

    //! Open an auxilliary file
    void OpenFile( AuxFileType fileType );

    //! create a pointer to a class for reading input-results, a derived class
    //! of the FileType according to the specification of the conf-file.
    FileType* CreateMeshFileHandler();

    //! create a pointer to a class for writing output-results, a derived class
    //! of the WriteResults according to the specification of the conf-file.
    WriteResults* Create_ptWriteResults( FileType * const aInFile );

  private:

    //! pointer to a class for reading input mesh-data
    FileType *infileType_;

    //! pointer to a class for writing output data
    WriteResults *ptWriteResults_;

  };

} // end of namespace
#endif

