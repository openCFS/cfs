#ifndef DEFINEFILES_2001
#define DEFINEFILES_2001

#include "filetype.hh"
#include "writeresults.hh"

namespace CoupledField
{ 

  //! Define trace, debug, info files
  /*!
    In this class we define auxiliary files, such as {\it trace }-file (file,where we list all methods and classes, that are used during running of program;
 this
    file should help developer to trace a mistake), {\it debug}-file (file, where intermediate results are stored), {\it info} - file (in this file we print specific information about methods and types of data, which were used in code)
  */
class DefineInOutFiles
{
public:

   //! constructor
  /*!
    \param name name of trace file without extension
  */
  DefineInOutFiles(const Char * name);

   //! deconstructure
   ~DefineInOutFiles();

   //! create a pointer to a class for reading input-reults, a derived class of the FileType according to the specification of the conf-file.
   FileType * Create_ptFileType();

   //! create a pointer to a class for writing output-results, a derived class of the WriteResults according to the specification of the conf-file.
   WriteResults * Create_ptWriteResults();

private:

   //! name of all auxialiary files
   Char * filename;

   //! pointer to a class for reading input mesh-data
   FileType * infiletype;

   //! pointer to a class for writing output data
   WriteResults * ptWriteResults;   
};

} // end of namespace
#endif

