#ifndef DEFINEFILES_2001
#define DEFINEFILES_2001

namespace CoupledField
{ 

  //! Define trace, debug, info files
  /*!
    In this class we define auxiliary files, such as {\it trace }-file (file,where we list all methods and classes, which were used during running of program; this
    file should help developer to trace a mistake), {\it debug}-file (file, where intermediate results are stored), {\it info} - file (in this file we print specific information about methods and types of data, which were used in code)
  */

class DefineInOutFiles
{
 public:

   //!
   DefineInOutFiles(const Char * NameOfInputFile);

   //!
   ~DefineInOutFiles();

   //!
   FileType * Create_ptFileType(Char * atype);
 
 private:

   //!
   Char * filename;

   //!
   FileType * infiletype;
};

} // end of namespace
#endif

