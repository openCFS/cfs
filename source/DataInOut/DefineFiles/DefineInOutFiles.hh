// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef DEFINEFILES_2001
#define DEFINEFILES_2001

#include "DataInOut/SimInput.hh"
#include "DataInOut/SimOutput.hh"
#include "DataInOut/MaterialHandler.hh"

namespace CoupledField
{ 

    class MeshIOModule;

  //! Define trace, debug, info files

  //! In this class we define auxiliary files, such as the <em>trace</em>-file
  //! (file,where we list all methods and classes, that are used during
  //! running of program; this file should help developer to trace a mistake),
  //! <em>debug</em>-file (file, where intermediate results are stored),
  //! <em>info</em>- file (in this file we print specific information about
  //! methods and types of data, which were used in code)
  class DefineInOutFiles
  {

  public:

    //! constructor
    DefineInOutFiles();

    //! destructor
    ~DefineInOutFiles();

    //! Open an auxilliary file
    void OpenFile( AuxFileType fileType );

    //! create input readers 
    void CreateSimInputFiles( std::map<std::string, shared_ptr<SimInput> >& inFiles,
                              std::map<std::string, 
                              StdVector<shared_ptr<SimInput> > >& gridInputs );

    //! create pointer to output class and their related ids
    void CreateSimOutputFiles( std::map<std::string, 
                               shared_ptr<SimOutput> >&  out );

    //! create pointer to Materialfile Handler
    MaterialHandler* CreateMaterialHandler();

  private:

    //! pointer to a class for reading input data
    SimInput* simInput_;

    //! pointer to MaterialHandler
    MaterialHandler * ptMaterialHandler_;

  };

} // end of namespace
#endif

