// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef DEFINEFILES_2001
#define DEFINEFILES_2001

#include <boost/shared_ptr.hpp>

#include "DataInOut/SimInput.hh"
#include "DataInOut/SimOutput.hh"
#include "DataInOut/ParamHandling/MaterialHandler.hh"

namespace CoupledField
{ 

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

    //! create input readers 
    void CreateSimInputFiles( PtrParamNode paramNode, 
                              PtrParamNode infoNode,
                              std::map<std::string, shared_ptr<SimInput> >& inFiles,
                              std::map<std::string, 
                              StdVector<shared_ptr<SimInput> > >& gridInputs );

    //! Create pointer to output classes and read the corresponding gridIds
    void CreateSimOutputFiles( PtrParamNode paramNode,
                               PtrParamNode infoNode,
                               std::map<std::string,  
                               shared_ptr<SimOutput> >&  out,
                               std::map<std::string, std::string>& gridIds );

    //! create pointer to Materialfile Handler
    shared_ptr<MaterialHandler> CreateMaterialHandler( PtrParamNode rootNode );

    //! Generic function to obtain a single SimInput object based on the given parameters
    static shared_ptr<SimInput>  CreateSingleInputFileObject(std::string fName,
                                                             std::string simName,
                                                             PtrParamNode configNode,
                                                             PtrParamNode infoNode);

    //! Generic function to obtain a single SimOutpu object based on the given parameters
    static shared_ptr<SimOutput> CreateSingleOutputFileObject(std::string fName,
                                                              PtrParamNode configNode,
                                                              PtrParamNode infoNode,
                                                              bool isRestart);

  private:

    //! pointer to MaterialHandler
    shared_ptr<MaterialHandler> ptMaterialHandler_;
    
  };

} // end of namespace
#endif

