// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <string>
#include <stdarg.h>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;

#include <def_cfs_stats.hh>
#include <def_cplreader.hh>

#include "General/exception.hh"
#include "ParamsInit.hh"
#include "Settings.hh"
#include "FileReader.hh"

#include "FileReaders/CFS++/FileReader_CfsHdf5.cc"

#ifdef CPLREADER_ANSYS
#include "FileReaders/ANSYS/FileReader_MKHDF5.hh"
#include "FileReaders/ANSYS/FileReader_NACS.hh"
#endif

#ifdef CPLREADER_FASTEST
#include "FileReaders/FASTEST/FileReader_FASTEST.hh"
#endif

#ifdef CPLREADER_FIELDVIEW
#include "FileReaders/FieldView/FileReader_fvuns.hh"
#endif

#ifdef CPLREADER_CFX
#include "FileReaders/CFX/FileReader_CFX.hh"
#include "FileReaders/CFXexport/FileReader_CFXexport.hh"
#endif

#ifdef CPLREADER_OPENFOAM
#include "FileReaders/VTKBased/OPENFOAM/FileReader_OPENFOAM.hh"
#endif

#ifdef CPLREADER_ENSIGHT
#include "FileReaders/VTKBased/EnSight/FileReader_EnSight.hh"
#endif

#ifdef CPLREADER_CGNS
#include "FileReaders/CGNS/FileReader_CGNS.hh"
#endif

#include "FileReaders/GENGRIDS/FileReader_GENGRIDS.hh"


#include "OutputWriters/hdf5/OutputWriter_HDF5.hh"

#include "Documentation/WriteDocumentation.hh"

#include "CouplingHandler.hh"


namespace CoupledField
{
  void FileReaderFactory(shared_ptr<FileReader>& fileReader) {
    Settings& settings = Settings::Instance();

    std::string type = settings.GetString("type");

    std::ostringstream sstr;
    if(type != "ANSYS" && type != "GENGRIDS" && type != "CFS++")
    {
      sstr << settings.GetString("basedir") << "/" << settings.GetString("name");
    }
    else
    {
      sstr << settings.GetString("basedir");
    }

    /* check if directory even exists */
    if(!fs::exists(sstr.str()) ||
        !fs::is_directory(sstr.str()))
      EXCEPTION("The directory '" << sstr.str()
          << "' doesn't exist!");

    // Create new file reader
    if(type == "ANSYS")
    {
#ifdef CPLREADER_ANSYS
      std::string readerType = FileReader_ANSYS::GetReaderType();

      if(readerType == "MKHDF5")
      {
        fileReader.reset(new FileReader_MKHDF5(settings.GetString("name"),
                                           settings.GetInt("dim"), 0, 0));
      }
      else if(readerType == "NACS")
      {
        fileReader.reset(new FileReader_NACS(settings.GetString("name"),
                                           settings.GetInt("dim"), 0, 0));
      }
      else
      {
        EXCEPTION("Could not find any suitable input files!");
      }
#else
      EXCEPTION("Reading of ANSYS files not supported!");
#endif
    }

    if(type == "FASTEST")
    {
#ifdef CPLREADER_FASTEST
      fileReader.reset(new FileReader_FASTEST(settings.GetString("name"),
                                          settings.GetInt("dim"),
                                          settings.GetInt("numsteps"),
                                          settings.GetInt("firststep")));
#else
      EXCEPTION("Reading of FASTEST files not supported!");
#endif
    }

    if(type == "FIELDVIEW")
    {
#ifdef CPLREADER_FIELDVIEW
      fileReader.reset(new FileReader_fvuns(settings.GetString("name"),
                                            settings.GetInt("dim"),
                                            settings.GetInt("numsteps"),
                                            settings.GetInt("firststep")));
#else
      EXCEPTION("Reading of FieldView files not supported!");
#endif
    }

    if(type == "GENGRIDS")
    {
      fileReader.reset(new FileReader_GENGRIDS(settings.GetString("name"),
                                               settings.GetInt("dim"),
                                               settings.GetInt("numsteps"),
                                               settings.GetInt("firststep")));
    }

 #if 0
    if(type == "Stanford")
    {
      fileReader.reset(new FileReader_Stanford(settings.GetString("name"),
                                           settings.GetInt("dim"),
                                           settings.GetInt("numsteps")));
    }
 #endif

    if(type == "CFX")
    {
#ifdef CPLREADER_CFX
      fileReader.reset(new FileReader_CFX(settings.GetString("name"),
                                      settings.GetInt("dim"),
                                      settings.GetInt("numsteps")));
#else
      EXCEPTION("Reading of CFX files not supported!");
#endif
    }

    if(type == "CFX_EXPORT")
    {
#ifdef CPLREADER_CFX
      fileReader.reset(new FileReader_CFXexport(settings.GetString("name"),
                                            settings.GetInt("dim"),
                                            settings.GetInt("numsteps"),
                                            settings.GetInt("firststep")));
#else
      EXCEPTION("Reading of CFX files not supported!");
#endif
    }

    if(type == "OPENFOAM")
    {
#ifdef CPLREADER_OPENFOAM
      fileReader.reset(new FileReader_OPENFOAM(settings.GetString("name"),
                                           settings.GetInt("dim"),
                                           settings.GetInt("numsteps")));
#else
      EXCEPTION("Reading of OPENFOAM files not supported!");
#endif
    }

    if(type == "ENSIGHT")
    {
#ifdef CPLREADER_ENSIGHT
      fileReader.reset(new FileReader_EnSight(settings.GetString("name"),
                                              settings.GetInt("dim"),
                                              settings.GetInt("numsteps")));
#else
      EXCEPTION("Reading of EnSight files not supported!");
#endif
    }

    if(type == "CGNS")
    {
#ifdef CPLREADER_CGNS
      fileReader.reset(new FileReader_CGNS(settings.GetString("name"),
                                           settings.GetInt("dim"),
                                           settings.GetInt("numsteps")));
      std::cerr << "test";
#else
      EXCEPTION("Reading of CGNS files not supported!");
#endif
    }

    if(type == "CFS++")
    {
      fileReader.reset(new FileReader_CfsHdf5(settings.GetString("name"),
                                              0,
                                              settings.GetInt("numsteps"),
                                              settings.GetInt("firststep")));
    }

    if(!fileReader)
    {
      EXCEPTION("ERROR: Could not initialize " << type
                << " filereader.");
    }
  }

  void OutputWriterFactory(OutputWriterVectorType& outputWriters) {
    shared_ptr<OutputWriter> p(new OutputWriter_HDF5());
    outputWriters.push_back(p);
  }
}


using namespace CoupledField;

int main(int argc, char *argv[])
{
  int ret = 0;
  shared_ptr<FileReader> fileReader;
  OutputWriterVectorType outputWriters;

  try
  {
    Settings& settings = Settings::Instance();

  std::cout << std::endl
            << "============================================================"
            << "===========" << std::endl;
  std::cout << " cplreader (CFS++) - A fluid data reader for CFS++/MpCCI coupling" << std::endl << std::endl
            << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'"
            << " (rev " << CFS_SUBVERSION_REV << ")" << std::endl
            << " compiled " << __DATE__
            << " as " << CMAKE_BUILD_TYPE << std::endl;
  std::cout << "============================================================"
            << "==========="
            << std::endl << std::endl;

    // Initialize settings singleton with command line parameters
    ParamsInit(argc, argv);

    // Set user defined exception behaviour
    Exception::segfault_ = settings.GetInt("segfault") != 0;

    // Set Enums from environment.hh
    SetEnvironmentEnums();

    // Write HTML help
    if(settings.GetString("docu") != "")
    {
      WriteDocumentation();
      return 0;
    }

    // Generate a file reader
    FileReaderFactory(fileReader);

    // Generate output writers
    OutputWriterFactory(outputWriters);

    // Initialize and perform coupling
    CouplingHandler cplHandler(fileReader, outputWriters);
    cplHandler.Init(argc, argv);
    if(!settings.GetInt("justinit"))
    {
      cplHandler.ConvertMesh();

      if(!settings.GetInt("justmesh"))
      {
        cplHandler.Couple();
      }
    }
    cplHandler.Finish();

  } catch (std::exception& ex)
  {
    std::cerr << std::endl << std::endl
              << "CAUGHT EXCEPTION:" << std::endl
              << std::endl
              << ex.what()
              << std::endl;
    ret = 1;
  }

  return ret;
}

