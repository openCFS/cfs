#include <iostream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>

namespace fs = boost::filesystem;

#include "General/exception.hh"
#include "ParamsInit.hh"
#include "Settings.hh"

namespace CoupledField
{

  namespace po = boost::program_options;

  std::auto_ptr<Settings> Settings::settingsInstance_;

  void ParamsInit(int argc, char* argv[])
  {
    Settings& settings = Settings::Instance();
    std::string param_name;
    std::string param_coupling;
    uint32_t param_dim;
    uint32_t param_numsteps;
    uint32_t param_firststep;
    uint32_t param_stepincr;
    std::string param_type;
    bool param_calcsrc;
    bool param_verbose;
    bool param_floatds;
    bool param_justinit;
    bool param_justmesh;
    bool param_degen;
    bool param_strict;
    bool param_deltemp;
    bool param_segfault;
    std::string param_deffile;
    double param_timestep;
    std::string param_dump;
    std::string param_quantities;
    std::string param_lhsrc;
    std::string param_vx;
    std::string param_vy;
    std::string param_vz;
    std::string param_pres;
    std::string param_outputfields;
    std::string param_outprec;
    std::string param_activeparts;
    std::string param_basedir;
    bool param_extfiles;


    if(argc == 2) {
      std::string arg = argv[1];
      if(arg == "-?" || arg == "--help")
      {
        std::cout << std::endl
                  << "cplreader (CFS++) - A fluid data reader for CFS++/MpCCI coupling"
                  << std::endl << std::endl
                  << "Compiled:" << std::endl << "  "
                  << __DATE__ << std::endl << std::endl;
      }
    }

    int myargc = argc;
    char** myargv = argv;
    char* dummy[32];

    for (int n=1; n<argc; n++)
    {
      std::string arg = argv[n];

      if(arg == "-p4amslave" ||
         arg == "-p4yourname" ||
         arg == "-p4rmrank")
      {
        myargc = argc - 8;

        dummy[0] = argv[0];
        for(int i=0; i<myargc; i++)
        {
          dummy[i+1] = argv[8+i];
        }

        myargv = dummy;
        myargc++;

        break;

      }
    }

    try
    {
      // Declare the supported options.
      po::options_description desc("Allowed options");
      desc.add_options()
        ("help", "Produce help message")
        ("name", po::value< std::string >(&param_name)->default_value("dataset"),
         "Name of dataset as well as of the directory containing the dataset")

        ("type", po::value< std::string >(&param_type)->default_value("CFX"),
         "Type of dataset (can be CFX | CFX_EXPORT | FASTEST | OPENFOAM)")

        ("coupling", po::value< std::string >(&param_coupling)->default_value("file"),
         "Specify kind of coupling MpCCI | file")

        ("dim", po::value< uint32_t >(&param_dim)->default_value(3),
         "Dimension of fluid data. (can be 2|3)")

        ("numsteps", po::value< uint32_t >(&param_numsteps)->default_value(0),
         "Number N of timesteps to read")

        ("firststep", po::value< uint32_t >(&param_firststep)->default_value(1),
         "Index of first time step (only for CFX_EXPORT, FASTEST)")

        ("stepincr", po::value< uint32_t >(&param_stepincr)->default_value(1),
         "Step increment for reading the files")

        ("timestep", po::value< double >(&param_timestep)->default_value(-1),
         "Time step length T in seconds")

        ("deffile", po::value< std::string >(&param_deffile)->default_value(""),
         "Definition file name. Only for CFX!")

        ("calcsrc", po::value< bool >(&param_calcsrc)->default_value(false),
         "Calculate the acoustic sources from velocity")

#ifdef DEBUG
        ("segfault", po::value< bool >(&param_segfault)->default_value(true),
         "Cause program to segfault when an exception is encountered.")
#endif
        ("deltemp", po::value< bool >(&param_deltemp)->default_value(true),
         "Delete temporary files.")

        ("floatds", po::value< bool >(&param_floatds)->default_value(true),
         "Do the transient files contain float or double values."
         "ATTENTION: If this flag is true, it overrides the setting of"
         "--outprec for CFX!")

        ("outprec", po::value< std::string >(&param_outprec)->default_value("double"),
         "Write data as single or as double precision (can be single|double)")

        ("justinit", po::value< bool >(&param_justinit)->default_value(false),
         "Just perform initialization step.")

        ("justmesh", po::value< bool >(&param_justmesh)->default_value(false),
         "Just convert the mesh.")

        ("degen", po::value< bool >(&param_degen)->default_value(false),
         "Allow degenerated element shapes.")

        ("strict", po::value< bool >(&param_strict)->default_value(false),
         "Be strict about holes in node list and connectivity.")

        ("basedir", po::value< std::string >(&param_basedir)->default_value("."),
         "Base directory where subdirectories are searched.")

        ("extfiles", po::value< bool >(&param_extfiles)->default_value(true),
         "Use external time step files (just like CFX .trn files).")
/*
        ("dump", po::value< std::string >(&param_dump)->default_value("GID"),
         "Dump the grid and data to a file (can be GID|GMV)")
*/
        ("verbose", po::value< bool >(&param_verbose)->default_value(false),
         "Be verbose")

        ("lhsrc", po::value< std::string >(&param_lhsrc)->default_value(""),
         "Column of Lighthill source term in FASTEST result files (e.g. col1).")

        ("vx", po::value< std::string >(&param_vx)->default_value(""),
         "Column of x-velocity in FASTEST result files (e.g. col2).")

        ("vy", po::value< std::string >(&param_vy)->default_value(""),
         "Column of y-velocity in FASTEST result files (e.g. col3).")

        ("vz", po::value< std::string >(&param_vz)->default_value(""),
         "Column of z-velocity in FASTEST result files (e.g. col4).")

        ("pres", po::value< std::string >(&param_pres)->default_value(""),
         "Column of pressure in FASTEST result files (e.g. col5).")

        ("activeparts", po::value< std::string >(&param_activeparts)->default_value("all"),
         "Values will only be output on partitions specified by activeparts"
         "(all | numbers seperated by SPACE or ; or |). ")

        ("outputfields", po::value< std::string >(&param_outputfields)->default_value("acouRhsLoad"),
         "Which physical fields should be output "
         "([acouRhsLoad | fluidMechDensity | fluidMechPressure | fluidMechVelocity | fluidMechTKE]). "
         "Values may be separated by SPACE or ; or |")


        ;

      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);

      if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(0);
      }
    }
    catch(std::exception &ex)
    {
      EXCEPTION("Error while parsing command line: " << ex.what());
    }

    std::string baseDir;
    std::string baseName;

    try
    {
      fs::path fn = fs::system_complete(argv[0]);
      fn.normalize();
      baseDir = fn.branch_path().native_directory_string();
      baseName = (fs::change_extension( fn.leaf(), "" )).native_directory_string();
    } catch (fs::filesystem_error& ex)
    {
      EXCEPTION("Received exception: " << ex.what());
      return;
    }

#ifndef DEBUG
    param_segfault = false;
#endif

    settings.SetString("name", param_name);
    settings.SetString("type", param_type);
    settings.SetString("coupling", param_coupling);
    settings.SetInt("dim", param_dim);
    settings.SetInt("numSteps", param_numsteps);
    settings.SetInt("firstStep", param_firststep);
    settings.SetInt("stepIncr", param_stepincr);
    settings.SetInt("calcSrc", param_calcsrc);
    settings.SetInt("floatDataset", param_floatds);
    settings.SetString("outprec", param_outprec);
    settings.SetInt("justinit", param_justinit);
    settings.SetInt("justmesh", param_justmesh);
    settings.SetInt("degen", param_degen);
    settings.SetInt("strict", param_strict);
    settings.SetInt("segfault", param_segfault);
    settings.SetInt("deltemp", param_deltemp);
    settings.SetInt("verbose", param_verbose);
    settings.SetString("defFile", param_deffile);
    settings.SetDouble("timeStep", param_timestep);
    settings.SetString("dump", param_dump);
    settings.SetString("exename", baseName);
    settings.SetString("baseExeDir", baseDir);
    settings.SetString("basedir", param_basedir);
    settings.SetString("lhsrc", param_lhsrc);
    settings.SetString("vx", param_vx);
    settings.SetString("vy", param_vy);
    settings.SetString("vz", param_vz);
    settings.SetString("pres", param_pres);
    settings.SetString("activeparts", param_activeparts);
    settings.SetString("outputfields", param_outputfields);
    settings.SetInt("extfiles", param_extfiles);
  }


}
