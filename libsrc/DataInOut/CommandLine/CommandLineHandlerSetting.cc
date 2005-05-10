#include <string>

#include "CommandLineHandlerSetting.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  CommandLineHandlerSetting::CommandLineHandlerSetting( Integer argc,
                                                        const char **argv ) {

    ENTER_FCN( "CommandLineHandlerSetting::CommandLineHandlerSetting" );


    // ----------------------------------------
    //  Name of simulation run goes separately
    // ----------------------------------------

    // We expect the last command line argument to be the name of
    // the current simulation run. Check that a proper name was
    // provided by the user.
    if ( argc > 1 ) {
      simName_ = argv[argc-1];
      std::cout << " --> Using '" << simName_ << "' as basename for the "
                << "simulation\n\n";
    }
    else {
      std::cerr << " Please specify a name for the simulation run\n\n";
      exit(-1);
    }


    // -----------------------------------------------------------
    //  Generate specification of allowed command line parameters
    // -----------------------------------------------------------

    const SettingDef commandLineParamsSpec[] = {

      // --paramFile
      SettingDef( markerParamFile_.c_str(),
                  markerLongParamFile_.c_str(),
                  Setting::STRING,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY | Setting::EXPLICIT_ASSIGNMENT,
                  helpParamFile_.c_str() ),

      // --meshFile
      SettingDef( markerMeshFile_.c_str(),
                  markerLongMeshFile_.c_str(),
                  Setting::STRING,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY | Setting::EXPLICIT_ASSIGNMENT,
                  helpMeshFile_.c_str() ),

      // --traceDepth
      SettingDef( markerTraceDepth_.c_str(),
                  markerLongTraceDepth_.c_str(),
                  Setting::INT,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY | Setting::EXPLICIT_ASSIGNMENT,
                  helpTraceDepth_.c_str() ),

      // --writeSkeleton
      SettingDef( markerWriteSkeleton_.c_str(),
                  markerLongWriteSkeleton_.c_str(),
                  Setting::FLAG,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY,
                  helpWriteSkeleton_.c_str() ),

      // -printGrid
      SettingDef( markerPrintGrid_.c_str(),
                  markerLongPrintGrid_.c_str(),
                  Setting::FLAG,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY,
                  helpPrintGrid_.c_str() ),

      // empty final one
      SettingDef()

    };


    // --------------------
    //  Parse command line
    // --------------------

    // Let the SettingDataBase object parse the command line
    // The simulation name does not belong to the declared parameters
    // for SettingDataBase, so hide it from the object by decrementing
    // the command line parameter counter
    bool pickyResult = true;
    commandLine_.readSettings( commandLineParamsSpec, argv, (int)argc-1,
                               false, false, &pickyResult );

    if ( pickyResult == true ) {
      std::cout << "\n\n Successfully parsed command line\n\n";
    }

    commandLine_.finalCheck( commandLineParamsSpec );

  };

}
