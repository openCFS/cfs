#include <string>

#include "CommandLineHandlerSetting.hh"


namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  CommandLineHandlerSetting::CommandLineHandlerSetting( Integer argc,
                                                        const char **argv ) {

    ENTER_FCN( "CommandLineHandlerSetting::CommandLineHandlerSetting" );


    // --------------------------------------
    //  Check, if any parameter was supplied
    // --------------------------------------
    if ( argc == 1 ) {
      std::cerr << "\n As HAL said in 2001\n\n"
                << "   \"I'm sorry Dave, I don't have enough information.\""
                <<"\n\n So please specify a name for the current "
                << "simulation run!\n\n";
      exit(-1);
    }


    // ------------------------------
    //  Check, if user asks for help
    // ------------------------------
    for ( Integer i = 0; i < argc; i++ ) {
      if ( std::strcmp( argv[i], markerHelp_.c_str() ) == 0 ||
           std::strcmp( argv[i], markerLongHelp_.c_str() ) == 0 ) {
        PrintUsage();
        exit(1);
      }
    }

    // No help wanted, so proceed normally
    Info->StartProgress( "Parsing command line" );


    // ----------------------------------------
    //  Name of simulation run goes separately
    // ----------------------------------------

    // We expect the last command line argument to be the name of
    // the current simulation run. Check that a proper name was
    // provided by the user.
    simName_ = argv[argc-1];


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

      // --schemaPath
      SettingDef( markerSchemaPath_.c_str(),
                  markerLongSchemaPath_.c_str(),
                  Setting::STRING,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY | Setting::EXPLICIT_ASSIGNMENT,
                  helpSchemaPath_.c_str() ),

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

      // --printGrid
      SettingDef( markerPrintGrid_.c_str(),
                  markerLongPrintGrid_.c_str(),
                  Setting::FLAG,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY,
                  helpPrintGrid_.c_str() ),

      // --showEqnMap
      SettingDef( markerShowEqnMap_.c_str(),
                  markerLongShowEqnMap_.c_str(),
                  Setting::FLAG,
                  1,
                  1,
                  Setting::COMMAND_LINE_ONLY,
                  helpShowEqnMap_.c_str() ),

      // --help
      // SettingDef( markerHelp_.c_str(),
      //             markerLongHelp_.c_str(),
      //             Setting::FLAG,
      //             1,
      //             1,
      //             Setting::COMMAND_LINE_ONLY,
      //             helpHelp_.c_str() ),

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
    commandLine_.finalCheck( commandLineParamsSpec );

    // Parallel jobs define numerous additional options for MPI and co
    // so we cannot afford to be picky
#ifndef MpCCI
    if ( pickyResult == false ) {
      PrintUsage();
      (*error) << "Encountered problems while parsing command line!";
      Error( __FILE__, __LINE__ );
    }
#endif

#ifdef DEBUG
    PrintParams( std::cerr );
#endif

    // That's it
    Info->FinishProgress();

  }

}
