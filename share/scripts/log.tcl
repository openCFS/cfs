# Sample tcl scripting file demonstrating how to switch on/off the logging
# streams defined in CFS++ at runtime.
# author: Andreas Hauck
# Last modified: 21.12.2006

proc CFS_Init { problemname } {

    # Configuration of log behavior can be influenced by method
    #   cfs logConf ....
    #
    # There are two different additional parameters:
    #
    # setLogLevel <streamName> <level>
    #
    #   Defines the level of detail
    #   to be activated for a logging stream.
    #
    #   <streamName> Name of the logging stream to be activated.
    #                Can also be "*" to capture all streams
    #   <level>      Defines the logging level to be captured.
    #                Possible values are:
    #                  trace
    #                  trace2
    #                  dbg
    #                  dbg2
    #                  dbg3
    #                  all
    #
    # addAppender <streamName> <output> <name>
    #
    #   Route the output of a logging stream to a given destination.
    #
    #   <streamName> Name of the logging stream to be activated.
    #                Can also be "*" to capture all streams
    #   <output>     Destination type for the stream. Can be one of
    #                  file: Output is appended in an ASCII-file
    #                  cout: Output is simply printed on the terminal
    #
    #   <name>      If <output> is file the value denotes the filename.
    #               Otherwise the value is omitted.

    # 1) Set log-level for eqnMap class to DBG3 and print it to file
    #    eqns.txt
    cfs logConf setLogLevel eqnMap dbg3    
    cfs logConf addAppender eqnMap file "eqns.txt"

    # 2) Set log-level for all pdes to TRACE and print to cout
    cfs logConf setLogLevel pde trace
    cfs logConf addAppender pde cout ""
    
}
