// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/**********************************************************
 * Separate file for the preprocessor flags and definitions
 * for AMG.
 **********************************************************/

#ifndef _AMG_PPFLAGS_HH_
#define _AMG_PPFLAGS_HH_

#define  AMG_USES_MPI_CLOCK

/*****************************/
#define AMG_IMPORT_PREFIX  "import/" // "lasexport/las"
#define AMG_EXPORT_PREFIX  "export/" // "tools/debug/import/"


/**********************************************************
 * some more detailed flags for import, export, split trace
 **********************************************************/

#ifdef AMG_IMPORT_HIERARCHY
//// more detailed import settings
#  define TOPOLOGY_IMPORT_CF_SPLITTING
#  define  PREMATRIX_ENABLE_IMPORT
#endif

#ifdef AMG_EXPORT_HIERARCHY
#  define AMG_EXPORT_COMPONENTS
#  define AMG_EXPORT_CYCLE
#  define TOPOLOGY_EXPORT_CF_SPLITTING
#endif

// some definitions for the CF-splitting trace
#ifdef AMG_TRACE_CF_SPLITTING
#  define AMG_CF_TRACE_TYPE    unsigned char
#  define AMG_CF_TRACE_BINARY
   // fix name of the trace file
#  ifdef  AMG_CF_TRACE_BINARY
#    define  AMG_CF_SPLITTING_TRACE_FILE_NAME  "CFtrace.bin"
#  else
#    define  AMG_CF_SPLITTING_TRACE_FILE_NAME  "CFtrace.ascii"
#  endif
#endif



/**********************************************************
 * DEBUGGING
 **********************************************************/

#ifdef DEBUG_MULTIGRID ///////////////
                                    //
#ifndef DEBUG_AMG                   //
#define  DEBUG_AMG                  //
#endif                              //
#ifndef DEBUG_HIERARCHYLEVEL        //
#define  DEBUG_HIERARCHYLEVEL       //
#endif                              //
#ifndef DEBUG_TOPOLOGY              //
#define  DEBUG_TOPOLOGY             //
#endif                              //
#ifndef DEBUG_DEPENDENCYGRAPH       //
#define  DEBUG_DEPENDENCYGRAPH      //
#endif                              //
#ifndef DEBUG_PREMATRIX             //
#define  DEBUG_PREMATRIX            //
#endif                              //
#ifndef DEBUG_SMOOTHER              //
#define  DEBUG_SMOOTHER             //
#endif                              //
#ifndef DEBUG_GAUSSSEIDEL           //
#define  DEBUG_GAUSSSEIDEL          //
#endif                              //
#ifndef DEBUG_JACOBI                //
#define  DEBUG_JACOBI               //
#endif                              //
                                    //
#endif // DEBUG_MULTIGRID ////////////



/**********************************************************/
#endif // _AMG_PPFLAGS_HH_
