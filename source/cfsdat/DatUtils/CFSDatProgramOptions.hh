// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ProgramOptions.hh
 *       \brief    <Description>
 *
 *       \date     Nov 19, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef PROGRAMOPTIONS_HH_
#define PROGRAMOPTIONS_HH_

#include <boost/program_options.hpp>
#include <filesystem>
#include "DataInOut/ProgramOptions.hh"


namespace CFSDat
{

class CFSDatProgramOptions : public CoupledField::ProgramOptions
{

public:

  //! Standard constructor
  CFSDatProgramOptions( CoupledField::Integer argc, const char **argv );

  virtual ~CFSDatProgramOptions();

  //! Gather information from commandline and environment

  //! This method triggers the reading of information from the command line
  //! and the environment.
  virtual void ParseData();

  /** This gives the head line of CFSDat printed to cout */
  virtual void GetHeaderString(std::ostream& out);

  /** Get number of threads for CFS supplied on the command line */
  virtual CoupledField::UInt GetNumThreads() const;


private:
   // =======================================================================
   // INTERNAL HELPER METHODS
   // =======================================================================

   //@{ \name Additional internal helper methods

   //! Helper function for mapping environment variables to internal names
   std::string EnvironmentNameMapper( const std::string& envVarName );

   //@}
};

}

#endif /* PROGRAMOPTIONS_HH_ */
