// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef HELPER_FUNCS_HEADER_2013
#define HELPER_FUNCS_HEADER_2013

#include <string>

namespace CFSTool
{

  // Here we define our hard-coded node and element limits for testsuite
  // tests. These limits come into effect, if the reference file of a
  // Diff() operation has the extension '.h5ref' or if the --checkLimits
  // command line option is set to true.
  const static UInt SOFT_NODE_LIMIT = 800;
  const static UInt HARD_NODE_LIMIT = 1000;
  const static UInt SOFT_ELEM_LIMIT = 1000;
  const static UInt HARD_ELEM_LIMIT = 6000;

  // Here we print out the hard-coded node and element limit in response
  // to the --printLimits command line option.
  void PrintLimits(const std::string& type);
  
  shared_ptr<SimInput> GetInputReader(const std::string& fileName,
                                      const PtrParamNode& param,
                                      const PtrParamNode& info);

  shared_ptr<SimOutput> GetOutputWriter( const std::string& fileName,
                                         const PtrParamNode& param,
                                         const PtrParamNode& info );

  PtrParamNode GetParamNodeFromHDF5( shared_ptr<SimInput>& inputHDF5,
                                     const std::string& xmlFile );

  double RadPhase( const Complex& c );

  bool CheckReaderCapabilities(const std::string& readerDescription,
                               shared_ptr<SimInput>& reader,
                               std::string& result);

  /** Initialize static Enums.
   * todo: do better once - Fabian */
  void InitEnums();

  void SetFreeCoord(const PtrParamNode& param,
                    std::string coordSysId="default",
                    std::string node_name="averageDomain");
  

}

#endif
