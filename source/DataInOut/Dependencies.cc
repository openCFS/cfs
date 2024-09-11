#include <def_use_arpack.hh>
#include <def_use_eigen.hh>
#include <def_use_blas.hh>
#include <def_use_cgal.hh>
#include <def_use_cgns.hh>
#include <def_use_embedded_python.hh>
#include <def_use_ensight.hh>
#include <def_use_feast.hh>
#include <def_use_flann.hh>
#include <def_use_ghost.hh>
#include <def_use_gidpost.hh>
#include <def_use_hwloc.hh>
#include <def_use_ipopt.hh>
#include <def_use_libfbi.hh>
#include <def_use_libxml2.hh>
#include <def_use_lis.hh>
#include <def_use_metis.hh>
#include <def_use_mpi.hh>
#include <def_use_openmp.hh>
#include <def_use_pardiso.hh>
#include <def_use_petsc.hh>
#include <def_use_phist_cg.hh>
#include <def_use_phist_ev.hh>
#include <def_use_sgp.hh>
#include <def_use_sgpp.hh>
#include <def_use_snopt.hh>
#include <def_use_scpip.hh>
#include <def_use_suitesparse.hh>
#include <def_use_superlu.hh>
#include <def_use_unv.hh>
#include <def_use_xerces.hh>
#include <def_xmlschema.hh>
#include <def_xmlschema.hh>
#include <def_config.hh>

#ifdef USE_MKL
#include <mkl_service.h>
#ifndef mkl_get_version
#define MKL_Get_Version MKLGetVersion
#define MKL_Free_Buffers MKL_FreeBuffers
#endif
#endif


#ifdef USE_OPENBLAS
#include <openblas/openblas_config.h>
#endif

#include <zlib.h>

#include <H5public.h>
#include <H5Ppublic.h>

#ifdef USE_IPOPT
#include <coin-or/IpoptConfig.h>
#endif

#ifdef USE_METIS
#include <metis.h>
#endif

#ifdef USE_SUITESPARSE
#include <cholmod.h>
#include <umfpack.h>
#include <amd.h>
#endif

#ifdef USE_CGAL
#include <CGAL/version.h>
#endif

#ifdef USE_EIGEN
#include <Eigen/src/Core/util/Macros.h>
#endif

#ifdef USE_FLANN
#include <flann/flann.hpp>
#endif

#ifdef USE_XERCES
#include <xercesc/util/XercesVersion.hpp>
#endif

#ifdef USE_LIBXML2
#include <libxml/xmlversion.h>
#endif

#ifdef USE_GIDPOST
#include <gidpost.h>
#endif

#ifdef USE_CGNS
#include <cgnslib.h>
#endif

#ifdef USE_EMBEDDED_PYTHON
#include <patchlevel.h>
#endif

#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

#include <boost/version.hpp>
#include <fstream>
#include <muParserBase.h>

#if _WIN32
  #undef max
#endif

#include "DataInOut/Dependencies.hh"
#include "General/Environment.hh"

namespace CoupledField {

Dependencies::Dependencies()
{
  ReadSetting();
}


bool Dependencies::IsDistributable() const
{
  for(auto& dep : data)
    if(dep.active && !IsDistributable(dep.lic))
      return false;

  return true;
}


bool Dependencies::HasLicense(Dependencies::License test) const
{
  for(auto& dep : data) {
    assert(dep.lic != NOT_SET);
    if(dep.active && dep.lic == test)
      return true;
  }
  return false;
}

void Dependencies::Dump() const
{
  for(auto& dep : data)
    std::cout << dep.name << ", " << dep.cmake << ", " << dep.lic << ", " << dep.active << ", " << dep.version
              << ", " << dep.comment << ", " << dep.date << std::endl;
}


bool Dependencies::WriteCMakeUSE(const string& filename)
{
  std::ofstream out(filename.c_str());
  if(!out.is_open()) {
    std::cerr << "cannot open '" << filename << "' to write CMake dependencies.";
    return false;
  }

  out << "# written by cfs via --dependencies option.\n";
  out << "# same information as with --version option.\n";
  for(auto& dep : data)
    if(dep.IsSwitchable())
      out << "set(" << dep.cmake << " " << (dep.active ? "ON" : "OFF") << ")\n";

  std::cout << "wrote CMake file '" << filename << "' with dependencies\n";

  return true; // auto close out
}

void Dependencies::ReadSetting()
{
  std::stringstream ss;

  // parallel

  Dependency omp("OpenMP", "USE_OPENMP", EASY);
#ifdef USE_OPENMP
  // for some cmake versions, OpenMP_CXX_VERSION is not set
  #ifdef OpenMP_CXX_VERSION
    omp.SetVersion(OpenMP_CXX_VERSION);
  #endif  
  ss.clear();
  ss << "OMP_NUM_THREADS=" << (getenv("OMP_NUM_THREADS") != NULL ? getenv("OMP_NUM_THREADS") : "-") << ", ";
  if(HasBlasThreadsEnvVariable()) { // false for openblas (would be again OMP_NUM_THREADS) or netlib (is serial)
    const char* otherenv = GetBlasThreadsEnvVariable(); // MKL_NUM_THREADS or VECLIB_MAXIMUM_THREADS
    ss << otherenv << "=" << (getenv(otherenv) != NULL ? getenv(otherenv) : "-") << ", ";
  }
  ss << "CFS_NUM_THREADS=" << CFS_NUM_THREADS;
  omp.comment = ss.str();
#endif
  data.Push_back(omp);

  Dependency mpi("mpi", "USE_MPI", NOT_KNOWN);
#ifdef USE_MPI
  mpi.active = true;
  ss.clear();
  ss << "Base=" << CFS_MPI_BASE << " ";
  ss << "CXX=" << CFS_MPI_CXX_COMPILER << " ";
  ss << "FC=" << CFS_MPI_Fortran_COMPILER << " ";
  ss << "mpiexec=" << CFS_MPI_BIN << "/mpiexec";
  mpi.comment = ss.str();
#endif
  data.Push_back(mpi);

  Dependency petsc("PETSc", "USE_PETSC", BSD);
#ifdef USE_PETSC
  petsc.SetVersion(CFS_PETSC_VERSION);
#endif
  data.Push_back(petsc);


  Dependency arpack("ARPACK", "USE_ARPACK", BSD);
#ifdef USE_ARPACK
  arpack.SetVersion(ARPACK_VER);
#endif
  data.Push_back(arpack);

  Dependency eigen("EIGEN", "USE_EIGEN", EASY);
#ifdef USE_EIGEN
  eigen.SetVersion(EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
#endif
  data.Push_back(eigen);

  Dependency feast("FEAST", "USE_FEAST", BSD);
#ifdef USE_FEAST
  feast.SetVersion(FEAST_VER);
#endif
  data.Push_back(feast);

  Dependency mkl("MKL", "USE_MKL", ISSL);
#ifdef USE_MKL
  CFSMKLVersion ver;
  MKL_Get_Version(reinterpret_cast<MKLVersion*>(&ver));
  MKL_Free_Buffers();
  mkl.SetVersion(ver.MajorVersion,ver.MinorVersion,ver.BuildNumber);
#endif
  data.Push_back(mkl);

  Dependency ob("OpenBLAS", "USE_OPENBLAS", BSD);
#ifdef USE_OPENBLAS
  ob.SetVersion(OPENBLAS_VER);
  ob.comment = "Core: " + string(OPENBLAS_CHAR_CORENAME);
#endif
  data.Push_back(ob);

  Dependency acc("Accelerate", "USE_ACCELERATE", EASY); // dynamically link to Apple's Accelerate Framework
#ifdef USE_ACCELERATE
  acc.active = true;
  acc.comment = "set VECLIB_MAXIMUM_THREADS";
#endif
  data.Push_back(acc);

  // All other blas have their own lapack. Netlib contains blas/lapack and is the slow reference implementation
  Dependency nl("Netlib", "USE_NETLIB", BSD);
#ifdef USE_NETLIB
  nl.SetVersion(NETLIB_VER);
#endif
  data.Push_back(nl);

  // liner system solver

  Dependency pardiso("PARDISO", "USE_PARDISO", ISSL); // would be wrong if not from MKL
#ifdef USE_PARDISO
  pardiso.active = true;
  pardiso.comment = CFS_PARDISO;
  if(!mkl.active)
    pardiso.lic = CLOSED;
#endif
  data.Push_back(pardiso);

  // we can configure variant with or without gpl. With gpl cholmod is faster for highly populated systems
  Dependency suite("SuiteSparse", "USE_SUITESPARSE", GPL);
 #ifdef USE_SUITESPARSE
  suite.SetVersion(SUITESPARSE_MAIN_VERSION,SUITESPARSE_SUB_VERSION,SUITESPARSE_SUBSUB_VERSION);
  suite.date = SUITESPARSE_DATE;

  data.Push_back(Dependency("AMD", "", "", LGPL, "part of SuiteSparse", AMD_DATE));
  data.Last().SetVersion(AMD_MAIN_VERSION, AMD_SUB_VERSION, AMD_SUBSUB_VERSION);
  data.Push_back(Dependency("CHOLMOD", "", "", GPL, "part of SuiteSparse", CHOLMOD_DATE));
  data.Last().SetVersion(CHOLMOD_MAIN_VERSION, CHOLMOD_SUB_VERSION, CHOLMOD_SUBSUB_VERSION);
  data.Push_back(Dependency("UMFPACK", "", "", GPL, "part of SuiteSparse", UMFPACK_DATE));
  data.Last().SetVersion(UMFPACK_MAIN_VERSION, UMFPACK_MAIN_VERSION, UMFPACK_MAIN_VERSION);
 #endif
  data.Push_back(suite);

  Dependency lis("LIS", "USE_LIS", BSD);
#ifdef USE_LIS
  lis.SetVersion(LIS_VER);
#endif
  data.Push_back(lis);

  Dependency superlu("SuperLU", "USE_SUPERLU", EASY);
#ifdef USE_SUPERLU
  superlu.SetVersion(SUPERLU_VER);
#endif
  data.Push_back(superlu);

// for ghost only
#ifdef BUILD_HWLOC
  Dependency hwloc("hwloc", "BUILD_HWLOC", BSD);
  hwloc.SetVersion(HWLOC_VER);
  data.Push_back(hwloc);
#endif

// only relevant for phist
#ifdef BUILD_GHOST
  Dependency ghost("ghost", "BUILD_GHOST", NOT_KNOWN);
  ghost.SetVersion(GHOST_REV);
  ghost.comment = GHOST_SOURCE;
  data.Push_back(ghost);
#endif

  Dependency phist_cg("phist_cg", "USE_PHIST_CG", NOT_KNOWN);
#ifdef USE_PHIST_CG
  phist_cg.SetVersion(PHIST_REV);
  phist_cg.comment = PHIST_SOURCE;
#endif
  data.Push_back(phist_cg);

  Dependency phist_ev("phist_ev", "USE_PHIST_EV", NOT_KNOWN);
#ifdef USE_PHIST_EV
  phist_ev.SetVersion(PHIST_REV);
  phist_ev.comment = PHIST_SOURCE;
#endif
  data.Push_back(phist_ev);


  // technical stuff

  Dependency metis("METIS", "USE_METIS", APACHE2);
#ifdef USE_METIS
  metis.SetVersion(METIS_VER_MAJOR, METIS_VER_MINOR, METIS_VER_SUBMINOR);
#endif
  data.Push_back(metis);

  Dependency cgal("CGAL", "USE_CGAL", GPL);
#ifdef USE_CGAL
  cgal.SetVersion(QUOTEME(CGAL_VERSION));
#endif
  data.Push_back(cgal);

  Dependency fbi("libfbi","USE_LIBFBI",NOT_KNOWN);
#ifdef USE_LIBFBI
  fbi.active = true;
#endif
  data.Push_back(fbi);

  Dependency flann("flann", "USE_FLANN", NOT_KNOWN);
#ifdef USE_FLANN
  flann.SetVersion(FLANN_VERSION_);
#endif
  data.Push_back(flann);

  Dependency expr("XPRT", "USE_EXPRESSION_TEMPLATES", NOT_KNOWN);
#ifdef USE_EXPRESSION_TEMPLATES
  expr.active = true;
#endif
  data.Push_back(expr);

  // I/O

  Dependency gid("GiDpost", "USE_GIDPOST", EASY);
#ifdef USE_GIDPOST
  metis.SetVersion(GIDPOST_VERSION);
#endif
  data.Push_back(gid);

  Dependency vtk("VTK", "USE_ENSIGHT", BSD);
#ifdef USE_ENSIGHT
  vtk.SetVersion(CFS_VTK_VERSION);
  vtk.comment = "VTK enabled by USE_ENSIGHT";
#endif
  data.Push_back(vtk);

  Dependency cgns("CGNS", "USE_CGNS", LGPL); // close but not exactly -> https://cgns.github.io/CGNS_docs_current/charter/license.html
#ifdef USE_CGNS
  cgns.SetVersion(CGNS_VERSION/1000, (CGNS_VERSION%1000)/100, (CGNS_VERSION%100)/10);
  ss.clear();
#endif
  data.Push_back(cgns);

  Dependency xerces("Xerces", "USE_XERCES", APACHE2);
#ifdef USE_XERCES
  xerces.SetVersion(XERCES_FULLVERSIONDOT);
#endif
  data.Push_back(xerces);

  Dependency libxml2("libxml2", "USE_LIBXML2", MIT);
#ifdef USE_LIBXML2
  libxml2.SetVersion(LIBXML_DOTTED_VERSION);
#endif
  data.Push_back(libxml2);
  assert(xerces.active || libxml2.active);


  // optimizers

  Dependency scpip("SCPIP", "USE_SCPIP",CLOSED);
#ifdef USE_SCPIP
  scpip.SetVersion(SCPIP_VER);
#endif
  data.Push_back(scpip);

  Dependency snopt("SNOPT", "USE_SNOPT", COMMERCIAL);
#ifdef USE_SNOPT
  snopt.SetVersion(SNOPT_VER);
#endif
  data.Push_back(snopt);

  // https://projects.coin-or.org/Ipopt/wiki/FAQ
  Dependency ipopt("IPOPT", "USE_IPOPT", ECLIPSE);
#ifdef USE_IPOPT
  ipopt.SetVersion(IPOPT_VERSION);
#endif
  data.Push_back(ipopt);

  // requires a system python to be linked to
  Dependency python("Python", "USE_EMBEDDED_PYTHON", EASY);
#ifdef USE_EMBEDDED_PYTHON
  python.SetVersion(PY_VERSION);
#endif
  data.Push_back(python);

  // was developed from FAU
  Dependency sgp("SGP", "USE_SGP", MIT);
#ifdef USE_SGP
  sgp.SetVersion(SGP_VER);
#endif
  data.Push_back(sgp);

  // to be replaced by a open source version when it is available
  Dependency sgpp("SGPP", "USE_SGPP", COMMERCIAL);
#ifdef USE_SGPP
  sgpp.SetVersion(SGPP_VER);
#endif
  data.Push_back(sgpp);

  // build options are not for the version output to be used via WriteCMakeUSE() in the testsuite
  Dependency cfsdat("cfsdat", "BUILD_CFSDAT", CFS);
#ifdef BUILD_CFSDAT
  cfsdat.active = true;
#endif
  data.Push_back(cfsdat);

  Dependency cfstool("cfstool", "BUILD_CFSTOOL", CFS);
#ifdef BUILD_TOOL
  cfstool.active = true;
#endif
  data.Push_back(cfstool);

  Dependency testing("TESTING", "BUILD_TESTING", CFS);
#ifdef BUILD_TESTING
  testing.active = true;
#endif
  data.Push_back(testing);

  Dependency unittests("Unit-Tests", "BUILD_UNIT_TESTS", CFS);
#ifdef BUILD_UNIT_TESTS
  unittests.active = true;
#endif
  data.Push_back(unittests);

  // https://www.boost.org/users/license.html
  Dependency boost("boost", "", BOOST);
  boost.SetVersion(BOOST_VERSION / 100000, BOOST_VERSION / 100 % 1000, BOOST_VERSION % 100);
  data.Push_back(boost);

  // https://de.wikipedia.org/wiki/Zlib-Lizenz
  data.Push_back(Dependency("zlib", "", zlibVersion(), EASY));

  // https://github.com/beltoforion/muparser
  data.Push_back(Dependency("muparser", "", MUP_VERSION, BSD, ""));

  Dependency hdf5("hdf5", "", BSD);
  hdf5.SetVersion(H5_VERS_MAJOR,H5_VERS_MINOR,H5_VERS_RELEASE);
  data.Push_back(hdf5);

}

// Allows the binary a free distribution? Not for GPL, COMMERICAL and CLOSED
bool Dependencies::IsDistributable(License lic)
{
  switch(lic)
  {
  case GPL:
  case COMMERCIAL:
  case CLOSED:
    return false;
  default:
    break;
  }
  return true;
}

Dependencies::Dependency::Dependency(const string& name, const string& cmake, License lic)
{
  this->name = name;
  this->cmake = cmake;
  this->lic = lic;
}

Dependencies::Dependency::Dependency(const string& name, const string& cmake, const string& version, License lic, const string comment, const string date)
{
  this->name = name;
  this->cmake = cmake;
  this->version = version;
  this->lic = lic;
  this->comment = comment;
  this->date = date;

  this->active = true; // when we have the version, we are active
}

bool Dependencies::Dependency::IsSwitchable() const
{
  return cmake.find("USE_") != string::npos || IsBuildOption();
}

bool Dependencies::Dependency::IsBuildOption() const
{
  return cmake.find("BUILD_") != string::npos;
}


void Dependencies::Dependency::SetVersion(const string& version, const string& sub, const string& minor)
{
  SetVersion(version + "." + sub + "." + minor);
}
void Dependencies::Dependency::SetVersion(int version, int sub, int minor)
{
  std::stringstream ss;
  ss << version;
  if(sub >= 0)
    ss << "." << sub;
  if(minor >= 0)
    ss << "." << minor;
  SetVersion(ss.str());
}

} // end of namespace

