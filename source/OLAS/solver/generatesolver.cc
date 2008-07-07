// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_ilupack.hh>
#include <def_use_lapack.hh>
#include <def_use_pardiso.hh>

#include "matvec/matvec.hh"
#include "solver/solver.hh"
#include "solver/generatesolver.hh"
#include "General/exception.hh"

#ifdef USE_LAPACK
#include "external/lapack/lapacklu.hh"
#include "external/lapack/lapackll.hh"
#endif

#ifdef USE_PARDISO
#include "external/pardiso/pardisosolver.cc"
#endif

#ifdef USE_ILUPACK
#include "external/ilupack/Ilupack.cc"
#endif

// include source code for templated solvers
#include "solver/basesolver.cc"
#include "solver/richardson.cc"
#include "solver/cgsolver.cc"
#include "solver/gmres.cc"
#include "solver/minres.cc"
#include "solver/lusolver.cc"
#include "solver/ldlsolver.cc"
#include "solver/diagsolver.cc"

#if 0

#if defined(__INTEL_COMPILER)
#include "algsys/olascomm_impl.hh"
#include "matvec/vector.cc"
#endif

#endif

using CoupledField::Exception;

namespace OLAS {



  // **********************************************************
  //   Macro for generation of a solver object intended for
  //   the case that the solver depends not only on the entry
  //   type, but also on the block size of the matrix.
  // *********************************************************
#define SOLVER_OBJ( entryType, blockSize, solverObjType ) \
if ( entryType == eType && blockSize == bSize ) { \
  retSolver = new solverObjType( params, report ); \
  (*cla) << " GenerateSolver: Generated " \
         << MACRO2STRING( solverObjType ) \
         << " solver" << std::endl;\
}



// ****************************
//   Generate a solver object
// ****************************
BaseSolver* GenerateSolverObject( const BaseMatrix &mat, SolverType solver,
                                  ParamNode* xml, OLAS_Params *params, OLAS_Report *report ){

  BaseSolver *retSolver = NULL;
  MatrixEntryType eType = mat.GetEntryType();
  bool eTypeUnknown = false;

  // extract a solver node if there is one 
  ParamNode* solver_xml = NULL;
  if(xml != NULL && xml->Has("solver")) solver_xml = xml->Get("solver");

  // Branch depending on desired solver
  switch( solver ) {

  case RICHARDSON:
    if ( eType == DOUBLE ) {
      retSolver = New RichardsonSolver<Double>( params, report );
      AssertMem( retSolver, sizeof(RichardsonSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Richardson solver"
      << std::endl;
    }
    else if ( eType == COMPLEX ) {
      retSolver = New RichardsonSolver<Complex>( params, report );
      AssertMem( retSolver, sizeof(RichardsonSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Richardson solver"
      << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case DIAGSOLVER:
    if ( eType == DOUBLE ) {
      retSolver = New DiagSolver<Double>( params, report );
      AssertMem( retSolver, sizeof(DiagSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Diagonal solver"
      << std::endl;
    }
    else if ( eType == COMPLEX ) {
      retSolver = New DiagSolver<Complex>( params, report );
      AssertMem( retSolver, sizeof(DiagSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Diagonal solver"
      << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case CG:
    if(eType == DOUBLE) {
      retSolver = new CGSolver<Double>(solver_xml, params, report );
      (*cla) << " GenerateSolver: Generated real CG solver" << std::endl;
    }
    if(eType == COMPLEX) {
      retSolver = new CGSolver<Complex>(solver_xml, params, report );
      (*cla) << " GenerateSolver: Generated complex CG solver" << std::endl;
    }
    break;

  case GMRES:
    if ( eType == DOUBLE ) {
      retSolver = New GMRESSolver<Double>( params, report );
      AssertMem( retSolver, sizeof(GMRESSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real GMRES solver" << std::endl;
    }
    else if ( eType == COMPLEX ) {
      retSolver = New GMRESSolver<Complex>( params, report );
      AssertMem( retSolver, sizeof(GMRESSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex GMRES solver"
      << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case MINRES:
    if ( eType == DOUBLE ) {
      retSolver = New MINRESSolver<Double>( params, report );
      AssertMem( retSolver, sizeof(MINRESSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real MINRES solver" << std::endl;
    }
    else if ( eType == COMPLEX ) {
      retSolver = New MINRESSolver<Complex>( params, report );
      AssertMem( retSolver, sizeof(MINRESSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex MINRES solver"
      << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case LU_SOLVER:
    if ( eType == DOUBLE ) {
      retSolver = New LUSolver<Double>( params, report );
      AssertMem( retSolver, sizeof(LUSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real "
      << Enum2String( LU_SOLVER ) << std::endl;
    }
    else if ( eType == COMPLEX ) {
      retSolver = New LUSolver<Complex>( params, report );
      AssertMem( retSolver, sizeof(LUSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex "
      << Enum2String( LU_SOLVER ) << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case LDL_SOLVER:

    // Get block size of matrix entries (we need an StdMatrix for this,
    // but, if it's not, an LDL_Solver makes no sense anyhow)
  {
    Integer bSize = 0;
    TRY_CAST {
      ConstRefCast( mat, StdMatrix, stdMat );
      bSize = stdMat.GetBlockSize();
    }
    CATCH_CAST;

    SOLVER_OBJ( DOUBLE,  1, LDLSolver<Double>      );
    SOLVER_OBJ( COMPLEX, 1, LDLSolver<Complex> );
  }
  break;

#ifdef USE_LAPACK
  case LAPACK_LU:
    if ( mat.GetStructureType() != STDMATRIX ) {
      Error( "LAPACK_LU only works with a LAPACK_GBMATRIX!", __FILE__,
             __LINE__ );
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != LAPACK_GBMATRIX ) {
        Error( "LAPACK_LU only works with a LAPACK_GBMATRIX!", __FILE__,
               __LINE__ );
      }
      else {
        retSolver = New Lapack_LU( params, report );
        AssertMem( retSolver, sizeof(Lapack_LU) );
        (*cla) << " GenerateSolver: Generated Lapack_LU solver"
        << std::endl;
      }
    }
    break;

  case LAPACK_LL:
    if ( mat.GetStructureType() != STDMATRIX ) {
      Error( "LAPACK_LL only works with an SCRS_MATRIX!", __FILE__,
             __LINE__ );
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != SPARSE_SYM ) {
        Error( "LAPACK_LL only works with an SCRS_MATRIX!", __FILE__,
               __LINE__ );
      }
      else {
        retSolver = New Lapack_LL( params, report );
        AssertMem( retSolver, sizeof(Lapack_LL) );
        (*cla) << " GenerateSolver: Generated Lapack_LL solver"
        << std::endl;
      }
    }
    break;
#else
  case LAPACK_LU:
    Error( "Compile with USE_LAPACK to enable support for LAPACK_LU solver",
           __FILE__, __LINE__ );
    break;
  case LAPACK_LL:
    Error( "Compile with USE_LAPACK to enable support for LAPACK_LL solver",
           __FILE__, __LINE__ );
    break;
#endif

  case PARDISO:

#ifdef USE_PARDISO

    // Check suitability of matrix
    if ( mat.GetStructureType() != STDMATRIX ) {
      Error( "PardisoSolver only works with (S)CRS_Matrix class!", __FILE__,
             __LINE__ );
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != SPARSE_NONSYM &&
          stdmat.GetStorageType() != SPARSE_SYM  ) {
        Error( "PardisoSolver only works with (S)CRS_Matrix class!",
               __FILE__, __LINE__ );
      }
    }

    if ( eType == DOUBLE ) {
      retSolver = New PardisoSolver<Double>( params, report );
      AssertMem( retSolver, sizeof(PardisoSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Pardiso solver"
      << std::endl;
    }
    if ( eType == COMPLEX ) {
      retSolver = New PardisoSolver<Complex>( params, report );
      AssertMem( retSolver, sizeof(PardisoSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Pardiso solver"
      << std::endl;
    }
#else

    Error( "Compile with USE_PARDISO to enable interface to Pardiso "
           "library", __FILE__, __LINE__ );
#endif
    break;


  case ILUPACK_SOLVER:

#ifdef USE_ILUPACK
  {  
    // Check suitability of matrix
    if (mat.GetStructureType() != STDMATRIX)
      throw Exception("Ilupack only works with (S)CRS_Matrix class!");

    const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
    if(stdmat.GetStorageType() != SPARSE_NONSYM && stdmat.GetStorageType() != SPARSE_SYM  )
      throw Exception("Ilupack only works with (S)CRS_Matrix class!"); 

    if(eType == DOUBLE) {
      retSolver = new Ilupack<Double>(solver_xml, report, eType);
      (*cla) << " GenerateSolver: Generated real ilupack solver" << std::endl;
    }
    if(eType == COMPLEX) {
      retSolver = new Ilupack<Complex>(solver_xml, report, eType);
      (*cla) << " GenerateSolver: Generated complex ilupack solver" << std::endl;
    }   

  }   
#else
  throw Exception("Compile with USE_ILUPACK to enable interface to ilupack");
#endif
  break;

  default:
    throw Exception("GenerateSolver: Request for unknown solver type!");
  }

  // Check for unsupported matrix entry type
  if (retSolver == NULL ) EXCEPTION("unhandled type " << eType);

  // Force instantiation of member functions of templated matrix classes
  // Note that the InstantiatePublicMethods() method itself should never
  // actually be called, but the compiler must not notice this. Thus, the
  // akward if condition.
  bool neverTrue = (eType == NOENTRYTYPE);
  if ( neverTrue ) {
    BaseMatrix * dummy;
    retSolver->InstantiatePublicMethods( *dummy );
  }

  return retSolver;
}

}
