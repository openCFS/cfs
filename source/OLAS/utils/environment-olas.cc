// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <sstream>

#include "utils/environment.hh"
#include "utils/tools.hh"

namespace OLAS {

  // --------------------------------------------
  //  Implementation of enum conversion routines 
  // --------------------------------------------

  // Default implementation
  template <class ENUMTYPE> 
  std::string Enum2String( const ENUMTYPE &in ) {
    ENTER_IFCN( "Enum2String" );
    (*error) << "Enum2String not implemented for the desired type"
         << " of enumeration";
    Error( __FILE__, __LINE__ );
  }

  // Specialisation for SolverType
  template<> 
  std::string Enum2String<SolverType>( const SolverType &in ) {

    ENTER_IFCN( "Enum2String" );

    std::string out;

    switch( in ) {
    case NOSOLVER:
      out = "no solver";
      break;
    case RICHARDSON:
      out = "Richardson";
      break;
    case DIAGSOLVER:
      out = "diagsolver";
      break;
    case CG:
      out = "cg";
      break;
    case GMRES:
      out = "gmres";
      break;
    case MINRES:
      out = "minres";
      break;
    case SYMMLQ:
      out = "symmlq";
      break;
    case LAPACK_LU:
      out = "lapackLU";
      break;
    case LAPACK_LL:
      out = "lapackLL";
      break;
    case LU_SOLVER:
      out = "directLU";
      break;
    case LDL_SOLVER:
      out = "directLDL";
      break;
    case LDL_SOLVER2:
      out = "directLDL2";
      break;
    case PARDISO:
      out = "pardiso";
      break;
    case ILUPACK_SOLVER:
      out = "ilupack";
      break;



    default:
      (*error) << "No string value found for the specified value of the "
           << "enumeration datatype SolverType.\n"
           << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }

    return out;
  }

  // Specialisation for EigenSolverType
  template<> 
  std::string Enum2String<EigenSolverType>( const EigenSolverType &in ) {
    
    ENTER_IFCN( "Enum2String" );
    
    std::string out;
    
    switch( in ) {
    case NOEIGENSOLVER:
      out = "no eigensolver";
      break;
    case ARPACK:
      out = "arpack";
      break;
    case SUBSPACE:
      out = "subspace";
      break;
    default:
      (*error) << "No string value found for the specified value of the "
               << "enumeration datatype EigenSolverType.\n"
               << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    
    return out;
  } 


  // Specialisation for PrecondType
  template<> 
  std::string Enum2String<PrecondType>( const PrecondType &in ) {

    ENTER_IFCN( "Enum2String" );

    std::string out;

    switch( in ) {
    case NOPRECOND:
      out = "no precond";
      break;
    case ID:
      out = "Id";
      break;
    case MG:
      out = "MG";
      break;
    case JACOBI:
      break;
    case SSOR:
      out = "SSOR";
      break;
    case ILU0:
      out = "ILU0";
      break;
    case ILUTP:
      out = "ILUTP";
      break;
    case ILUK:
      out = "ILUK";
      break;
    case ILDL0:
      out = "ILDL0";
      break; 
    case ILDLK:
      out = "ILDLK";
      break;
    case ILDLTP:
      out = "ILDLTP";
      break;
    case ILDLCN:
      out = "ILDLCN";
      break;
    case IC0:
      out = "IC0";
      break;

    default:
      (*error) << "No string value found for the specified value of the "
           << "enumeration datatype PrecondType.\n"
           << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }


  // Specialisation for MatrixStorageType
  template<> 
  std::string Enum2String<MatrixStorageType>( const MatrixStorageType &in ) {

    ENTER_IFCN( "Enum2String" );

    std::string out;

    switch( in ) {
    case NOSTORAGETYPE:
      out = "noStorageType";
      break;
    case SPARSE_SYM:
      out = "sparseSym";
      break;
    case SPARSE_NONSYM:
      out = "sparseNonSym";
      break;
    case SKYLINE_SYM:
      out = "skylineSym";
      break;
    case SKYLINE_NONSYM:
      out = "skylineNonSym";
      break;
    case DIAG:
      out = "diag";
      break;
    case LAPACK_GBMATRIX:
      out = "lapackGBMatrix";
      break;

    default:
      (*error) << "No string value found for the specified value of the "
           << "enumeration datatype MatrixStorageType.\n"
           << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }


  // Specialisation for MatrixEntryType
  template<> 
  std::string Enum2String<MatrixEntryType>( const MatrixEntryType &in ) {

    ENTER_IFCN( "Enum2String" );

    std::string out;

    switch( in ) {
    case NOENTRYTYPE:
      out = "noEntryType";
      break;
    case INTEGER:
      out = "integer";
      break;
    case DOUBLE:
      out = "double";
      break;
    case COMPLEX:
      out = "complex";
      break;
    case F77REAL8:
      out = "F77REAL8";
      break;
    case F77COMPLEX16:
      out = "F77COMPLEX16";
      break;

    default:
      (*error) << "No string value found for the specified value of the "
           << "enumeration datatype MatrixEntryType.\n"
           << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }


  // Specialisation for StopCritType
  template<> 
  std::string Enum2String<StopCritType>( const StopCritType &in ) {

    ENTER_IFCN( "Enum2String" );

    std::string out;

    switch( in ) {
    case NOSTOPCRITTYPE:
      out = "no stopping criterion";
      break;
    case ABSNORM:
      out = "absNorm";
      break;
    case RELNORM_RHS:
      out = "relNormRHS";
      break;
    case RELNORM_RES0:
      out = "relNormRes0";
      break;

    default:
      (*error) << "No string value found for the specified value of the "
           << "enumeration datatype StopCritType.\n"
           << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }

  // Specialisation for ReorderingType
  template<>
  std::string Enum2String<ReorderingType>( const ReorderingType &in ) {

    std::string out;
 
    switch( in ) {
    case NOREORDERING:
      out = "noReordering";
      break; 
    case SLOAN:
      out = "Sloan";
      break;
    case METIS:
      out = "Metis";
      break;
    case MINIMUM_DEGREE:
      out = "minimumDegree";
      break;
    case NESTED_DISSECTION:
      out = "nestedDissection";
      break;
    default:  
      (*error) << "No string value found for the specified value of the "
           << "enumeration datatype ReorderingType.\n"
           << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }

  // Specialisation for AMG interpolation type
  template<>
  std::string Enum2String<AMGInterpolationType>(
                    const AMGInterpolationType &in ) {

    std::string out;

    switch( in ) {
    case AMG_INTERPOLATION_CONSTANT:
      out = "constant";
      break;
    case AMG_INTERPOLATION_SIMPLE_WEIGHTED:
      out = "simpleWeighted";
      break;
    case AMG_INTERPOLATION_SMOOTHED_SCALING:
      out = "smoothedScaling";
      break;
    case AMG_INTERPOLATION_DEVELOP:
      out = "develop";
      break;
    default:
      (*error) << "No string value found for the specified value of the "
               << "enumeration datatype AMGInterpolationType.\n"
               << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }

  // specialisation for AMG smoother type
  template<>
  std::string Enum2String<AMGSmootherType>(
                    const AMGSmootherType &in ) {

    std::string out;

    switch( in ) {
        case AMG_SMOOTHER_GAUSSSEIDEL:
          out = "GaussSeidel";
          break;
        case AMG_SMOOTHER_DAMPED_JACOBI:
          out = "Jacobi";
          break;
        default:
          (*error)
          << "No string value found for the specified value of the "
          << "enumeration datatype AMGSmootherType.\n"
          << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }

  // Specialisation for FEMatrixType
  template<>
  std::string Enum2String<FEMatrixType>(
                    const FEMatrixType &in ) {

    std::string out;

    switch( in ) {
    case NOTYPE:
      out = "no_fe_matrix";
      break;
    case SYSTEM:
      out = "system";
      break;
    case STIFFNESS:
      out = "stiffness";
      break;
    case DAMPING:
      out = "damping";
      break;
    case CONVECTION:
      out = "convection";
      break;
    case MASS:
      out = "mass";
      break;
    default:
      (*error) << "No string value found for the specified value of the "
               << "enumeration datatype FEMatrixTypeType.\n"
               << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }

  // Specialisation for IDBCType
  template<>
  std::string Enum2String<IDBCType>( const IDBCType &in ) {

    std::string out;

    switch( in ) {

    case IDBC_NOTYPE:
      out = "notype";
      break;
    case IDBC_ELIMINATION:
      out = "elimination";
      break;
    case IDBC_PENALTY:
      out = "penalty";
      break;
    default:
      (*error) << "No string value found for the specified value of the "
               << "enumeration datatype IDBCType.\n"
               << "Seems to indicate a missing case implementation!";
      Error( __FILE__, __LINE__ );
    }
    return out;
  }

  // -----------------------------------------------
  //  Implementation of string conversion routines 
  // -----------------------------------------------

  // Default implementation
  template <class ENUMTYPE> 
  void String2Enum( const std::string &in, ENUMTYPE &out ) {
    ENTER_IFCN( "String2Enum" );
    (*error) << "String2Enum not implemented for the desired type"
         << " of enumeration";
    Error( __FILE__, __LINE__ );
  }

  // Specialisation for StopCritType
  template<> 
  void String2Enum<StopCritType>( const std::string &in, StopCritType &out ) {

    if ( in == "noStopCritType" ) {
      out = NOSTOPCRITTYPE;
    }
    else if ( in == "absNorm" ) {
      out = ABSNORM;
    }
    else if ( in == "relNormRHS" ) {
      out = RELNORM_RHS;
    }
    else if ( in == "relNormRes0" ) {
      out = RELNORM_RES0;
    }
    else {
      (*error) << "No enumeration value found in StopCritType for '"
           << in << "'\n A missing case implementation?";
      Error( __FILE__, __LINE__ );
    }
  }

  // Specialisation for MatrixStorageType
  template<> 
  void String2Enum<MatrixStorageType>( const std::string &in,
                       MatrixStorageType &out ) {

    if ( in == "noStorageType" ) {
      out = NOSTORAGETYPE;
    }
    else if ( in == "sparseSym" ) {
      out = SPARSE_SYM;
    }
    else if ( in == "sparseNonSym" ) {
      out = SPARSE_NONSYM;
    }
    else if ( in == "skylineSym" ) {
      out = SKYLINE_SYM;
    }
    else if ( in == "skylineNonSym" ) {
      out = SKYLINE_NONSYM;
    }
    else if ( in == "diag" ) {
      out = DIAG;
    }
    else if ( in == "lapackGBMatrix" ) {
      out = LAPACK_GBMATRIX;
    }
    else {
      (*error) << "No enumeration value found in MatrixStorageType for '"
           << in << "'\n A missing case implementation?";
      Error( __FILE__, __LINE__ );
    }
  }

  // Specialisation for SolverType
  template<> 
  void String2Enum<SolverType>( const std::string &in, SolverType &out ) {

    if ( in == "no solver" ) {
      out = NOSOLVER;
    }
    else if ( in == "Richardson" ) {
      out = RICHARDSON;
    }
    else if ( in == "diagsolver" ) {
      out = DIAGSOLVER;
    }
    else if ( in == "cg" ) {
      out = CG;
    }
    else if ( in == "gmres" ) {
      out = GMRES;
    }
    else if ( in == "minres" ) {
      out = MINRES;
    }
    else if ( in == "symmlq" ) {
      out = SYMMLQ;
    }
    else if ( in == "lapackLL" ) {
      out = LAPACK_LL;
    }
    else if ( in == "directLU" ) {
      out = LU_SOLVER;
    }
    else if ( in == "directLDL" ) {
      out = LDL_SOLVER;
    }
    else if ( in == "directLDL2" ) {
      out = LDL_SOLVER2;
    }
    else if ( in == "pardiso" ) {
      out = PARDISO;
    }
    else if ( in == "ilupack" ) {
      out = ILUPACK_SOLVER;
    }
    else {
      (*error) << "No enumeration value found in SolverType for '"
           << in << "'\n A missing case implementation?";
      Error( __FILE__, __LINE__ );
    }
  }

 // Specialisation for EigenSolverType
  template<> 
  void String2Enum<EigenSolverType>( const std::string &in, EigenSolverType &out ) {

    if ( in == "no eigensolver" ) {
      out = NOEIGENSOLVER;
    }
    else if ( in == "arpack" ) {
      out = ARPACK;
    }
    else if ( in == "subspace" ) {
      out = SUBSPACE;
    } 
    else {
      (*error) << "No enumeration value found in EigenSolverType for '"
               << in << "'\n A missing case implementation?";
      Error( __FILE__, __LINE__ );
    }
  }
  // Specialisation for PrecondType
  template<> 
  void String2Enum<PrecondType>( const std::string &in, PrecondType &out ) {

    if ( in == "noPrecond" ) {
      out = NOPRECOND;
    }
    else if ( in == "Id" ) {
      out = ID;
    }
    else if ( in == "MG" ) {
      out = MG;
    }
    else if ( in == "Jacobi" ) {
      out = JACOBI;
    }
    else if ( in == "SSOR" ) {
      out = SSOR;
    }
    else if ( in == "ILU0" ) {
      out = ILU0;
    }
    else if ( in == "ILUTP" ) {
      out = ILUTP;
    }
    else if ( in == "ILUK" ) {
      out = ILUK;
    }
    else if ( in == "ILDL0" ) {
      out = ILDL0;
    }
    else if ( in == "ILDLK" ) {
      out = ILDLK;
    }
    else if ( in == "ILDLTP" ) {
      out = ILDLTP;
    }
    else if ( in == "ILDLCN" ) {
      out = ILDLCN;
    }
    else if ( in == "IC0" ) {
      out = IC0;
    }
    else {
      (*error) << "No enumeration value found in PrecondType for '"
           << in << "'\n A missing case implementation?";
      Error( __FILE__, __LINE__ );
    }
  }

  // map specialisation for interpolation types
  template <>
  void String2Enum<AMGInterpolationType>( const std::string    &in,
                                          AMGInterpolationType &out )
  {
    if ( in == "constant" ) {
      out = AMG_INTERPOLATION_CONSTANT;
    } else if( in == "simpleWeighted" ) {
      out = AMG_INTERPOLATION_SIMPLE_WEIGHTED;
    } else if( in == "develop" ) {
      out = AMG_INTERPOLATION_DEVELOP;
    } else {
      (*error)
       << "No enumeration value found in AMGInterpolationType "
       << "for '" << in << "'\n A missing case implementation?";
      Error( __FILE__, __LINE__ );
    }
  }

  // map specialisation for interpolation types
  template <>
  void String2Enum<AMGSmootherType>( const std::string &in,
                                     AMGSmootherType   &out )
  {
    if ( in == "GaussSeidel" ) {
      out = AMG_SMOOTHER_GAUSSSEIDEL;
    } else if( in == "Jacobi" ) {
      out = AMG_SMOOTHER_DAMPED_JACOBI;
    } else {
      (*error)
       << "No enumeration value found in AMGSmoothertype "
       << "for '" << in << "'\n A missing case implementation?";
      Error( __FILE__, __LINE__ );
    }
  }

  // Specialisation for ReorderingType
  template<>
  void String2Enum<ReorderingType>( const std::string &in,
                    ReorderingType &out ) {

    if ( in == "noReordering" )
      out = NOREORDERING;
    else if ( in == "Sloan" )
      out = SLOAN;
    else if ( in == "Metis" )
      out = METIS;
    else if ( in == "minimumDegree" )
      out = MINIMUM_DEGREE;
    else if ( in == "nestedDissection" )
      out = NESTED_DISSECTION;
    else {
      (*error) << "String '" << in << "' cannot be converted to item of "
           << "'ReorderingType'!";
      Error(  __FILE__, __LINE__ );
    }
  }

  // Specialisation for FEMatrixType
  template<>
  void String2Enum<FEMatrixType>( const std::string &in,
                                  FEMatrixType &out ) {

    if ( in == "nomatrixtype" )
      out = NOTYPE;
    else if ( in == "system" )
      out = SYSTEM;
    else if ( in == "stiffness" )
      out = STIFFNESS;
    else if ( in == "damping" )
      out = DAMPING;
    else if ( in == "convection" )
      out = CONVECTION;
    else if ( in == "mass" )
      out = MASS;
    else {
      (*error) << "String '" << in << "' cannot be converted to item of "
           << "'FEMatrixType'!";
      Error(  __FILE__, __LINE__ );
    }
  }

  // Specialisation for IDBCType
  template<>
  void String2Enum<IDBCType>( const std::string &in, IDBCType &out ) {

    if ( in == "notype" )
      out = IDBC_NOTYPE;
    else if ( in == "elimination" )
      out = IDBC_ELIMINATION;
    else if ( in == "penalty" )
      out = IDBC_PENALTY;
    else {
      (*error) << "String '" << in << "' cannot be converted to item of "
           << "'IDBCType'!";
      Error(  __FILE__, __LINE__ );
    }
  }

  // Function that constitutes as sort of factory
  void EnumConversionFactory() {

    ENTER_FCN( "Function: EnumConversionFactory" );

    std::string auxString = "dummyString";

    // SolverType
    {
      SolverType auxType = NOSOLVER;
      auxString = Enum2String( auxType );
      String2Enum( auxString, auxType );
    }

    // PrecondType
    {
      PrecondType auxType = NOPRECOND;
      auxString = Enum2String( auxType );
      String2Enum( auxString, auxType );
    }

    // MatrixStorageType
    {
      MatrixStorageType auxType = NOSTORAGETYPE;
      auxString = Enum2String( auxType );
      String2Enum( auxString, auxType );
    }

    // MatrixEntryType
    {
      MatrixEntryType auxType = NOENTRYTYPE;
      auxString = Enum2String( auxType );
      // String2Enum( auxString, auxType );
    }

    // StopCritType
    {
      StopCritType auxType = NOSTOPCRITTYPE;
      auxString = Enum2String( auxType );
      String2Enum( auxString, auxType );
    }

    // ReorderingType
    {
      ReorderingType auxType = NOREORDERING;
      auxString = Enum2String( auxType );
      String2Enum( auxString, auxType );
    }

    // FEMatrixType
    {
      FEMatrixType auxType = NOTYPE;
      auxString = Enum2String( auxType );
      String2Enum( auxString, auxType );
    }

    // IDBCType
    {
      IDBCType auxType = IDBC_NOTYPE;
      auxString = Enum2String( auxType );
      String2Enum( auxString, auxType );
    }
  }

}
