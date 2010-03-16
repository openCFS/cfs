// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFSOLASPARAMS
#define FILE_CFSOLASPARAMS

#include "OLAS/solver/basesolver.hh"
#include "OLAS/solver/baseEigensolver.hh"
#include "OLAS/precond/baseprecond.hh"
#include "OLAS/graph/baseordering.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

  // Forward declaration of classes
  class Assemble;

  //! Class for passing steering parameters from CFS++ to OLAS

  //! This class constitutes the interface between CFS++ and OLAS with respect
  //! to setting parameters for the solution of the linear system. For clarity
  //! it offers a single public method which will handle the setting of
  //! parameters in the OLAS parameter object depending on the parameters
  //! obtained from CFS++.
  class CFSOLASParams {

  public:

    //! Central method for passing steering parameters to OLAS

    //! This method can be used to pass steering parameters from CFS++ to OLAS.
    //! The method queries the parameter handler for the parameters related to
    //! the matrix, solver and preconditioner types, converts them to the
    //! respective OLAS settings and inserts them into the OLAS parameter
    //! object. By default the Expert() module is called to check the
    //! consisteny and sensibility of the parameters and adapt them, if
    //! necessary. This behaviour can be disabled by setting overrideExpert
    //! to no.
    //! \param pdename        symbolic name of PDE; the value is matched
    //!                       against the 'name' attribute of the 'system'
    //!                       element in the 'linearSystems' section of the
    //!                       XML file
    //! \param cfs            pointer to %CFS parameter object
    //! \param olas           pointer to %OLAS parameter object
    //! \param analysistype   specifies the type of analysis (required e.g.
    //!                       for correcting the matrix entry type)
    //! \param overrideExpert if set to 'false' the method will call the
    //!                       Expert() module in order to validate or
    //!                       improve the solution parameters or set those
    //!                       no specified by the user
    static void SetParams( std::string pdename, PtrParamNode cfs,
                           BasePDE::AnalysisType analysistype,
                           Assemble * assemble,
                           bool overrideExpert = false );

  private:

    //! Set parameters for solver for linear system

    //! This routine queries the ParamNode object for the parameters
    //! belonging to the solver for the linear system associated to the
    //! PDE specified via the pdename input parameter and inserts them into
    //! The olasParams object.
    //! \note The solver parameters in the XML file are optional. This routine
    //!       relies on GetList() to determine, whether a parameter is set.
    //!       If there is exactly one occurence, i.e. GetList() returns a
    //!       vector of length 1, then the parameter is set to the value found
    //!       in the XML file. If there is no match or a multiple match, then
    //!       the parameter is not set. The case of multiple matches can only
    //!       occur, if a non-validating parser is used, since the Schema
    //!       definitions require the parameters to be unique.
    static void SetSolverParams( std::string pdename, PtrParamNode cfs,
                                 BaseSolver::SolverType sType );

    //! Set parameters for eigenvalue solver for linear system

    //! This routine queries the ParamNode object for the parameters
    //! belonging to the eigenvalue solver for the linear system associated to
    //! the PDE specified via the pdename input parameter and inserts them 
    //! into rhe olasParams object.
    //! \note The solver parameters in the XML file are optional. This routine
    //!       relies on GetList() to determine, whether a parameter is set.
    //!       If there is exactly one occurence, i.e. GetList() returns a
    //!       vector of length 1, then the parameter is set to the value found
    //!       in the XML file. If there is no match or a multiple match, then
    //!       the parameter is not set. The case of multiple matches can only
    //!       occur, if a non-validating parser is used, since the Schema
    //!       definitions require the parameters to be unique.
    static void SetEigenSolverParams( std::string pdename, 
                                      PtrParamNode cfs,
                                      BaseEigenSolver::EigenSolverType sType );    

    //! Set parameters for preconditioner for linear system

    //! This routine queries the ParamNode object for the parameters
    //! belonging to the preconditionerr for the linear system associated to
    //! the PDE specified via the pdename input parameter and inserts them into
    //! The olasParams object.
    //! \note The solver parameters in the XML file are optional. This routine
    //!       relies on GetList() to determine, whether a parameter is set.
    //!       If there is exactly one occurence, i.e. GetList() returns a
    //!       vector of length 1, then the parameter is set to the value found
    //!       in the XML file. If there is no match or a multiple match, then
    //!       the parameter is not set. The case of multiple matches can only
    //!       occur, if a non-validating parser is used, since the Schema
    //!       definitions require the parameters to be unique.
    static void SetPrecondParams( std::string pdename, PtrParamNode cfs,
                                  BasePrecond::PrecondType pType );

    //! Expert routine for correcting parameter inconsistencies
    static void Expert( PtrParamNode cfs,
                        std::string pdename, 
                        BaseEigenSolver::EigenSolverType &esType,
                        BaseSolver::SolverType &sType,
                        BasePrecond::PrecondType &pType, 
                        BaseMatrix::StorageType &mType,
                        BaseMatrix::EntryType &eType, 
                        BaseOrdering::ReorderingType &rType,
                        BasePDE::AnalysisType analysisType,
                        Assemble * assemble,
                        bool allowChangeOfReordering );

  };

}

#endif
