#ifndef FILE_CFSOLASPARAMS
#define FILE_CFSOLASPARAMS

namespace CoupledField
{

  // Forward declaration of classes
  class BaseParamHandler;
  class OLAS_Params;

  //! Class for passing steering parameters from CFS++ to OLAS
  class CFSOLASParams {

  public:
    //! Document me!!!
    static void SetParams( std::string pdename, BaseParamHandler *cfs,
                           OLAS_Params *olas, bool overrideExpert = false );

  private:

    //! Set parameters for solver for linear system

    //! This routine queries the BaseParamHandler object for the parameters
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
    static void SetSolverParams( std::string pdename, BaseParamHandler *cfs,
                                 OLAS_Params *olas, SolverType sType );

    //! Set parameters for preconditioner for linear system

    //! This routine queries the BaseParamHandler object for the parameters
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
    static void SetPrecondParams( std::string pdename, BaseParamHandler *cfs,
                                  OLAS_Params *olas, PrecondType pType );

    //! Expert routine for correcting parameter inconsistencies
    static void Expert( BaseParamHandler *cfs,
                        std::string pdename, SolverType &sType,
                        PrecondType &pType, MatrixStorageType &mType,
                        MatrixEntryType &eType, ReorderingType &rType );

  };

}

#endif
