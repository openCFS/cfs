#ifdef USE_OLAS

#ifndef FILE_CFSOLASPARAMS
#define FILE_CFSOLASPARAMS

#ifdef OLAS
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
    static void SetSolverParams( std::string pdename, BaseParamHandler *cfs,
				 OLAS_Params *olas, SolverType sType );

    //! Set parameters for preconditioner for linear system

    //! This routine queries the BaseParamHandler object for the parameters
    //! belonging to the preconditionerr for the linear system associated to
    //! the PDE specified via the pdename input parameter and inserts them into
    //! The olasParams object.
    static void SetPrecondParams( std::string pdename, BaseParamHandler *cfs,
				  OLAS_Params *olas, PrecondType pType );

    //! Expert routine for correcting parameter inconsistencies
    static void Expert( SolverType &sType, PrecondType &pType,
			MatrixStorageType &mType, MatrixEntryType &eType );

  };

}

#endif

#endif
