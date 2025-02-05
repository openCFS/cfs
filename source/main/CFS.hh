#ifndef CFS_HH_
#define CFS_HH_

#include "Utils/StdVector.hh"
#include "DataInOut/SimInput.hh"
#include "DataInOut/DefineInOutFiles.hh"

namespace CoupledField
{
  class ResultHandler;
  class Timer;
  class MaterialHandler;
  class LogConfigurator;
  class SimState;

  /** This is the base class of CFS.
   * It basically gives more structure for the original main() call. */
  class CFS
  {
  public:
    /** The constructor initialized the most important "global" objects
     * based on the command line parameters */
    CFS(int argc, const char **argv);

    /** Call this after the constructor. That is all.
     * @return the return value for main*/
    int Run();

    /** Initialized global Enum<> objects.
     * To be called also from cfstool.cc */
    static void SetGlobalEnums();

    /** Print exception, also to info.xml.
     * Export to be used by serial and mpi main */
    static void HandleException(const std::exception& ex);

    /** delete the "global" objects */
    ~CFS();
  private:

    /** Read XML file */
    void ReadXMLFile();

    /** Setup I/O without info.xml and problem.xml */
    void SetupIO(PtrParamNode rootNode );

    /** Write the grid to result but do not calculate */
    void PrintGrid();

    /** Solve the simulation/ optimition problem */
    void SolveProblem();

    //! Object for file I/O creation
    DefineInOutFiles fileHandler;

    //! Global result handler
    ResultHandler* resultHandler;
    
    //! Object for saving the simulation state
    shared_ptr<SimState> simState;

    //! Root node of parameter xml file
    PtrParamNode paramNode_;
    
    /** Timer for the whole CFS runtime */
    shared_ptr<Timer> timer;
    
    /** this is our hostname. Empty if it cannot be determined */
    std::string hostname_;

    /** This is a string for output with the start time */
    std::string start_time_;

    /** The object itself in kept in the fileHandler */
    shared_ptr<MaterialHandler> materialHandler;
    
    std::map<std::string, StdVector<shared_ptr<SimInput> > > gridInputs;
  };
}
#endif /* CFS_HH_ */
