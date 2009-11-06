#ifndef CFS_HH_
#define CFS_HH_

#include "Utils/StdVector.hh"
#include "DataInOut/simInput.hh"
#include "DataInOut/DefineFiles/definefiles.hh"

namespace CoupledField
{
  class ResultHandler;
  class Timer;
  class MaterialHandler;

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

    /** Initialized global Enum<> objectes.
     * To be called also from cfstool.cc */
    static void SetGlobalEnums();

    /** delete the "global" objects */
    ~CFS();
  private:

    /** This is the internal main run.
     * To be split for more overview
     * @exception in an error case exceptions are forwarded */
    void Process();

    /** Write the skeleton file only */
    void WriteXMLSkeleton();

    /** Read XML file */
    void ReadXMLFile();

    /** Setup I/O without info.xml and problem.xml */
    void SetupIO();

    /** Write the grid to result but do not calculate */
    void PrintGrid();

    /** Solve the simulation/ optimition problem */
    void SolveProblem();

    DefineInOutFiles fileHandler;

    ResultHandler* resultHandler;

    /** Timer for the whole CFS runtime */
    Timer* timer;

    /** this is our hostname. Empty if it cannot be determined */
    std::string hostname_;

    /** This is a string for output with the start time */
    std::string start_time_;

    /** The object itself in kept in the fileHandler */
    MaterialHandler* materialHandler;

    std::map<std::string, StdVector<shared_ptr<SimInput> > > gridInputs;
  };
}
#endif /* CFS_HH_ */
