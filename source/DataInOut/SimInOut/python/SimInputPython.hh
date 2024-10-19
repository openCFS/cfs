#ifndef SIM_INPUT_PYTHON_HH
#define SIM_INPUT_PYTHON_HH

#include <DataInOut/SimInput.hh>
#include <Utils/PythonKernel.hh>

namespace CoupledField
{
  /** execute python code (e.g. by loading a file) which sets the mesh.
   * USE_EMBEDDED_PYTHON is required as cfs starts it's own interpreter, offering a cfs package
   * to be imported where cfs functions wait to be called. */
  class SimInputPython: virtual public SimInput
  {
    /** the PythonKernel may call our private SetNodes(), SetRegions() from python */
    friend class PythonKernel;

  public:
    /** fileName might be empty, then we get the stuff from inputNode */
    SimInputPython(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode);

    virtual ~SimInputPython();

    void InitModule() override;

    /** expects a python function set_cfs_mesh() which is called with a dict of options */
    virtual void ReadMesh(Grid* grid) override;

    virtual UInt GetDim() override;

  private:

    /** to be called first as cfs.set_nodes with a numpy array with three columns of doubles with the nodes coordinates.
     * The nodes are referred to as 1-based ids.
     * Called from PythonKernel, which is our friend (can call private functions) */
    void SetNodes(pyObject* args);

    /** to be called second as cfs.set_regions with a list of region names. The region id for the other functions is the 0-based index */
    void SetRegions(pyObject* args);

    /** to be called after SetRegions() as cfs.add_elements with the following parameters. Can be called many times for different fe types
     * - integer of the total number of elements
     * - integer with the Elem::FEType
     * - numpy array with the columns element number (1-based), region id, node numbers (1-based) corresponding to fe type */
    void AddElements(pyObject* args);

    /** to be called optionally after SetNodes as cfs.add_named_nodes , can be called multiple times.
     * first parameter is a string with a name, second is a numpy array of dtype=numpy.uintc with 1-based node ids */
    void AddNamedNodes(pyObject* args);

    /** to be called optionally after the last AddElements() call as cfs.add_named_elements, can be called multiple times.
     * first parameter is a string with a name, second is a numpy array of dtype=numpy.uintc with 1-based element ids */
    void AddNamedElements(pyObject* args);

    /** helper for AddNamedNodes() and AddNamedElements() */
    void AddNamedNodesElements(pyObject* args, bool nodes);

    /** all AddElements() must match and sum up with this number. */
    unsigned int total_elements_ = 0;

    /** given filename. Possibly via file and path which path = cfs:share:python key */
    std::string givenname;

    /** the python options from the xml file */
    StdVector<std::pair<std::string, std::string> > options;

    /** this represents the python script */
    pyObject* module = NULL;

    /** set by ReadMesh() */
    Grid* grid = NULL;

    /** maximal dimension obtained by AddElements() */
    unsigned int dim_ = 0;

    /** child element of base class myInfo_ in header/domain */
    PtrParamNode info_;
  }; // end class

} // end name space

#endif // SIM_INPUT_PYTHON_HH
