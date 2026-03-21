#ifndef FILE_MOD_ANSYSFILE_2002
#define FILE_MOD_ANSYSFILE_2002

#include <fstream>

#include <DataInOut/SimInput.hh>

namespace CoupledField 
{
  /** Class for reading in a mesh created by the ANSYS mkmesh-extension
   * or via our own create_mesh.py or mesh_tool.py.
   * 
   * It seems that his mesh format a native openCFS format and does not follow
   * any standard, yet it is quite nice and easiy to read and write format. */
  class SimInputMESH: public SimInput 
  {
  public:

    SimInputMESH(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode);
    
    virtual ~SimInputMESH();

    void InitModule() override;

    void ReadMesh(Grid *mi) override;
  
    unsigned int GetDim() override { return dim_; }

    /** Maps Ansys-element-ids to CFS-element-ids (BaseFE) */
    static Elem::FEType AnsysType2ElemType(unsigned int itype);

  private:  
    /** helper which processes num [xD Elements] with x is 1,2,3.
    * Old test meshes have the order 3D, 2D, 1D, hence we need to identify the case automatically :( */ 
    void ParseElementBlock();

    /** helper which parses a named entity block
    * @param label is [Node BC], [Save Nodes] or [Save Elements] */
    StdVector<std::pair<std::string, StdVector<unsigned int>>> ParseNamedEntities(const std::string& label, unsigned int num);

    //! Dimension of the mesh
    UInt dim_ = 0;

    //! Total number of elements
    UInt maxNumElems_ = 0;

    //! Total number of nodes
    UInt maxNumNodes_ = 0;
    
    //! Pointer to input file
    std::ifstream inFile_;

    /** the [Info] header, containing the first 9 key-value pairs
     * Version, Dimension, NumNodes, Num3DElements, Num2DElements, Num1DElements, 
     * NumNodeBC, NumSaveNodes,NumSavedElements  */
    std::map<std::string, unsigned int> header_;
  };
}

#endif
