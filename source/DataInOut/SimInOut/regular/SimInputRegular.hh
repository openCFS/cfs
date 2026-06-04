#ifndef SIM_INPUT_REGULAR_HH
#define SIM_INPUT_REGULAR_HH

#include <set>

#include "DataInOut/SimInput.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
  class GridCFS;
  
  /** class that generates a mesh inside cfs for regular rectangular or cubic domains */
  class SimInputRegular : virtual public SimInput 
  {
  public:
    /** @param fileName is rarely used and optional. 
        @param pn is null if loaded as an external file via -m */
    SimInputRegular(string fileName, PtrParamNode pn, PtrParamNode info);
    
    virtual ~SimInputRegular() {}

    /** we need to do the major initialization late with ReadMesh() */
    void InitModule() override {};

    /** This is the actual constructor as we need to be called after domain initialization.
     @see ParseXML() */
    void ReadMesh(Grid *grid) override;
  
    //! Return dimension of the mesh
    unsigned int GetDim() override { return dim_; }
   
  private:    

    /** input is either from 'regular' in problem.xml or from the external file.
     * Needs to be called later and after domain initialization as the domain mathParser is used. */
    void ParseXML(const PtrParamNode& pn);

    /** process region element, is like reading named elems but replacing the region for them. 
     * Easy before Grid::FinishInit() is called */
    void SetRegion(GridCFS* grid, const PtrParamNode& pn);

    void Set2DNamedNodes(Grid* grid);

    void Set3DNamedNodes(Grid* grid);

    /** add a named node set as a Python-like range start..end (node numbers are made
     * 1-based internally). With only 'start' given a single node is added; default step 1 */
    void AddNamedNodes(Grid* grid, const std::string& name,
                       unsigned int start, unsigned int end = 0, unsigned int step = 1);

    //! Dimension of the mesh
    unsigned int dim_ = 0;

    /** this is the number of nodes for the current grid, before we add stuff.
     * when we have only one mesh per grid (the common case), this is zero */
    unsigned int base_ = 0;

    //! Number of elements in x, y and z-direction
    Vector<int> n_;

    /** corner with minimal coordinates. Max = min_ + n_ * h_ */
    Vector<double> min_;

    /** spacing, see min_ */
    Vector<double> h_;

    /** this is the xml description of the regular mesh. Either from problem.xml or file */ 
    PtrParamNode xml_;
    
    /** is the name of the single region */
    std::string region_;

    /** Info Node base */
    PtrParamNode info_;
  };

}

#endif // SIM_INPUT_REGULAR_HH
