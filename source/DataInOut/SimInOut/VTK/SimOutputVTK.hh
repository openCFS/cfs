/*
 * SimOutputVTK.hh
 *
 *  Created on: Nov 17, 2008
 *      Author: fwein
 */

#ifndef SIMOUTPUTVTK_HH_
#define SIMOUTPUTVTK_HH_

#include <string>
#include <map>
#include <vector>

#include <DataInOut/simOutput.hh>
#include <Utils/result.hh>
#include <Domain/grid.hh>
#include <Domain/resultInfo.hh>


class vtkPointSet;
class vtkXMLDataElement;
class vtkUnstructuredGrid;
class vtkDoubleArray;
class vtkCellType;
class vtkXMLWriter;

namespace CoupledField
{
  //! This class provides an interface for writing files in VTK-format
  class SimOutputVTK : public SimOutput
  {
  public:

    //! Constructor
    SimOutputVTK(const std::string& fileName, ParamNode* outputNode);
  
    //! Deconstructor
    virtual ~SimOutputVTK();
  
    //! Initialize class
    void Init( Grid* ptGrid, bool printGridOnly );

    //! write grid definition in file
    void WriteGrid( bool printGridOnly );
    
    //! Begin multisequence step
    void BeginMultiSequenceStep( UInt step,
                                 BasePDE::AnalysisType type,
                                 UInt numSteps );
    
    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                          UInt saveBegin, UInt saveInc,
                          UInt saveEnd,
                          bool isHistory );
    
    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );
    
    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );
    
    //! End single analysis step
    void FinishStep( );
    
    //! End multisequence step
    void FinishMultiSequenceStep( );
    
    //! Finalize the output
    void Finalize();

 
  private:

    class RegionData
    {
    public:
      RegionData()
      {
        regionId = NO_REGION_ID;
      }
      
      /** This is our unstructured grid representing a region */
      vtkUnstructuredGrid* grid; 

      /** this are all solutions for our region */
      StdVector<shared_ptr<BaseResult> > results;
      
      /** our region id == index */
      RegionIdType regionId;
    };
    
    
    /** Create a unstructured grid for a region.
     * Adds all elements/cells and nodes/points. Makes refeneces in the
     * node_point_mapping and element_cell_mapping for the specified regionId.
     * Within a unstructured grid we always start from 0 */
    void CreateMesh(RegionData& data);

    
    /** Writes a specific result for a single region. Harmonic stuff is splitted
     * into imag/real or ampl/phase producing a further VTK result */
    void WriteData(RegionData& data, shared_ptr<BaseResult> sol, double actStepVal, bool harmonic);

    /** Writes the current region data to a single file. We add ourself to singleFiles */
    void WriteSingleXMlFile();
    
    /** Helper for WriteSingleXMlFile() which creates the proper xml writer and sets the input.
     * @return the user has to delete it! */
    vtkXMLWriter* CreateXMlWriter();
    
    /** Helper for WriteData. For complex data we write an on vtk result for the real/ampl part and one
     * for the imag/phase part. Or the complex data is interpreted as real (@see IgnoreImaginaryPart) */
    template<class TYPE>
    void WriteSingleVTKData(RegionData& data, shared_ptr<BaseResult> sol, double actStepVal, std::string& label, bool harmonic, bool real);

    
    /** Helper method for WriteData. Adjusts the dimension (e.g. displacement is always xyz
     * such that Paraview can show the deformation (warp filter)
     * @param vtk_out allocates space and fills with values data
     * @param harmonic does values contain complex numbers?
     * @param real only the real part? Depending of output complex format this is the amplitude
     *        if false the for complex data the imaginary/phase information is extracted from values */
    template<class TYPE>
    void TransferData(ResultInfo& info, vtkDoubleArray* vtk_out, CFSVector* values, bool harmonic, bool real);
    
    /** For some results there is never an imaginary part (e.g. optimization design variables)
     * TODO: better create real results in this case and use also for other writers - Fabian
     * @return false does not mean that we are complex.*/
    static bool IgnoreImaginaryPart(ResultInfo& info);
    
    /** translate from CFS element type to VTK cell type.
     * @throws exception for incompatible types.
     * @return vtkCellType is a enum - we just don't want to include the vtk stuff here */
    int GetCellIdType(FEType cfsType); 
    
    
    /** do we write single files for the steps? Then we also write into a subdirectory and a pvd file is created */
    bool writeSingleFiles_;
    
    struct SingleVTUFile
    {
       std::string name;
       int         id;        // the glocal counter -> index within singleFiles  
       int         an_step;   // the current step within analyis
       int         ms_step;   // current multisequence step
       double      num_label; // iteration, frequency or timestep
    };
    
    /** in case of single files, this are the written filenames */
    StdVector<SingleVTUFile> singleFiles;
    
    //! current multiSequence step
    Integer currMsStep_;
    
    //! current analysis type
    BasePDE::AnalysisType currAnalysis_;

    //!  Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
        
    //! Offset for step value in case of multisequence analysis
    Double stepValOffset_;
    
    // fast access to the regions
    StdVector<RegionIdType> regions;
    
    /** An unstructured grid for each region. The index is the regionId - hence data_ might be larger than regions. */
    StdVector<RegionData> data_; 
  };

} // end of namespace


#endif /* SIMOUTPUTVTK_HH_ */
