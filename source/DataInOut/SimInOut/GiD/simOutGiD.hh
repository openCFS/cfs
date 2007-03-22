// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_OUTGID_HH
#define FILE_CFS_SIM_OUTGID_HH

#include <Domain/grid.hh>
#include <DataInOut/simOutput.hh>
#include <Domain/resultInfo.hh>

// Header of postprocessing library
#include <gidpost.h>

namespace CoupledField
{

  //! Class for writing results in the GiD postprocessing format
  class SimOutputGiD: virtual public SimOutput {
    
  public:

    //! Constructor
    SimOutputGiD( const std::string& fileName, ParamNode * outputNode );
  
    //! Destructor
    virtual ~SimOutputGiD();
  
    //! Initialize class
    void Init( Grid* grid );

    //! Write grid definition in file
    void WriteGrid();

    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd );
    
    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );

  private:

    //! Write nodal coordinates
    void WriteNodes();
    
    //! Write element declarations
    void WriteElements();
    
    //! Write transient results on nodes or elements
    void WriteNodeElemDataTrans( const Vector<Double> & var, 
                                 const StdVector<std::string> & dofNames,
                                 const std::string& name, 
                                 ResultInfo::EntryType entryType,
                                 ResultInfo::EntityUnknownType entityType,
                                 Double time );

    //! Write harmonic results on nodes or elements
    void WriteNodeElemDataHarm( const Vector<Complex>  & var, 
                                const StdVector<std::string> & dofNames,
                                std::string name, 
                                ResultInfo::EntryType entryType,
                                ResultInfo::EntityUnknownType entityType,
                                Double freq, 
                                ComplexFormat outputFormat );
    
    //! Write element to mesh file

    //! This method takes a pointer to an element, replaces it with 
    //! the degenerated element if necessary and writes this element
    //! to the mesh file.
    void WriteElement( Elem* ptEl );
    
    //! Map with result objects for each result type
    ResultMapType resultMap_;

    //! Dimension of grid
    GiD_Dimension dim_;

    //! Flag for binary file format
    bool isAscii_;

    //! Flag for use of degenerated elements
    bool degen3DElems_;

    //! Check if class is initialized
    bool isInitialized_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SimOutputGiD
  //! 
  //! \purpose 
  //! This class provides the interface for writing meshes and elements to the
  //! GiD own postprocessing files .post.msh and .post.res.
  //! This format is capable of visualizing volume elements, surface elements
  //! and named nodes. 
  //! /Note GiD has some restrictions regarding the visualization of element
  //! resuls on different element types. Therefore, if different types of 
  //! elements are contained in the grid (e.g. Quad/Tria, Tet/Hex/Pyra), all
  //! elements are treated as degenerated quadrilaterals or hexahedras. 
  //! This might cause some distortions in the resulting visualization, 
  //! especially if quadratic elements are used.
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
