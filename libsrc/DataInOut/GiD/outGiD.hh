#ifndef FILE_OUTGID_HH
#define FILE_OUTGID_HH

#include "Domain/grid.hh"
#include "DataInOut/writeresults.hh"

// Header of postprocessing library
#include "gidpost.h"

namespace CoupledField
{

  //! Class for writing results in the GiD postprocessing format
  class WriteResultsGiD: virtual public WriteResults
  {
  public:

    //! Constructor
    WriteResultsGiD( const Char * const filename );
  
    //! Deconstructor
    virtual ~WriteResultsGiD();
  
    //! Initialization with grid
    //! \param ptgrid pointer to grid object
    virtual void Init( Grid * aptgrid );
  
    //! Write grid definition in file
    virtual void WriteGrid();


    //! Write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteNodeSolutionTransient( const NodeStoreSol<Double>& sol, 
                                             const UInt step, 
                                             const Double time );
  
    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteElemSolutionTransient( const ElemStoreSol<Double>& data, 
                                             const UInt step, 
                                             const Double time );
  
    //! Write element solution vector 
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step        step of calculation
      \param frequency   frequency of exciting function
      \param format      format for writing complex solution
      (real-imag/amplitude-phase)
    */
    virtual void WriteNodeSolutionHarmonic( const NodeStoreSol<Complex>& sol, 
                                            const UInt step,
                                            const Double frequency,
                                            const ComplexFormat format );
  
    //! Write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step      step of calculation
      \param frequency frequency of exciting function
      \param format    format for writing complex solution
      (real-imag/amplitude-phase)
    */
    virtual void WriteElemSolutionHarmonic( const ElemStoreSol<Complex>& data, 
                                            const UInt step,
                                            const Double frequency,
                                            const ComplexFormat format );

  private:

    //! Write nodal coordinates
    void WriteNodes();
    
    //! Write element declarations
    void WriteElements();

    //! Write transient results on nodes or elements
    void WriteNodeElemDataTrans( const Vector<Double> & var, 
                                 const UInt dof,
                                 const std::string name, 
                                 const GiD_ResultLocation entType, 
                                 const Double time );

    //! Write harmonic results on nodes or elements
    void WriteNodeElemDataHarm( const Vector<Complex>  & var, 
                                const UInt dof,
                                const std::string name, 
                                const GiD_ResultLocation entType, 
                                const Double freq, 
                                const ComplexFormat outputFormat );
    
    //! Write element to mesh file

    //! This method takes a pointer to an element, replaces it with 
    //! the degenerated element if necessary and writes this element
    //! to the mesh file.
    void WriteElement( Elem* ptEl );
    
    //! Dimension of grid
    GiD_Dimension dim_;

    //! Flag for binary file format
    Boolean isAscii_;

    //! Flag for use of degenerated elements
    Boolean degen3DElems_;
  };

#endif#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class WriteResultsGiD
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
  //! /Note When writing binary files, GiD shows all regions in the same 
  //! color. Therefore, by default only Ascii format is used.
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
  //! - Add support for tensor results
  //! - Add support for named elements
  //! 

#endif

} // end of namespace

#endif
