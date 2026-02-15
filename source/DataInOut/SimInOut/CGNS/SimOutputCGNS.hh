// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_OUT_CGNS_HH
#define FILE_CFS_SIM_OUT_CGNS_HH

#include <Domain/Mesh/Grid.hh>
#include <Domain/Results/ResultInfo.hh>
#include <DataInOut/SimOutput.hh>

#include <cgnslib.h>

namespace fs = std::filesystem;

#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

namespace CoupledField
{

  //! Class for writing simulation result data (grid, results) in CGNS files.

  //! This class is for writing (C)FD (G)eneral (N)otation (S)ystem file(s).
  //!
  //! CGNS files can be opened in the following post-processors:
  //! ParaView: Off-the-shelf ParaView binaries from www.paraview.org may be used.
  //!           Since the VisItBridge inside ParaView does not support the
  //!           selection of subdomains according to CGNS sections, one zone
  //!           per region is written. Due to this, the ExtractBlock filter can
  //!           can be used to extract regions from openCFS by their name.
  //! CFDPost: The ANSYS CFD post processor can open only nodal transient results.
  //! EnSight: Has the best support next to ParaView.
  //! Tecplot: Only tested with Tecplot 360 2010 which needs CGNS files written
  //!          with CGNS library versions &lt;= 3.0.
  //! VisIt: Another free post-processor based on VTK. Up until version &lt;= 2.7
  //!        VisIt is linked against the CGNS library v3.0 which has some
  //!        incompatibilities due to a different definition of the element type
  //!        enum with version 3.1 (cf. cgnslib.h). Since openCFS is linked to
  //!        CGNS 3.1 there may occur problems!
  //!
  //! By default separate files for nodal and element results will be written.
  //! This controlled by the 'separateFile' attribute in the XML file.
  //! It has been found, that some post-processors for CGNS (VisIt &lt;= 2.7,
  //! EnSight &lt;= 10.0) have problems handling nodal and element results 
  //! defined under the same base node. ParaView (&lt;= 4.1) may go to 100% CPU 
  //! load and even crash! Even if the different result types are defined with
  //! respect to different base nodes in the same .cgns file, which would be okay
  //! for ParaView, VisIt and EnSight have problems handling this situation.
  //! Therefore, the default is to write node and element results to different
  //! files.
  //!
  //! As the name CGNS states, the format is mostly used for exchange of CFD data.
  //! Since most CFD codes are using linear cell types, most pre- (ANSYS WB) and
  //! post-processors have problems handling quadratic element types often used in
  //! FE codes, even if the CGNS standard properly defines such element types.
  //! Therefore, this class writes only linear represantions of quadratic elements.
  //! This behavior may be changed by the 'writeQuadElems' parameter.
  //!
  //! The CGNS standard defines two different so-called mid-level libraries, which
  //! may be used to actually write the binary data to files. These are HDF5 and ADF.
  //! Since older post-processors may have problems handling HDF5 .cgns files, we
  //! define ADF as standard. This behavior may be changed using the 'mll' attribute.
  //!
  //! \note ATTENTION: This class has been implemented for the exchange of CFD data
  //!                  during the E+H harmonic FSI project by Simon Triebenbacher.
  //!                  In this project we only need to transfer results for a single
  //!                  harmonic step. Therefore, no support for more time/freq steps
  //!                  has been implemented and only one multisequence step is supported!
  class SimOutputCGNS : virtual public SimOutput {
    
  public:
    
    //! Constructor
    SimOutputCGNS( std::string& fileName,
                   PtrParamNode outputNode,
                   PtrParamNode infoNode,
                   bool isRestart=false );
    
    //! Destructor
    virtual ~SimOutputCGNS();
  
    //! Initialize class
    void Init( Grid* grid, bool printGridOnly );

    //! Write grid definition in file
    void WriteGrid();

    void BeginMultiSequenceStep( UInt step,
                                 BasePDE::AnalysisType type,
                                 UInt numSteps  );

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
    
    //! Finish multisequence step
    void FinishMultiSequenceStep( );

    //! Finalize the output
    virtual void Finalize();

  private:

    void WriteNodesAndElements(UInt baseIdx);
    
    void WriteMixedSection(UInt baseIdx, const StdVector<Elem*>& elems);
    void WritePureSection(UInt baseIdx, const StdVector<Elem*>& elems);

    void TranslateConnectivity(Elem::FEType feType,
                               cgsize_t* cgnsConn,
                               Elem* elem);

    void InitElemTypeMap();
    
    std::map<Elem::FEType,CGNSLIB_H::ElementType_t> elemTypeMap_;

    StdVector<int> indexFile_;
    StdVector<int> indexBase_;
    std::map<UInt, StdVector<int> > idxZone_;
    StdVector<int> idxNodeSol_;
    StdVector<int> idxElemSol_;
    int cellDim_;
    bool outputFileOK_;
    bool writeAmplPhase_;
    bool separateFiles_;

    std::map< RegionIdType, std::map<UInt, UInt> > regionNodeMap_;

    //! Map with result objects for each result type
    ResultMapType resultMap_;
    
    //! Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
            
    //! Offset for step value in case of multisequence analysis
    Double stepValOffset_;

    bool writeQuadElems_;

    typedef std::map<std::string, std::map< RegionIdType, Vector<Double> > > RegionSolsType;
    

    //! for printing nodal results of simulation (static/transient)
    /*!
      \param definedOnNode is data defined on nodes?
      \param title title of the results.
      \param x array with nodal results
      \param step number of the step of the calculation
      \param time time of the calculation
    */
    void NodeElemDataTransient(const bool definedOnNode,
                               RegionSolsType& regionSols,
                               const UInt step, 
                               const Double time);
  
    void FillRegionSols(RegionSolsType& regionSols, 
                        const StdVector<shared_ptr<BaseResult> > & solList,
                        ResultInfo::EntityUnknownType entityType );

    //! Write parameter and xml file
    void WriteXmlFiles( fs::path simFile, fs::path matFile );

  };

} // end of namespace
 
#endif
