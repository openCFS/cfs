// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_OUT_CGNS_HH
#define FILE_CFS_SIM_OUT_CGNS_HH

#include <Domain/Mesh/Grid.hh>
#include <Domain/Results/ResultInfo.hh>
#include <DataInOut/SimOutput.hh>

#ifdef UNV_NONSTD_WIDTH
#define UNV_WIDTH 14
#else
#define UNV_WIDTH 13
#endif

namespace CoupledField
{

  //! Class for writing output-information(grid,results) in unverg-format file
  
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

  private:

    void writeNodesAndElements();
    
    int indexFile_, indexBase_, indexZone_;
    bool outputFileOK_;
    char baseName_[33], zoneName_[33];

    //! pointer to the ofstream with output -data
    std::ofstream * output;

    //! Map with result objects for each result type
    ResultMapType resultMap_;
    
    //! Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
            
    //! Offset for step value in case of multisequence analysis
    Double stepValOffset_;

    //! for printing nodal results of simulation (static/transient)
    /*!
      \param dataSetNr number of dataset (55/56)
      \param title title of the results.
      \param x array with nodal results
      \param step number of the step of the calculation
      \param time time of the calculation
    */
    void NodeElemDataTransient(const UInt dataSetNr,
                               const std::string & title, 
                               const Vector<Double> & x, 
                               const UInt step, 
                               const Double time, 
                               const UInt nrNodes,
                               const UInt nrDofs=1);
  
    //! for printing nodal results of simulation (harmonic)
    /*!
      \param dataSetNr number of dataset (55/56)
      \param title title of the results.
      \param x array with nodal results
      \param freuqncy exciting frequency of current result
      \param format output format for complex numbers
    */
    void NodeElemDataHarmonic(const UInt dataSetNr,
                              const std::string & title, 
                              const Vector<Complex> & x, 
                              const UInt step,
                              const Double frequency,
                              const ComplexFormat format,
                              const UInt nrNodes,
                              const UInt nrDofs=1);
  
    //! Convertes enum SolutionType to string
    std::string SolutionTypeToString(const SolutionType type) const;
  
    //! Re-sort stresses according to vector with dofnames
    template<class TYPE>
    void SortStresses( Vector<TYPE>& vec,
                       const StdVector<std::string>& dofNames );
  };

} // end of namespace
 
#endif
