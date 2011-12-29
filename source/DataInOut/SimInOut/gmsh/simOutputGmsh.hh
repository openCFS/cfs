  // -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GMSHWRITER_2009
#define FILE_GMSHWRITER_2009

#include <iosfwd>
#include <map>
#include <string>
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/simOutput.hh"
#include "Domain/resultInfo.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "GmshHelper.hh"
#include "MatVec/vector.hh"
#include "PDE/basePDE.hh"
#include "Utils/StdVector.hh"
// #include <boost/bimap/set_of.hpp>
#include "boost/bimap.hpp"

namespace CoupledField {

  /**
   **/
class BaseResult;
class Grid;

  class SimOutputGmsh : public SimOutput
  {
  public:
    SimOutputGmsh(std::string fileName, PtrParamNode inputNode);
    virtual ~SimOutputGmsh();

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
    //! pointer to ofstream of current file
    std::ofstream * output_;

    //! current multiSequence step
    Integer currMsStep_;
    
    //! current analysis type
    BasePDE::AnalysisType currAnalysis_;

    //! only output grid
    bool printGridOnly_;

    //! indicator of type for data
    bool ascii_;

    //! indicator of endianess for binary data
    bool bigEndian_;

    //!  Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
        
    //! Offset for step value in case of multisequence analysis
    Double stepValOffset_;

    // Converters for endianess
    EndianSwapper<UInt> uiSwap_;
    EndianSwapper<Integer> iSwap_;
    EndianSwapper<Double> dSwap_;
    
    //! Map with result objects for each result type
    ResultMapType resultMap_;

    UInt numRegions_;

    //! open file for output
    std::ofstream*  OpenFile( const std::string& fileName );

    //! write number of nodes and coordinates of them
    void WriteNodes( std::ofstream * gridFile );

    //! write cell description 
    void WriteCells( std::ofstream * gridFile ); 

    //! write region description 
    void WriteRegions();
    
    void WriteNodeElemDataTrans( const Vector<Double> & var,
                                 const StdVector<std::string> & dofNames,
                                 const std::string& name,
                                 ResultInfo::EntryType entryType,
                                 ResultInfo::EntityUnknownType entityType,
                                 Double time );
    void WriteNodeElemDataHarm( const Vector<Complex> & var,
                                const StdVector<std::string> & dofNames,
                                const std::string name,
                                const ResultInfo::EntryType entryType,
                                ResultInfo::EntityUnknownType entityType,
                                const Double freq,
                                const ComplexFormat outputFormat );

  }; 

} 

#endif
