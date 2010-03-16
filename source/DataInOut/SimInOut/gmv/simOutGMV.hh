// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GMVWRITER_2002
#define FILE_GMVWRITER_2002

#include <string>

#include <DataInOut/simOutput.hh>
#include <Domain/grid.hh>
#include <Domain/resultInfo.hh>

namespace CoupledField
{
  //! This class provides an interface for writing files in gmv-format
  class SimOutputGMV : public SimOutput
  {
  public:

    //! Constructor
    SimOutputGMV(const std::string fileName, PtrParamNode outputNode);
  
    //! Deconstructor
    virtual ~SimOutputGMV();
  
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

    //! open file for output
    std::ofstream*  OpenFile( const std::string& fileName );

    //! print end statement to file
    void PrintFileEpilog( std::ofstream * outFile,
                          bool printStepInfo );

    //! pointer to ofstream of current file
    std::ofstream * output;

    //! current multiSequence step
    Integer currMsStep_;
    
    //! current analysis type
    BasePDE::AnalysisType currAnalysis_;

    //! indicator of type for data
    bool ascii_;

    //!  Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
        
    //! Offset for step value in case of multisequence analysis
    Double stepValOffset_;
    
    //! size of character string to be written out
    UInt charOutSize_;

    //! Map with result objects for each result type
    ResultMapType resultMap_;

    //! indicator of adaptive grid or not
    bool fixedgrid_; 

    //! True, if grid was already written one time
    bool firstGridWritten_;

    //! True, if unit of result should be printed to result name
    bool printUnit_;

    //! write number of nodes and coordinates of them
    void WriteNodes( std::ofstream * gridFile );

    //! write cell description 
    void WriteCells( std::ofstream * gridFile ); 

    //! write named entities (nodes, elements)
    void WriteNamedEntities( std::ofstream * gridFile );

    void ElemType2GMVElemId(Elem::FEType et, std::string & id);

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

    void WriteVariable( const std::string name,
                        const UInt dataType,
                        const Vector<Double> & var);

    void WriteVector( const std::string name,
                      const UInt dataType,
                      const Vector<Double> & var,
                      const StdVector<std::string> & dofNames);

  
 
    //! Truncate a std::string and copy the result into a *char array.
    //! This function is needed, since gmv permits only to write 8 characters
    //! (binary format up to gmv 3.4) or 32 characters (ascii format, binary 
    //! format since version 3.5)
    /*!
      \param name (input) title
      \param suffix (input) suffix for name
      \param result (output) result
    */
    void TruncateString(std::string & str,
                        UInt maxLen);
  
  };

} // end of namespace

#endif
