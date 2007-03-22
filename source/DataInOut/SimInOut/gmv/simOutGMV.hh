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
    SimOutputGMV(const std::string fileName, ParamNode * outputNode);
  
    //! Deconstructor
    virtual ~SimOutputGMV();
  
    //! write grid definition in file
    virtual void WriteGrid();

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step, AnalysisType type );
    
    //! Register result (within one multisequence step)
    virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                 UInt saveBegin, UInt saveInc,
                                 UInt saveEnd );

    //! Begin single analysis step
    virtual void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    virtual void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    virtual void FinishStep( );

    //! End multisequence step
    virtual void FinishMultiSequenceStep( );

    //! Finalize the output
    virtual void Finalize();

    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */

    //virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& sol, 
    //                                        const UInt step, 
    //                                        const Double time);
  
    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    //virtual void WriteElemSolutionTransient(const ElemStoreSol<Double>& data, 
    //                                        const UInt step, 
    //                                        const Double time);
  
    //! write element solution vector 
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step        step of calculation
      \param frequency   frequency of exciting function
      \param format      format for writing complex solution
      (real-imag/amplitude-phase)
    */
    //virtual void WriteNodeSolutionHarmonic(const NodeStoreSol<Complex>& sol, 
    //                                       const UInt step,
    //                                       const Double frequency,
    //                                       const ComplexFormat format);
  
    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step      step of calculation
      \param frequency frequency of exciting function
      \param format    format for writing complex solution
      (real-imag/amplitude-phase)
    */
    //virtual void WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& data, 
    //                                       const UInt step,
    //                                       const Double frequency,
    //                                       const ComplexFormat format);
  


    //! write comments
    /*!
      \param comments string with comments
    */
    virtual void WriteComments(const std::string comments){;}

    //! function for open file with number num 
    void OpenFile(const Integer num);

  private:
    //! pointer to ofstream with history information
    std::ofstream * output;

    // number of step
    UInt currStep_;

    //! number of last step
    UInt lastStep_;

    // current time 
    Double currTime_;

    // previous time 
    Double lastTime_;

    //! indicator of type for data
    bool ascii_;

    //! size of character string to be written out
    UInt charOutSize_;

    //! Map with result objects for each result type
    ResultMapType resultMap_;

    //! indicator of adaptive grid or not
    bool fixedgrid_; 

    //! True, if grid was already written one time
    bool firstGridWritten_;

    bool printGridOnly_;

    //! True, if unit of result should be printed to result name
    bool printUnit_;

    //! name of gridfile
    std::string nameGridFile_;


    //! write number of nodes and coordinates of them
    void WriteNodes();

    //! write cell description 
    void WriteCells(); 

    //! write named entities (nodes, elements)
    void WriteNamedEntities();

    void ElemType2GMVElemId(FEType et, std::string & id);

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

    //! write variable information
    /*!
      \param dataType data type of the var: 0.. cell data, 1.. node data,
      2.. face data
      \param var      vector with data
      \param name     name of output-data
    */
    //void WriteNodeVariableTransient(const Vector<Double> & var, 
    //                                const std::string name, 
    //                                const UInt dataType);
  
    //! write variable information
    /*!
      \param dataType   data type of the var: 0.. cell data, 1.. node data,
      2.. face data
      \param var vector with data
      \param name name  of output-data
      \param outFormat  format of complex numbers
    */
    //void WriteNodeVariableHarmonic(const Vector<Complex> & var, 
    //                               const std::string name, 
    //                               const UInt dataType,
    //                               const ComplexFormat outFormat);

    //! write vector-variable information
    /*!
      \param dataType data type of the var: 0.. cell data, 1.. node data,
      2.. face data
      \param var      pointer to vector with output data
      \param name     name of output-data
    */
    //void WriteVelocity(const Vector<Double>* var, const std::string name,
    //                   const UInt dataType);
 
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
