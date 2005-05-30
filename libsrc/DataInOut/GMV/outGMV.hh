#ifndef FILE_OUTGMV_2002
#define FILE_OUTGMV_2002

#include "Domain/grid.hh"
#include "DataInOut/writeresults.hh"

namespace CoupledField
{

  //! This class provides an interface for writing files in gmv-format
  class WriteResultsGMV: virtual public WriteResults
  {
  public:

    //! Constructor
    WriteResultsGMV(const Char * const filename);
  
    //! Deconstructor
    virtual ~WriteResultsGMV();
  
    //! initialization with grid
    //! \param ptgrid pointer to grid object
    virtual void Init(Grid * aptgrid);
  
    //! write grid definition in file
    virtual void WriteGrid();


    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& sol, 
                                            const UInt step, 
                                            const Double time);
  
    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteElemSolutionTransient(const ElemStoreSol<Double>& data, 
                                            const UInt step, 
                                            const Double time);
  
    //! write element solution vector 
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step        step of calculation
      \param frequency   frequency of exciting function
      \param format      format for writing complex solution
      (real-imag/amplitude-phase)
    */
    virtual void WriteNodeSolutionHarmonic(const NodeStoreSol<Complex>& sol, 
                                           const UInt step,
                                           const Double frequency,
                                           const ComplexFormat format);
  
    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step      step of calculation
      \param frequency frequency of exciting function
      \param format    format for writing complex solution
      (real-imag/amplitude-phase)
    */
    virtual void WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& data, 
                                           const UInt step,
                                           const Double frequency,
                                           const ComplexFormat format);
  


    //! write comments
    /*!
      \param comments string with comments
    */
    virtual void WriteComments(const std::string comments){;}

    //! check, is it the gmv-output file
    virtual Boolean IsGMV(){ return TRUE;}

    //! function for open file with number num 
    void OpenFile(const Integer num);

  private:
    //! pointer to ofstream with history information
    std::ofstream * output;

    //! name of dir for output results
    Char * namedir_;

    // number of step
    UInt currStep_;

    //! number of last step
    UInt lastStep_;

    // current time 
    Double currTime_;

    // previous time 
    Double lastTime_;

    //! indicator of type for data
    Boolean ascii_;

    //! indicator of adaptive grid or not
    Boolean fixedgrid_; 

    //! True, if grid was already written one time
    Boolean firstGridWritten_;

    //! name of gridfile
    std::string nameGridFile_;

    //! write number of nodes and coordinates of them
    void WriteNodes();

    //! write cell description 
    void WriteCells(); 

    //! write cell materials
    void WriteMaterials();

    //! write variable information
    /*!
      \param dataType data type of the var: 0.. cell data, 1.. node data,
      2.. face data
      \param var      vector with data
      \param name     name of output-data
    */
    void WriteNodeVariableTransient(const Vector<Double> var, 
                                    const std::string name, 
                                    const UInt dataType);
  
    //! write variable information
    /*!
      \param dataType   data type of the var: 0.. cell data, 1.. node data,
      2.. face data
      \param var vector with data
      \param name name  of output-data
      \param outFormat  format of complex numbers
    */
    void WriteNodeVariableHarmonic(const Vector<Complex> var, 
                                   const std::string name, 
                                   const UInt dataType,
                                   const ComplexFormat outFormat);

    //! write vector-variable information
    /*!
      \param dataType data type of the var: 0.. cell data, 1.. node data,
      2.. face data
      \param var      pointer to vector with output data
      \param name     name of output-data
    */
    void WriteVelocity(const Vector<Double>* var, const std::string name,
                       const UInt dataType);
 
    //! transform string to 8 characters. we need it, because name in gmv,
    //! in binary format, should be from 8 characters
    /*!
      \param name (input) title
      \param result (output) result
    */
    void to8Char(const std::string name, char * result);
  
  };

} // end of namespace

#endif
