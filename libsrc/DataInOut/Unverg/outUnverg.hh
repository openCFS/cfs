#ifndef FILE_OUTRESULTUNVERG_2001
#define FILE_OUTRESULTUNVERG_2001

#include "Domain/grid.hh"
#include "DataInOut/writeresults.hh"

namespace CoupledField
{

  //! Class for writing output-information(grid,results) in unverg-format file

  class WriteResultsUnverg: virtual public WriteResults
  {

  public:
    //! constructor with name of a file for results
    WriteResultsUnverg(const Char * const filename);

    //! deconstructor
    virtual ~WriteResultsUnverg();

    //! initialization with grid
    //! \param ptgrid pointer to grid object
    virtual void Init(Grid * aptgrid);

    //! write information about grid  in file
    virtual void WriteGrid();

    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& data, 
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
      \param step step of calculation
      \param frequency frequency of exciting function
      \param format format for writing complex solution (real-imag/amplitude-phase)
    */
    virtual void WriteNodeSolutionHarmonic(const NodeStoreSol<Complex>& data, 
                                           const UInt step,
                                           const Double frequency,
                                           const ComplexFormat format);

    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param frequency frequency of exciting function
      \param format format for writing complex solution (real-imag/amplitude-phase)
    */
    virtual void WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& data, 
                                           const UInt step,
                                           const Double frequency,
                                           const ComplexFormat format);

  private:
    //! pointer to the ofstream with output -data
    std::ofstream * output;

    //! dataset 666. 
    void Dataset666();

    //! dataset 781
    void Dataset781();   

    //! dataset 780
    void Dataset780();

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
  
  };

} // end of namespace
 
#endif
