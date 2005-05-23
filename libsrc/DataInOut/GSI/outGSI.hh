#ifndef FILE_OUTGSI_2004
#define FILE_OUTGSI_2004

#include "Domain/grid.hh"
#include "DataInOut/writeresults.hh"

// forward declarations
namespace GridlibSocketInterface
{
    class ServerSocket;
    class BaseIO;
}


namespace CoupledField
{

  //! This class provides an interface for writing files in GSI-format

  class WriteResultsGSI: virtual public WriteResults {

  public:

    //! Constructor
    WriteResultsGSI( const Char *const filename );

    //! Deconstructor
    virtual ~WriteResultsGSI();
  
    //! initialization with grid
    /*!
      \param aptgrid pointer to class Grid
    */
    virtual void Init(Grid * aptgrid);

    //! write information about grid in file
    virtual void WriteGrid();

    //! write node solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
    */
    virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& data, 
					    const Integer step, 
					    const Double time);

    //! write node solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
      \param title name for the data
    */
    virtual void WriteElemSolutionTransient(const ElemStoreSol<Double>& data, 
					    const Integer step, 
					    const Double time);
 
    //! write element solution vector 
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param frequency frequency of exciting function
      \param format format for writing complex solution
             (real-imag/amplitude-phase)
    */
    virtual void WriteNodeSolutionHarmonic(const NodeStoreSol<Complex>& data, 
					   const Integer step,
					   const Double frequency,
					   const ComplexFormat format);

    //! write element solution vector
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param frequency frequency of exciting function
      \param format format for writing complex solution
             (real-imag/amplitude-phase)
    */
    virtual void WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& data, 
					   const Integer step,
					   const Double frequency,
					   const ComplexFormat format);

    //! \todo document me
    Integer DoCompute();

    //! \todo document me
    Integer TransferGrid();

  private:

    //! number of step
    Integer currstep_;

    //! \todo document me
    Integer maxnumnodes_;

    //! dataset 666. 
    void Dataset666();

    //! dataset 781
    void Dataset781();   

    //! dataset 780
    void Dataset780();

    //! for printing nodal results of simulation (static/transient)
    /*!
      \param title title of the results.
      \param x array with nodal results
      \param step number of the step of the calculation
      \param time time of the calculation
    */
    void Dataset55_Transient(const std::string & title, 
			     const Vector<Double> & x, 
			     const Integer step, 
			     const Double time, 
			     const Integer nrNodes,
			     const Integer nrDofs=1);
  
    //! for printing nodal results of simulation (harmonic)
    /*!
      \param title title of the results.
      \param x array with nodal results
      \param freuqncy exciting frequency of current result
      \param format output format for complex numbers
    */
    void Dataset55_Harmonic(const std::string & title, 
			    const Vector<Complex> & x, 
			    const Integer step,
			    const Double frequency,
			    const ComplexFormat format,
			    const Integer nrNodes,
			    const Integer nrDofs=1);
  
    //! for printing cell results of simulation (transient / static)
    /*!
      \param title title of the results.
      \param x array with cell results
      \param step number of the step of the calculation
      \param time time of the calculation
    */
    void Dataset56_Transient(const std::string & title, 
			     const Vector<Double> & x, 
			     const Integer step, 
			     const Double time,
			     const Integer numElems,
			     const Integer nrDofs=1);
  
    //! for printing cell results of simulation (harmonic)
    /*!
      \param title title of the results.
      \param x array with cell results
      \param step number of the step of the calculation
      \param time time of the calculation
    */
    void Dataset56_Harmonic(const std::string & title, 
			    const Vector<Complex> & x, 
			    const Integer step,
			    const Double frequency, 
			    const ComplexFormat format, 
			    const Integer numElems,
			    const Integer nrDofs=1);
  
    //! Convertes enum SolutionType to string
    std::string SolutionTypeToString(const SolutionType type) const;
  
    void WriteMsg(const std::string msg);

  private:

    //GSIServerSocket *ss_;
    //GSIServerSocket *sock_;

    //! \todo document me
    GridlibSocketInterface::BaseIO *io_;

    //! \todo document me
    FILE* fp_;

  };

} // end of namespace

#endif
