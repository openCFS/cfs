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
  WriteResultsUnverg(const Char * const filename, 
		     FileType * const aInFile=NULL); 

  //! deconstructor
  virtual ~WriteResultsUnverg();

  //! initialization with grid
  //! \param ptgrid pointer to grid object
  //! \param aptbcs pointer to BCs object
  virtual void Init(Grid * aptgrid, BCs * aptbcs);

  //! write information about grid on level i in file
   /*!
    \param level in: level of grid
  */
  virtual void WriteGrid(const Integer level);

  //! write element solution vector
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
  */
  virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& data, 
					  const Integer step, 
					  const Double time);

 //! write element solution vector
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
  */
  virtual void WriteElemSolutionTransient(const ElemStoreSol<Double>& data, 
					  const Integer step, 
					  const Double time);
 
  //! write element solution vector 
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param frequency frequency of exciting function
    \param format format for writing complex solution (real-imag/amplitude-phase)
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
    \param format format for writing complex solution (real-imag/amplitude-phase)
  */
  virtual void WriteElemSolutionHarmonic(const ElemStoreSol<Complex>& data, 
					 const Integer step,
					 const Double frequency,
					 const ComplexFormat format);

  //!  check, is it the gmv-output file
  virtual Boolean IsGMV() { return FALSE;}

private:
  //! pointer to the ofstream with output -data
  std::ofstream * output;

  //! dataset 666. 
  /*!
    \param level level of the Grid
  */
  void Dataset666(const Integer level);

  //! dataset 781
  /*!
    \param level level of the Grid
  */
  void Dataset781(const Integer level);   

  //! dataset 780
  /*!
     \param level level of the Grid
  */
  void Dataset780(const Integer level);

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
  
};

} // end of namespace
 
#endif
