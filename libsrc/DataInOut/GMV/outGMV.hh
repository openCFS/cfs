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
  WriteResultsGMV(const Char * const filename,
		  FileType * const aInFile=NULL);
  
  //! Deconstructor
  virtual ~WriteResultsGMV();
  
  //! initialization with grid
  //! \param ptgrid pointer to grid object
  //! \param aptbcs pointer to BCs object
  virtual void Init(Grid * aptgrid, BCs *aptbcs);
  
  //! write information about grid with level in file
  /*!
    \param level level of the grid
  */
  virtual void WriteGrid(const Integer level);


  //! write element solution vector
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
  */
  virtual void WriteNodeSolutionTransient(const NodeStoreSol<Double>& sol, 
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
  virtual void WriteNodeSolutionHarmonic(const NodeStoreSol<Complex>& sol, 
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
  Integer currstep_;

  //! pointer to Grid
  Grid * ptgrid;

  //! indicator of type for data
  Boolean ascii_;

  //! indicator of adaptive grid or not
  Boolean fixedgrid_; 

  //! write header of gmv-file: only ascii is implemented
  void WriteHeader();

  //! write number of nodes and coordinates of them
  void WriteNodes(const Integer level);

  //! write cell description 
  void WriteCells(const Integer level); 

  //! write variable information
  /*!
    \param dataType data type of the var: 0.. cell data, 1.. node data, 2.. face data
     \param var vector with data
     \param name name of output-data
  */
  void WriteNodeVariable(const Vector<Double> var, const std::string name, const Integer dataType);
  
   //! write vector-variable information
  /*!
    \param dataType data type of the var: 0.. cell data, 1.. node data, 2.. face data
    \param var pointer to vector with output data
    \param name name of output-data
  */
  void WriteVelocity(const Vector<Double>* var, const std::string name, const Integer dataType);
 
  //! transform string to 8 characters. we need it, because name in gmv, in binary format, should be from 8 characters
  /*!
    \param name (input) title
    \param result (output) result
  */
  void to8Char(const std::string name, char * result);
  
  //! Convertes enum SolutionType to string
  std::string SolutionTypeToString(const SolutionType type) const;
  
};

} // end of namespace

#endif
