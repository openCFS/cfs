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
  WriteResultsUnverg(const Char * const filename,Boolean withHistory=FALSE, FileType * const aInFile=NULL); 

  //! deconstructor
  virtual ~WriteResultsUnverg();
  
   //! initialization with grid
  virtual void Init(Grid * aptgrid);

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
    \param title name for the data    
    \param nrDofs dimension of solution
  */
  virtual void WriteNodeSolution(const NodeStoreSol<Double>& sol, const Integer step, const Double time, const std::string title);

 //! write element solution vector
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
    \param title name for the data    
    \param nrDofs dimension of solution
  */
  virtual void WriteElemSolution(const ElemStoreSol<Double>& data, const Integer step, const Double time, const std::string title);
  
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

  //! for printing nodal results of simulation
  /*!
    \param title title of the results.
    \param x array with nodal results
    \param step number of the step of the calculation
    \param time time of the calculation
  */
  void Dataset55(const std::string & title, 
		 const Vector<Double> & x, 
		 const Integer step, 
		 const Double time, 
		 const Integer nrNodes,
		 const Integer nrDofs=1);

  //! for printing cell results of simulation
   /*!
    \param title title of the results.
    \param x array with cell results
    \param step number of the step of the calculation
    \param time time of the calculation
  */
  void Dataset56(const std::string & title, 
		 const Vector<Double> & x, 
		 const Integer step, 
		 const Double time, 
		 const Integer numElems,
		 const Integer nrDofs=1);

};

} // end of namespace
 
#endif
