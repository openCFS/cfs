#ifndef FILE_OUTRESULTDATABASE_2004
#define FILE_OUTRESULTDATABASE_2004



#include "database.hh"
#include "dbMatrix.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/Database/dbLineData.hh"
#include "DataInOut/writeresults.hh"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace CoupledField
{

//! Class for writing output-information(grid,results) into a database

class WriteResultsDatabase: virtual public WriteResults
{

public:
//! Default constructor
WriteResultsDatabase(const Char * const filename, FileType * const aInFile=NULL);

//! Destructor 
~WriteResultsDatabase();

//! Initialization with grid
virtual void Init (Grid *aptgrid, BCs * aptbcs);

//! Write basis information (Calculation parameter, configuration file..)
void WriteBasisData();

//! write information about grid on level i in database
/*!
  param level Level of the grid
*/
virtual void WriteGrid (const Integer level);

//! write nodal solution vector
/*!
  \param data vector with data (ex. value of an error for the cell)
  \param step step of calculation
  \param time time of calculation
  \param title name for the data    
*/
virtual void WriteNodeSolutionTransient (const NodeStoreSol<Double> &sol, const Integer step, const Double time);

//! write element solution vector
/*!
  \param data vector with data (ex. value of an error for the cell)
  \param step step of calculation
  \param time time of calculation
  \param title name for the data    
*/
virtual void WriteElemSolutionTransient (const ElemStoreSol<Double>&sol, const Integer step, const Double time);

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
					 const ComplexFormat format)
{Warning("WriteNodeSolutionHarmonic not yet implemented",__FILE__,__LINE__);};

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
					 const ComplexFormat format)
{Warning("WriteElemSolutionHarmonic not yet implemented",__FILE__,__LINE__);};

  //! write node history vector (transient/static)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
  */
  void WriteNodeHistoryTransient(const NodeStoreSol<Double>& data, 
				 const Integer step, 
				 const Double time)
{Error("WriteNodeHistoryTransient not yet implemented.",__FILE__,__LINE__);};


  //! write nodal history vector (harmonic)
  /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param frequency frequency of exciting function
    \param frequencyStep step of calculation
    \param format format for writing complex solution (real-imag/amplitude-phase)
  */
  void WriteNodeHistoryHarmonic(const NodeStoreSol<Complex>& data, 
				const Integer step,
				const Double frequency,
				const ComplexFormat format)
{Error("WriteNodeHistoryHarmonic not yet implemented.",__FILE__,__LINE__);};

//! check, is it the gmv-output file
Boolean IsGMV();

//! Store config-file in database
void WriteConfFile();

protected:

//! Connection to Database
Database Db_;

//! Dummy function
void InitHistoryFiles();

//! Convertes enum SolutionType to string
std::string SolutionTypeToString(const SolutionType type) const;

//! Save information of finite/boundary elements
/*!
  param level Level of the grid
*/
long int Dataset780(const Integer level);

//! Save coordinates of the used grid
/*!
  param level Level of the grid
*/
long int Dataset781(const Integer level);

//! Save nodal results
/*!
  \param title Title of the results
  \param x Array with cell results
  \param step Number of the step of the calculation
  \param time Time of the calculation
  \param nrDofs Number of degrees of freedom
*/
void Dataset55(const std::string & title, 
	       const Vector<Double> & x, 
	       const Integer step, 
	       const Double time, 
               const Integer nrNodes,
	       const Integer nrDofs);

//! Save elements results
/*!
  \param title Title of the results
  \param x Array with cell results
  \param step Number of the step of the calculation
  \param time Time of the calculation
  \param nrDofs Number of degrees of freedom
*/
void Dataset56(const std::string &title, 
	       const Vector<Double> & x, 
	       const Integer step, 
	       const Double time, 
               const Integer nrNodes,
	       const Integer nrDofs);

//! Index of parameter file
unsigned long int InputParamIdx_;

//! Index in entity "result"
unsigned long int ResultIdx_;

//! Index in entity "Calculation"
unsigned long int CalculationIdx_;

//void AddInHistory(const Double time, const Double val, const int nodeidx);

//void AddVecInHistory(const Double time, const std::vector<Double> val, const int nodeidx);

//std::vector<unsigned long int> nodehistoryindex;

//std::vector<unsigned long int> elementhistoryindex;

}; // End class WriteResultsDatabase

} // End of namespace

#endif

