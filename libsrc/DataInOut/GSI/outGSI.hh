#ifndef FILE_OUTGSI_2002
#define FILE_OUTGSI_2002

#include "Domain/grid.hh"
#include "DataInOut/writeresults.hh"

class CFS_ServerSocket;
class CFS_BaseIO;

namespace CoupledField
{
//! This class provides an interface for writing files in gmv-format
class WriteResultsGSI: virtual public WriteResults
{
public:

  //! Constructor
  WriteResultsGSI(const Char * const filename,Boolean withHistory=FALSE);

  //! Deconstructor
  virtual ~WriteResultsGSI();
  
  int waitForConnection();
  //! initialization with grid
  /*!
    \param aptgrid pointer to class Grid
  */
  virtual void Init(Grid * aptgrid);

  //! write information about grid with level in file
  /*!
    \param level level of the grid
  */
  virtual void WriteGrid(const Integer level);

  //! write information about the solution
  /*!
    \param sol solution
    \param step number of step of the calculation
    \param time time of the calculation
     \param title name for the solution
  */
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title);

 //! write cell data
 /*!
    \param data vector with data (ex. value of an error for the cell)
    \param step step of calculation
    \param time time of calculation
    \param title name for the data
  */
  virtual void WriteDataOnCell(const Vector<Double> & data, const Integer step, const Double time, const std::string title);

  //! write vectors-data on cells
  /*!
    \param vec pointer to vector with data
     \param step step of calculation
    \param time time of calculation
    \param title name for the data
  */
  virtual void WriteResultsGSI::WriteVecDataOnCell(const Vector<Double>*vec,const Integer step, const Double time, const std::string title);

 //! write comments
  /*!
    \param comments string with comments
  */
 virtual void WriteComments(const std::string comments){;}

  //! check, is it the gmv-output file
  virtual Boolean IsGMV(){ return FALSE; }

 //! function for open file with number num 
  void OpenFile(const Integer num) {};

  void WriteMsg(const std::string msg);
  Integer DoCompute();
  Integer TransferGrid();

private:
  // number of step
  Integer currstep_;
  Integer maxnumnodes_;

  //! pointer to Grid
  Grid * ptgrid;

  //! write header of gmv-file: only ascii is implemented
  //  void WriteHeader();

  //! write number of nodes and coordinates of them
  //  void WriteNodes(const Integer level);

  //! write cell description 
  //  void WriteCells(const Integer level); 

  //! write variable information
  /*!
    \param dataType data type of the var: 0.. cell data, 1.. node data, 2.. face data
     \param var vector with data
     \param name name of output-data
  */
  //  void WriteVariable(const Vector<Double> var, const std::string name, const Integer dataType);
  
   //! write vector-variable information
  /*!
    \param dataType data type of the var: 0.. cell data, 1.. node data, 2.. face data
    \param var pointer to vector with output data
    \param name name of output-data
  */
  //  void WriteVelocity(const Vector<Double>* var, const std::string name, const Integer dataType);
 
  //! transform string to 8 characters. we need it, because name in gmv, in binary format, should be from 8 characters
  /*!
    \param name (input) title
    \param result (output) result
  */
  //  void to8Char(const std::string name, char * result);

  void  Dataset666(const Integer);
  void  WriteResultsGSI::Dataset781(const Integer);
  void  WriteResultsGSI::Dataset780(const Integer);
  void  WriteResultsGSI::Dataset55(const std::string &, const Vector<Double> &, const Integer step, const Double time);
  void  WriteResultsGSI::Dataset56_Scalar(const std::string &, const Vector<Double> &, const Integer, const Double);
  void  WriteResultsGSI::Dataset56_Vec(const std::string &, const Vector<Double> *, const Integer, const Double);
  

private:
  CFS_ServerSocket *ss_;
  CFS_ServerSocket *sock_;
  CFS_BaseIO *io_;
};

} // end of namespace

#endif
