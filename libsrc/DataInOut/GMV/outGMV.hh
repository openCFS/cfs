#ifndef FILE_OUTGMV_2002
#define FILE_OUTGMV_2002

#include "grid.hh"
#include "writeresults.hh"

namespace CoupledField
{

//! This class provides an interface for writing files in gmv-format
class WriteResultsGMV: virtual public WriteResults
{
public:

  //! Constructor
  WriteResultsGMV(const Char * const filename,Boolean withHistory=FALSE);
  
  //! Deconstructor
  virtual ~WriteResultsGMV();
  
  //! initialization with grid
  virtual void Init(Grid * aptgrid);

  //! write information about grid with level in file
  virtual void WriteGrid(const Integer level);

  /// write information about the solution
  virtual void WriteSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title);

 //! write cell data
  virtual void WriteDataOnCell(const Vector<Double> & data, const Integer step, const Double time, const std::string title);

  //! write vectors-data
  void WriteResultsGMV::WriteVecDataOnCell(const Vector<Double>*vec,const Integer step, const Double time, const std::string title);

 //! write comments
 virtual void WriteComments(const std::string comments){;}

  //!
  virtual Boolean IsGMV(){ return TRUE;}

 //! function for open file with number num 
  void OpenFile(const Integer num);

private:
  ///
  std::ofstream * output;

  //! name of dir for output results
  Char * namedir_;

  // number of step
  Integer currstep_;

  //!
  Grid * ptgrid;

  //! indicator of type for data
  Boolean ascii_;

  //! write header of gmv-file: only ascii is implemented
  void WriteHeader();

  //! write number of nodes and coordinates of them
  void WriteNodes(const Integer level);

  //! write cell description 
  void WriteCells(const Integer level); 

  //! write variable information
  /*!
    \param dataType data type of the var: 0.. cell data, 1.. node data, 2.. face data
  */
  void WriteVariable(const Vector<Double> var, const std::string name, const Integer dataType);
  
   //! write vector-variable information
  /*!
    \param dataType data type of the var: 0.. cell data, 1.. node data, 2.. face data
  */
  void WriteVelocity(const Vector<Double>* var, const std::string name, const Integer dataType);
 
};

} // end of namespace

#endif
