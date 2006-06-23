#ifndef FILE_OUTANSYS_2003
#define FILE_OUTANSYS_2003

#include "Domain/grid.hh"
#include "DataInOut/writeresults.hh"

namespace CoupledField
{

  //! This class provides an interface for writing files in ansys-format
  class WriteResultsAnsys: virtual public WriteResults
  {
  public:

    //! Constructor
    WriteResultsAnsys(const Char * const filename);
  
    //! Deconstructor
    virtual ~WriteResultsAnsys();
  
    //! initialization with grid
    /*!
      \param aptgrid pointer to class Grid
    */
    virtual void Init(Grid * aptgrid);

    //! write information about grid with level in file
    virtual void WriteGrid();

    //! write information about the solution
    /*!
      \param sol solution
      \param step number of step of the calculation
      \param time time of the calculation
      \param title name for the solution
    */
    virtual void WriteNodeSolution(const Vector<Double> & sol, const Integer step, const Double time, const std::string title);

    //! write cell data
    /*!
      \param data vector with data (ex. value of an error for the cell)
      \param step step of calculation
      \param time time of calculation
      \param title name for the data
    */
    virtual void WriteElemSolution(const Vector<Double> & data, const Integer step, const Double time, const std::string title);

    //! write vectors-data on cells
    /*!
      \param vec pointer to vector with data
      \param step step of calculation
      \param time time of calculation
      \param title name for the data
    */
    void WriteVecDataOnCell(const Vector<Double>*vec,const Integer step, const Double time, const std::string title);

    //! write comments
    /*!
      \param comments string with comments
    */
    virtual void WriteComments(const std::string comments){;}

    //! check, is it the gmv-output file
    virtual bool IsGMV(){ return false;}

    //! open files for writing grid
    void OpenFile(const Integer numstep);

  private:
  
    //! pointer to ofstream with nodes-information
    std::ofstream * out_nodes_;
  
    //! pointer to ofstream with cells-information
    std::ofstream * out_cells_;

    //! pointer to ofstream with results-information
    std::ofstream * out_res_;

    //! name of dir for output results
    Char * namedir_;

    // counter
    Integer currnum_;

    //! pointer to Grid
    Grid * ptgrid;

    // open files with results output
    void OpenFileRes(const std::string title);

    //! write number of nodes and coordinates of them
    void WriteNodes();

    //! write cell description 
    void WriteCells(); 

    //! write variable information
    /*!
      \param var vector with data
    */
    void WriteVariableNode(const Vector<Double> var);

    //! write variable information
    /*!
      \param var vector with data
    */
    void WriteVariableCell(const Vector<Double> var, const std::string title);

  
    //! write vector-variable information
    /*!
      \param dataType data type of the var: 0.. cell data, 1.. node data, 2.. face data
      \param var pointer to vector with output data
      \param name name of output-data
    */
    void WriteVelocity(const Vector<Double>* var, const std::string name, const Integer dataType);
 
  };

} // end of namespace

#endif
