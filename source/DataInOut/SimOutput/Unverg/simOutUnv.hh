#ifndef FILE_CFS_SIM_OUT_UNV_HH
#define FILE_CFS_SIM_OUT_UNV_HH

#include "Domain/grid.hh"
#include "DataInOut/simOutput.hh"
#include "Domain/resultInfo.hh"

namespace CoupledField
{

  //! Class for writing output-information(grid,results) in unverg-format file
  
  class SimOutputUnv : virtual public SimOutput {
    
  public:
    
    //! Constructor
    SimOutputUnv( const std::string& fileName, ParamNode * outputNode );
    
    //! Destructor
    virtual ~SimOutputUnv();
  
    //! Initialize class
    void Init( Grid* grid );

    //! Write grid definition in file
    void WriteGrid();

    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd );

    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );

  private:
    
    //! pointer to the ofstream with output -data
    std::ofstream * output;

    //! Map with result objects for each result type
    ResultMapType resultMap_;

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
  
    //! Re-sort stresses according to vector with dofnames
    template<class TYPE>
    void SortStresses( Vector<TYPE>& vec,
                       const StdVector<std::string>& dofNames );
  };

} // end of namespace
 
#endif
