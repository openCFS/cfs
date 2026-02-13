// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_OUT_UNV_HH
#define FILE_CFS_SIM_OUT_UNV_HH

#include <fstream>
#include <string>

#include "Domain/Mesh/Grid.hh"
#include "DataInOut/SimOutput.hh"
#include "Domain/Results/ResultInfo.hh"

namespace CoupledField
{

  //! Class for writing output-information(grid,results) in unverg-format file
  
  class SimOutputUnv : virtual public SimOutput {
    
  public:
    
    //! Constructor
    SimOutputUnv( const std::string& fileName, PtrParamNode outputNode,
                  PtrParamNode infoNode, bool isRestart );
    
    //! Destructor
    virtual ~SimOutputUnv();
  
    //! Initialize class
    void Init( Grid* grid, bool printGridOnly );

    //! Write grid definition in file
    void WriteGrid();

    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd,
                         bool isHistory );

    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );
    
    //! Finish multisequence step
    void FinishMultiSequenceStep( );

  private:
    
    //! pointer to the ofstream with output -data
    std::ofstream * output;

    //! Map with result objects for each result type
    ResultMapType resultMap_;
    
    //! Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
            
    //! Output for CAPA .unverg if true, otherwise output standard I-DEAS .unv files.
    bool capaOut_;

    //! Dataset 151 - Header containing information about simulation run.
    void Dataset151();

    //! Dataset 164 - Header containing information about units.
    void Dataset164();

    //! Dataset 666. - CAPA specific information
    void Dataset666();

    //! Dataset 781 - Nodes in old CAPA format
    void Dataset781();   

    //! dataset 780 - Elements in old CAPA format
    void Dataset780();

    //! Dataset 2411 - Nodes in new I-DEAS format
    void Dataset2411();   

    //! Dataset 2412 - Elements in new I-DEAS format
    void Dataset2412();

    //! Dataset 2467 - Regions and Named Elements
    void Dataset2467();   

    void ElemType2UnvElemId(Elem::FEType et, UInt & id);

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
    void SolutionTypeToString(const SolutionType type,
                              std::string& capaSolType,
                              std::string& ideas5xSolType,
                              UInt& ideas2414SolType) const;
  
    //! Re-sort stresses according to vector with dofnames
    template<class TYPE>
    void SortStresses( Vector<TYPE>& vec,
                       const StdVector<std::string>& dofNames );
  };

} // end of namespace
 
#endif
