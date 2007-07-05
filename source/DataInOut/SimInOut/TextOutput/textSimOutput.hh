// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_TEXT_SIM_OUTPUT_HH
#define FILE_CFS_TEXT_SIM_OUTPUT_HH

#include "Domain/grid.hh"
#include "DataInOut/simOutput.hh"
#include "Domain/entityList.hh"

namespace CoupledField {


  // forward class declaration
  class CoordSystem;

  class SimOutputText : public SimOutput {

  public:

    //! Enum defining type of collection
    typedef enum { ENTITY, TIMEFREQ} CollectionType;

    //! Constructor
    SimOutputText( const std::string& fileName, ParamNode * outputNode );

    //! Destructor
    ~SimOutputText();

    //! Initialize class
    void Init( Grid * ptGrid, bool printGridOnly );

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

    //! Create new output stream
    void CreateFiles( shared_ptr<BaseResult> res,
                      UInt step,
                      Double stepVal );
    
    //! Write step for collectionType ENTITY
    void WriteStepCollectEntity();
    
    //! Write step for collecitonType TIMEFREQ
    void WriteStepCollectTimeFreq();

    //! Create header information depending on result
    std::string ResultDofString( shared_ptr<BaseResult> res );

    // =======================================================================
    //  Helper functions to iterate over a list of entities
    // =======================================================================    

    //! Return for an entityiterator the node number and coordinates
    std::string GetNodeInfo( EntityIterator & it ) const;
    
    //! Return for an entityiterator the node number and coordinates
    std::string GetElemInfo( EntityIterator & it ) const;
    
    //! Return for an entityiterator the node number
    std::string GetRegionInfo( EntityIterator & it ) const;

    // =======================================================================
    //  Helper function to extract complex values in correct format
    // =======================================================================

    //! Extract complex value as Amplitude-Phase
    std::string ComplexAsAmplPhase( const Complex& val ) const;

    //! Extract complex value as Real-Imag
    std::string ComplexAsRealImag( const Complex& val ) const;

    //! type of collection
    CollectionType collecType_;

    //! Map with result objects for each result type
    ResultMapType resultMap_;
    
    //! Map of result and vector of ofstreams
    std::map<shared_ptr<BaseResult>,
             StdVector<std::ofstream*> > outFiles_;

    //! Comment char
    char cmChar_;

    //! Delimiter string
    static std::string delim_;

    //! Coordinate system
    CoordSystem * coordSys_;

    //! 
  };
}



#endif
