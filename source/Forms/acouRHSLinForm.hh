// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOURHSLINFORM_HH
#define FILE_ACOURHSLINFORM_HH

#include "Utils/mathParser/mathParser.hh"

#include "linearForm.hh"

namespace CoupledField
{

  class ParamNode;
  class NodeList;
  class ElemList;
  
  // =============================================================================
  //  Read / Calculate RHS for acoustic PDE
  // =============================================================================

  class AcouRHSLinForm : public LinearForm
  {
  public:
    ///
    AcouRHSLinForm(ParamNode* rhsValuesNode);

    ///
    virtual ~AcouRHSLinForm();

    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

    //! Calculation of RHS vector for double entries, i.e. harmonic 
    virtual void CalcElemVector( Vector<Complex> & result,
                                 EntityIterator& ent );

  private:

    typedef enum {
      CI_WARN_NO = 0,
      CI_WARN_YES = 1,
      CI_WARN_VERBOSE = 2,
      CI_WARN_LIST = 4
    } ciWarnFlags;
    
    std::string id_; //!< input ID of destination grid
    std::string regionName_; //!< name of destination region
    bool interpolate_; //!< is interpolation performed or not?
    std::string srcInputId_; //!< input ID of sources
    StdVector<std::string> srcRegions_; //!< names of source regions
    std::string restartFileMode_; //!< usage of restart file: r, w, or rw
    std::string coordSysId_; //!< ID of coordinate system of source grid
    Grid::ciTolerance globalEpsilon_; //!< absolute interpolation tolerance
    Grid::ciTolerance localEpsilon_; //!< relative interpolation tolerance
    Double z_; //!< z coordinate of xy plane in 2d simulations
    Double zEpsilon_; //!< tolerance for interpolation from 3D to 2D
    std::string asyncSteps_; //!< do we use asynchronous time steps?
    ciWarnFlags node_warnings_; //!< flags for verbosity of interpolation warnings
    
    //! pointer to destination grid
    Grid* ptGrid_;
    //! pointer to global coordinate system
    CoordSystem* globCoSy_;
    //! type of analysis (transient/harmonic)
    BasePDE::AnalysisType analysistype_;
    //! number of time/frequency step
    UInt step_;
    //! number of previous time step
    UInt lastStep_;

    //! node numbers of source regions 
    StdVector<NodeList*> sourceNodeLists_;
    //! list of destination nodes
    NodeList* destNodeList_;
    //! list of destination elements
    ElemList* destElemList_;
    //! vector of nodal sources
    StdVector<Double> rhsValues_;
    //! 2nd vector of nodal sources (for interpolation in time)
    StdVector<Double> rhsValues2_;
    //! vector of complex nodal sources (for harmonic analysis)
    StdVector<Complex> rhsValuesComplex_;
    
    //! flag to remember if interpolation weights have been computed already
    bool isFirstStep_;
    //! type definition for storage of interpolation weights
    typedef std::vector< std::map<UInt, Double> > ConsInterpWeightsType;
    //! vector of interpolation weights
    StdVector< ConsInterpWeightsType > consInterpWeights_;
    //! mapping from step number to time/frequency
    std::map<UInt, Double> stepValues_;

    //! math parser handle to current time/frequency
    MathParser::HandleType mph_tf_;
    //! math parser handle to factor
    MathParser::HandleType mphFactor_;
  };

} // end of namespace

#endif // FILE_ACOURHSLINFORM_HH
