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

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

    /// Calculation of RHS vector for double entries, i.e. harmonic 
    virtual void CalcElemVector( Vector<Complex> & result,
                                 EntityIterator& ent );

  private:

    std::string id_;
    std::string regionName_;
    bool interpolate_;
    std::string srcInputId_;
    StdVector<std::string> srcRegions_;
    std::string restartFileMode_;
    std::string coordSysId_;
    Grid::ciTolerance globalEpsilon_;
    Grid::ciTolerance localEpsilon_;
    std::string asyncSteps_;
    UInt node_warnings_;
    StdVector<NodeList*> sourceNodeLists_;
    NodeList* destNodeList_;
    ElemList* destElemList_;
    BasePDE::AnalysisType analysistype_;
    Double z_;
    Double zEpsilon_;
    
    StdVector<Double> rhsValues_;
    StdVector<Double> rhsValues2_;
    StdVector<Complex> rhsValuesComplex_;
    bool isFirstStep_;
    UInt step_;
    UInt lastStep_;

    typedef std::vector< std::map<UInt, Double> > ConsInterpWeightsType;
    StdVector< ConsInterpWeightsType > consInterpWeights_;
    std::map<UInt, Double> stepValues_;

    MathParser::HandleType mHandle2_;
    MathParser::HandleType mph_tf_;
  };

} // end of namespace

#endif // FILE_ACOURHSLINFORM_HH
