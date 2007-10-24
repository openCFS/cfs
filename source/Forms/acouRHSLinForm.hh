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
    AcouRHSLinForm(std::string id,
                   std::string regionName,
                   std::string factor,
                   bool interpolate,
                   std::string srcInputId,
                   std::string srcRegion,
                   std::string coordSysId,
                   Double globalEpsilon,
                   Double localEpsilon,
                   std::string restartFileMode);

    ///
    virtual ~AcouRHSLinForm();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );

    /// Calculation of RHS vector for double entries, i.e. harmonic 
    virtual void CalcElemVector( Vector<Complex> & result,
                                 EntityIterator& ent );

  private:

    std::string regionName_;
    std::string id_;
    bool interpolate_;
    std::string srcInputId_;
    StdVector<std::string> srcRegions_;
    std::string restartFileMode_;
    std::string coordSysId_;
    Double globalEpsilon_;
    Double localEpsilon_;
    StdVector<NodeList*> sourceNodeLists_;
    NodeList* destNodeList_;
    ElemList* destElemList_;
    
    std::vector<Double> rhsValues_;
    std::vector<Complex> rhsValuesComplex_;
    UInt step_;
    UInt lastStep_;

    typedef std::vector< std::map<UInt, Double> > ConsInterpWeightsType;
    StdVector< ConsInterpWeightsType > consInterpWeights_;

    MathParser::HandleType mHandle2_;
  };

} // end of namespace

#endif // FILE_ACOURHSLINFORM_HH
