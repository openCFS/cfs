// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOU_NEUMANN_INT
#define FILE_ACOU_NEUMANN_INT

#include <string>

#include "Forms/linSurfForm.hh"
#include "General/defs.hh"
#include "General/environment.hh"

namespace CoupledField {
class EntityIterator;
struct Elem;
template <class TYPE> class Vector;
}  // namespace CoupledField

namespace CoupledField 
{
  
  //! Class implementing an inhomogeneous Neumann integrator for acoustic field
  class LinNeumannInt : public LinearSurfForm {

  public:
    
    //! Standard constructor
    LinNeumannInt( std::string amplitudeStr, std::string phaseStr,
                   SolutionType quantity, MaterialType materialParam,
                   bool isaxi );
    
    //! Destructor
    ~LinNeumannInt();

    /// Calculation of RHS vector for double entries, i.e. transient and static 
    void CalcElemVector( Vector<Double> & elemVec,
                         EntityIterator& ent );

    /// Calculation of RHS vector for complex entries, i.e. harmonic
    void CalcElemVector( Vector<Complex> & elemVec,
                         EntityIterator& ent );
    
    void SetAmplitude(const std::string& value) { amplitude_ = value; }
    
  protected:

    /// comprises part of CalcElemVector, which is equal for double and complex
    void PrepareElemVector( Vector<Double> & elemVec,
                            EntityIterator& ent );
    
    //! amplitude string
    std::string amplitude_;
    
    //! phase string
    std::string phase_;

    SolutionType quantity_;
    
    MaterialType materialParam_;

    Elem * ptVolElem_;

  };

} // end of namespace
#endif
