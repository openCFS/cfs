// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PIERCEINT_1
#define FILE_PIERCEINT_1

#include "baseForm.hh"

namespace CoupledField
{

  // forward class declaration
  class SimpleFlow;
  class ParamNode;

  //! Class for calculation  element damping matrix for Pierce formulation
  class PierceDampInt : public BaseForm
  {
  public:

    // Constructor
    PierceDampInt(const Double aFactor, SimpleFlow * flow,
		  bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~PierceDampInt();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
      

  protected:
    
  private:

    Double dampFactor_;           //!< multiplicative value for damping integrator

    SimpleFlow* flowHandler_;
  };


  //! Class for calculation  element stiffness matrix for Pierce formulation
  class PierceStiffInt : public BaseForm
  {
  public:

    // Constructor
    PierceStiffInt(const Double aFactor, SimpleFlow* flow,
		   bool axi=false, bool coordUpdate = false );

    // Destructor
    virtual ~PierceStiffInt();

    // Calculation of element matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
      

  protected:
    
  private:

    Double stiffFactor_;           //!< multiplicative value for stiff integrator
    SimpleFlow* flowHandler_;
  };

  //! Class for simple flow computation
  class SimpleFlow
  {
  public:

    //! Constructor
    SimpleFlow();

    //! Destructor
    virtual ~SimpleFlow();

    //! reads in flow data
    void ReadFlowData( ParamNode * flowNode, UInt dim );

    //! computes the flow velocity and derivatives
    void ComputeActFlow( const Vector<Double> coord, 
			 Vector<Double>& velocity, 
			 Vector<Double>& derivVel);

  protected:
    
  private:
    //! for flow Data
    Directions flowDir_; //!< direction of flow:  X, Y, Z
    UInt   flowOrder_;   //!< power factor in flow formula
    Double flowVelocity_; //!< maximal velocity
    Double flowRadius_;  //!< radius of flow
    Double flowLength_;  //!< length of flow
    Double flowDecay_;   //!< length value for smoothing flow to zero
    Vector<Double> flowCenter_; //!< center of flow

  };


}

#endif // FILE_PIERCEINT
