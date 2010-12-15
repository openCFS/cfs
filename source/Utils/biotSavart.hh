// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include <map>
#include <string>
#include "Utils/StdVector.hh"
#include "Domain/grid.hh"
#include "PDE/eqnMap.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/mathParser/mathParser.hh"

#ifndef FILE_CFS_BIOTSAVART_HH
#define FILE_CFS_BIOTSAVART_HH

namespace CoupledField{

  //! Class implementing fundamental solution of magnetic field in air.
  
  //! This class can calculate either the magnetic vector potential A or the
  //! magnetic field strength H using the law of Biot-Savart. 
  //! The definition of the coils is specified using so-called current sticks,
  //! i.e. piece-wise straight, infinite thin current carrying lines.
  //! The definition of the coil is read in from an external text file. 
  class BiotSavart {

public:
   
    //! Typedef for formulation
    typedef enum {UNDEF, VEC_POT, MAG_FIELD_STRENGTH} FormulationType;
    
    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================

    //! Constructor
    BiotSavart( );

    //! Desturctor
    ~BiotSavart();

    //! Initialize from parameter node
    void Init( PtrParamNode currentNode, Grid *ptGrid,
               shared_ptr<EqnMap> eqnMap );
    
    //! Generate visual representation of current sticks in grid
    void GenGridRepresentation();
    
    //! Set formulation for the Biot Savart class (can be changed at any time)
    void SetFormulation( FormulationType formulation );
    
    //! Return type of formulation
    FormulationType GetFormulation() {
      return formulation_;
    }
    
    // =======================================================================
    //  CALCULATION / ACCESS METHODS
    // =======================================================================
    
    //! Get field for all equations
    Vector<Double>& CalcFieldAllEqns( bool normalized );
    
    //! Get derivative of field on all local nodes
    Vector<Double>& CalcFieldDeriv1AllEqns(FormulationType type);
    
    //! Get field for one single equation
    Double CalcFieldSingleEqn( UInt eqn ); 
    
    //! Get field for arbitrary observer point
    void CalcFieldAtPoint( Vector<Double>& hField, 
                           const Vector<Double>& observer );
    
  private:

    //! Helper structure for one biot savart coil (current stick)
    struct BsCoil{

      //! Name of coil
      std::string name;

      //! Store points of coil 
      StdVector<Vector<Double> > points_;

      //! Normalized field for all nodes 
      Vector<Double> fieldNormalized_;

      //! Handle for MathParser object
      MathParser::HandleType mHandle_;
    };

    
    //! Read in data from file
    void ReadFile(std::string fileName, BsCoil& coil ); 
    
    //! Computes the field (A/H) according to Biot-Savart for one point / coil
    
    //! This method computes the contribution of ONE single coil for the 
    //! field (A/H) for one observer point. 
    void ComputeFieldNormalized( Vector<Double>& field, 
                                 const Vector<Double>& observer,
                                 const BsCoil& coil );
    
    // =======================================================================
    // VECTOR POTENTIAL FORMULATION
    // =======================================================================
    
    //! Compute vector potential according to Biot-Savart along an arc
    void ArcIntegralVecPot( Vector<Double>& partMagVec, 
                            const Vector<Double>& observer,
                            const Vector<Double>& start, 
                            const Vector<Double>& end );

    //! Compute vector potential according to Biot-Savart along a straight line

    //! Compute contribution of Biot-Savart field of a current stick, i.e. a
    //! straight line
    void LineIntegralVecPot( Vector<Double>& partMagVec, 
                             const Vector<Double>& observer,
                             const Vector<Double>& start, 
                             const Vector<Double>& end );

       
    //! Kernel for vector potential evaluation
    void KernelVecPot( const Vector<Double>& observer, Vector<Double>& p,
                       Vector<Double>& dir, Vector<Double>& kernelVP);
    
    // =======================================================================
    // SCALAR POTENTIAL FORMULATION
    // =======================================================================

    //! Compute H-field according to Biot-Savart along a straight line

    //! Compute contribution of Biot-Savart field of a current stick, i.e. a
    //! straight line to the magnetic field strength H for a given observer 
    //! point
    void LineIntegralHField( Vector<Double>& partMagVec, 
                             const Vector<Double>& observer,
                             const Vector<Double>& start, 
                             const Vector<Double>& end );
   
    //! Computes the new time loading value (needed?)
    void SetTimeFncValue();


    //! Formulation of magnetic field
    FormulationType formulation_;
    
    //! Pointer to function for evaluating the H/A-field of a single coil-segment
    
    //! This function pointer is used to point to a generic function, which 
    //! calculates on of the following quantities per coil-segment:
    //! - line integral for vector potential A
    //! - arc integral for vector potential A
    //! - line integral for magnetic field strength H
    //! - arc integral for magnetic field strength H
    void (BiotSavart::*ptSegmentFunc_)( Vector<Double>&,
                                        const Vector<Double>&,
                                        const Vector<Double>&, 
                                        const Vector<Double>& );
    
    //! Dimension of problem
    UInt dim_;// for 2D case , it will not call the view function currently
        
    
    //! Number of vector components
    
    //! This variable will be set depending on the formulation and the grid
    //! dimension.   
    UInt numVecComponents_;
    
    //! Flag if normalized field is already computed
    bool fieldIsComputed_;
    
    //! Flag for axi-symmetry
    bool isAxi_;
    
    //! List of current sticks
    StdVector<BsCoil> coils_;    
    
    //! Field (A / H) for current time step 
    Vector<Double> field_;
    
    //! Parameter node
    PtrParamNode myParam_;

    //! Pointer to eqnMap
    shared_ptr<EqnMap> eqnMap_;
    
    //! Pointer to grid
    Grid * ptGrid_;

  };


} //end of namespace
#endif
