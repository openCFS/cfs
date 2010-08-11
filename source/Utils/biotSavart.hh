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
  class BiotSavart {

public:
   
    //! Typedef for formulation
    typedef enum {VEC_POT, SCALAR_POT} FormulationType;
    
    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================

    //! Constructor
    BiotSavart( FormulationType type );

    //! Desturctor
    ~BiotSavart();

    //! Initialize from parameter node
    void Init( PtrParamNode currentNode, Grid *ptGrid,
               shared_ptr<EqnMap> eqnMap );
    
    //! Generate visual representation of current sticks in grid
    void GenGridRepresentation();
    
    //! Return type of formulation
    FormulationType GetFormulation() {
      return formulation_;
    }
    
    // =======================================================================
    // VECTOR POTENTIAL FORMULATION
    // =======================================================================
    
    //! Returns the current magnetic vector potential value
    Double GetMagVec( UInt eqn, BasePDE::AnalysisType analysis );

    //! Returns the normalized vector potential
    Vector<Double>& GetMagVec(bool normalized );
    
    //! Return 1st time derivative of magnetic vector potentential
    Vector<Double>& GetMagVecDeriv1();
    
    // =======================================================================
    // SCALAR POTENTIAL FORMULATION
    // =======================================================================
    
    //! Calculates the magnetic field H at the given midpoints.
    void CalcH( Vector<Double>& HField,Vector<Double>& intPoint);
    

    
  private:

    //! Helper structure for one biot savart coil (current stick)
    struct BsCoil{

      //! Name of coil
      std::string name;

      //! Store points of coil 
      StdVector<Vector<Double> > points_;

      //! Normalized vector potential for all nodes 
      Vector<Double> magVecNormalized_;

      //! Handle for MathParser object
      MathParser::HandleType mHandle_;
    };

    
    //! Read in data from file
    void ReadFile(std::string fileName, BsCoil& coil ); 
    
    //! compute the magnetic vector potential
    void ComputeBiotSavartField();
    
    //! computes the magnetic vector potential according to Biot-Savart
    void ComputeMagVecNormalized( Vector<Double>& magVec, 
                                  const Vector<Double>& observer,
                                  const BsCoil& coil );

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

   
    //! Computes the new time loading value (needed?)
    void SetTimeFncValue();


    //! Formulation of magnetic field
    FormulationType formulation_;
    
    //! Dimension of problem
    UInt dim_;// for 2D case , it will not call the view function currently
        
    //! Flag if normalized field is already computed
    bool fieldIsComputed_;
    
    //! Flag for axi-symmetry
    bool isAxi_;
    
    //! List of current sticks
    StdVector<BsCoil> coils_;    
    
    //! Magnetic vector potential for current time step 
    Vector<Double> magVec_;
    
    //! Parameter node
    PtrParamNode myParam_;

    //! Pointer to eqnMap
    shared_ptr<EqnMap> eqnMap_;
    
    //! Pointer to grid
    Grid * ptGrid_;
    
    

  };


} //end of namespace
#endif
