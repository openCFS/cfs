// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PMLBASICS
#define FILE_PMLBASICS

#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
namespace CoupledField
{

  //! forward class declarations
  class CoordSystem;

  //! class containing basic PML methods
  class PMLBasics
  {
  public:

    //! Constructor
    PMLBasics( std::string dampingTypePML, Double damp, std::string formsType);
    
    /// 
    virtual ~PMLBasics();
    
    //! calculates position and values
    void ComputeFactorPML( Vector<Complex>& factorsPML, Complex& pmlDet,
                           Vector<Double>& coordAtIP, Double omega );

   //! calculates position and values
    void ComputeFactorAPML( Vector<Complex>& factorsPML, 
                            Vector<Double>& coordAtIP, Double omega );
    
    //! calulates position and values for time domain PML
    void ComputeTimeFactorPML(Vector<Double>& factorsPML, 
                              Vector<Double>& pos );

    //! calulates position and values for time domain PML Depending on the center of the element
    //! Usefull if Lobatto integration is used
    void ComputeTimeFactorPML(Vector<Double>& factorsPML, 
                              Vector<Double>& pos,
                              Vector<Double> center);

    //! calculates the damping factor
    void ComputeDampingFactor( Vector<Double>& factors, 
                               const Vector<Double>& globPos,
                               const Vector<Double>& locPos,
                               UInt dir );
    
    //! set min/max of x,y,z coordinates form where PML starts and ends
    void SetPosPML( Matrix<Double> & inner, Matrix<Double> & outer,
                    const std::string& coordSysId );
    
    //! returns the type of bilinear form
    std::string GetFormsType() {
      return formsType_;
    };

  private:
    
    //! type of PML damping
    std::string dampingTypePML_;
    
    //! type of bilinear form
    std::string formsType_;
    
    //! damping factor 
    Double dampingFactor_;
    
    //! layer thickness (given in local coordinates)
    Matrix<Double> layerThickness_;
    
    //! id of coordinate system
    std::string coordSysId_;
    
    //! pointer to reference coordinate system
    CoordSystem * coordSys_;
    
    //! minimum coordinates of PML box (in local coordinateS)
    Vector<Double> minLoc_;
    
    //! maximum coordinates of PML box (in local coordinateS)
    Vector<Double> maxLoc_; 
    
  };

}

#endif // FILE_PMLBASICS
