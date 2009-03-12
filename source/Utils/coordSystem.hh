// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COORD_SYSTEM_HH
#define COORD_SYSTEM_HH

#include <string>

#include "Domain/grid.hh"
#include "Utils/tools.hh"
#include "MatVec/vector.hh"
#include "MatVec/matrix.hh"
#include "General/environment.hh"

namespace CoupledField {

  // forward class declarations
  class ParamNode;

  //! Base class for describing a local coordinate system
  class CoordSystem {

  public:

    //! Default constructor
    //! \param name (in) name of the coordinate system which is used to 
    //!                  find the additional parameters in the .xml-file
    //! \param ptGrid (in) pointer to finite element grid object
    //! \param myParamNode (in) pointer to parameter node of current coosy
    CoordSystem(const std::string & name, Grid * ptGrid,
                ParamNode * myParamNode );

    //! Destructor
    virtual ~CoordSystem();
    
    //! Return the name of the coordinate system
    const std::string & GetName() const {
      return name_;
    }

    //! Return number of dimensions
    UInt GetDim( ) const { return dim_; }
    
    //! Transform local into global coordinate

    //! This method transforms a point given in local coordinates into a
    //! point w.r.t to the global cartesian x,y,z coordinates.
    //! \param glob (out) point w.r.t. to global cartesian coordinate system
    //! \param loc (in) point w.r.t. to local coordinate system
    virtual void Local2GlobalCoord( Vector<Double> & glob, 
                                    const Vector<Double> & loc ) const = 0;
    
    //! Transform global into local coordinate

    //! This method transforms a point given in global cartesian coordinates
    //! into a representation w.r.t to the local coordinate system.
    //! \param loc (out) point w.r.t. to local coordinate system
    //! \param glob (in) point w.r.t. to global cartesian coordinate system
    virtual void Global2LocalCoord( Vector<Double> &loc,  
                                    const Vector<Double> & glob ) const = 0;


    //! Return the global rotation angles for a given point

    //! This method returns the rotation angles about the x,y, and z axis,
    //! by which the global coordinate system has to be rotated, so that
    //! it represents the current one in that point.
    virtual void 
    GetGlobRotationAngles( Vector<Double> & angles,
                           const Vector<Double>& point ) const = 0;

    //@{
    //! Transform local vector into global one for a given global model point

    //! This method transforms a vector with a local coordinate representation
    //! and a given global model point (german: "Aufpunkt") into a vector with 
    //! a representation in global cartesian coordinates.
    //! This method can be used for example to get for a global element 
    //! mid-point and a given local load vector the global representation of
    //! the load vector.
    //! \param globVec (out) vector w.r.t. to global cartesian coordinate 
    //!                      system
    //! \param locVec (in) vector w.r.t. to local coordinate system
    //! \param globModelPoint (in) global model point (german: "Aufpunkt") of
    //!                            of the local/global vector
    virtual void 
    Local2GlobalVector( Vector<Double> & globVec, 
                        const Vector<Double> & locVec, 
                        const Vector<Double> & globModelPoint ) const = 0;

    virtual void 
    Local2GlobalVector( Vector<Complex> & globVec, 
                        const Vector<Complex> & locVec, 
                        const Vector<Double> & globModelPoint ) const = 0;
    //@}
    

    //! Returns for a given coordinate name the according index

    //! This method returns for a given coordinate name (x,y,z,rad,...)
    //! the according index in the local vector representation.
    //! \param dof (in) name of a coordinate component
    //! \return index of the coordinate component
    virtual UInt GetVecComponent( const std::string & dof ) const = 0;

    //! Returns for a given coordinate index the according name

    //! This method returns for a given coordinate index (1,2,3)
    //! the according name (x,y,z,rad,...).
    //! \param dof (in) index of the coordinate component
    //! \return name of the coordinate component
    virtual const std::string GetDofName( const UInt dof ) const = 0;

  protected:

    //! Helper function for obtaining a point vector
    void GetPoint(Vector<Double> & vec, 
                  ParamNode * pointNode );

    //! Calculate rotation matrix and inverse
    virtual void CalcRotationMatrix() = 0;
    
    //! Extract kardan angles from 3x3 rotation matrix
    void CalcKardanAngles( Vector<Double>& angles,
                           Matrix<Double>& rotMat );

    //! Return the correct angle (rad) for given sin and cos value
    Double GetAngle( Double sinAlpha, Double cosAlpha );

    //! Name of the coordinate system
    std::string name_;

    //! Dimension of the coordinate system
    UInt dim_;

    //! Origin of coordinate system
    Vector<Double> origin_;
    
    //! Rotation matrix (global -> local)
    Matrix<Double> rotationMat_;

    //! Rotation angles (kardan / Bryant angles) (global -> local)
    Vector<Double> rotationAng_;

    //! Inverse of the rotation matrix (local -> global)
    Matrix<Double> invRotationMat_;

    //! Rotation angles (kardan / Bryant angles) (local -> global)
    Vector<Double> invRotationAng_;
    
    //! Pointer to grid object
    Grid * ptGrid_;

    //! Pointer to parameter node defining the current coordinate system
    ParamNode * myParam_;
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class CoordSystem
  //! 
  //! \purpose 
  //! This class serves as base class for all representations of a local 
  //! coordinate system. In CFS++, all methods for integration and geometry
  //! handling expect by default a representation of points and vectors w.r.t.
  //! the \e global cartesian system with its x-,y- and z-axis. 
  //! In order to allow different representations of vectors and points, they
  //! can be described w.r.t. to a \e local coordinate system. Local coordinate
  //! systems are defined by their main axes and a point representing the origin
  //! of the local coordinate system.
  //! 
  //! \collab 
  //! Objects of this class get instantiated in the class Domain and are also
  //! administrated by it. They are used for example in the implenetation of
  //! SinglePDE::DefineIntegrators() to prescribe loads or forces.
  //!
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
