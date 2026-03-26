#ifndef COORD_SYSTEM_HH
#define COORD_SYSTEM_HH

#include <string>

#include "Domain/Mesh/Grid.hh"
#include "Utils/tools.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  // forward class declarations
  

  //! Base class for describing a local coordinate system
  class CoordSystem {

  public:

    //! Default constructor
    //! \param name (in) name of the coordinate system which is used to 
    //!                  find the additional parameters in the .xml-file
    //! \param ptGrid (in) pointer to finite element grid object
    //! \param myParamNode (in) pointer to parameter node of current coosy
    CoordSystem(const std::string & name, Grid * ptGrid,
                PtrParamNode myParamNode );

    //! Destructor
    virtual ~CoordSystem();
    
    //! Print information about coordinate system to info node
    virtual void ToInfo( PtrParamNode in ) {};
    
    //! Return the name of the coordinate system
    const std::string & GetName() const {
      return name_;
    }

    //! Return number of dimensions
    UInt GetDim( ) const { return dim_; }
    
    //! Return the origin of the coordinate system
    const Vector<Double>* GetOrigin( ) const { return &origin_; }

    //! Return the rotation matrix of the coordinate system
    const Matrix<Double>* GetRotationMatrix( ) const { return &rotationMat_; }

    //! Return the paramNode defining the coordinate system
    const PtrParamNode GetParamNode( ) const { return myParam_; }

    //! Transform local into global coordinate

    //! This method transforms a point given in local coordinates into a
    //! point w.r.t to the global cartesian x,y,z coordinates.
    //! \param glob (out) point w.r.t. to global cartesian coordinate system
    //! \param loc (in) point w.r.t. to local coordinate system
    virtual void Local2GlobalCoord( Vector<Double> & glob, 
                                    const Vector<Double> & loc ) const = 0;
    
    //! Transform global into local coordinate

    //! This method transforms a point given in global Cartesian coordinates
    //! into a representation w.r.t to the local coordinate system.
    //! \param loc (out) point w.r.t. to local coordinate system
    //! \param glob (in) point w.r.t. to global Cartesian coordinate system
    virtual void Global2LocalCoord( Vector<Double> &loc,  
                                    const Vector<Double> & glob ) const = 0;


    //! Return the global rotation matrix for a given point

    //! This method returns the rotation matrix defining defining a rotation,
    //! by which the global coordinate system has to be rotated, so that
    //! so that it represents the current one in that point.
    //! \param rotMatrix rotation matrix for global point
    //! \param point point w.r.t. to global Cartesian coordinate system
    virtual void 
    GetGlobRotationMatrix( Matrix<Double> & rotMatrix,
                           const Vector<Double>& point ) const = 0;
    
    //! Return the full 3x3 global rotation matrix for a given point

    //! This method returns the full 3x3 rotation matrix defining defining 
    //! a rotation, by which the global coordinate system has to be rotated, 
    //! it represents the current one in that point.
    //! \param rotMatrix rotation matrix for global point
    //! \param point point w.r.t. to global Cartesian coordinate system
    virtual void 
    GetFullGlobRotationMatrix( Matrix<Double> & rotMatrix,
                               const Vector<Double>& point ) const = 0;
    
    //! Return the full 3x3 global rotation matrix for a given point

    //! This method returns the full 3x3 rotation matrix defining
    //! a rotation, by which the global coordinate system has to be rotated,
    //! it represents the current one in that point.
    //! \param rotMatrix rotation matrix for global point
    //! \param lpm local point w.r.t. to an element
    virtual void GetFullGlobRotationMatrix(Matrix<Double> & rotMatrix,
                                           const LocPointMapped &lpm) const;

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
    

    
    static constexpr unsigned int INVALID_DOF = 123456; 
    //! This method returns for a given coordinate name (x,y,z,rad,...)
    //! the according index in the local vector representation.
    //! \param dof  name of a coordinate component
    //! \param silent no exception is thrown for invalid dof but INVALID_DOF is returned
    //! \return index of the coordinate component or INVALID_DOF
    virtual UInt GetVecComponent( const std::string & dof, bool silent = false ) const = 0;


    //! Returns for a given coordinate index the according name

    //! This method returns for a given coordinate index (1,2,3)
    //! the according name (x,y,z,rad,...).
    //! \param dof (in) index of the coordinate component
    //! \return name of the coordinate component
    virtual const std::string GetDofName( const UInt dof ) const = 0;

  protected:

    //! Helper function for obtaining a point vector
    void GetPoint(Vector<Double> & vec, 
                  PtrParamNode pointNode );

    //! Calculate rotation matrix and inverse
    virtual void CalcRotationMatrix() = 0;
    
    //! Extract kardan angles from 3x3 rotation matrix
    void CalcKardanAngles( Vector<Double>& angles,
                           Matrix<Double>& rotMat );

    //! Return the correct angle (rad) for given sin and cos value
    Double GetAngle( Double sinAlpha, Double cosAlpha );

    //! Check, if initial definition of the rotation matrix is correct
    void CheckRotationMat(const Matrix<Double>& rotMat );
    
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
    
    //! Full 3x3 inverse of the rotation matrix (local -> global)
    Matrix<Double> invRotationMatFull_;
    
    //! Rotation angles (kardan / Bryant angles) (local -> global)
    Vector<Double> invRotationAng_;

    //! Pointer to grid object
    Grid * ptGrid_;

    //! Pointer to parameter node defining the current coordinate system
    PtrParamNode myParam_;
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class CoordSystem
  //! 
  //! \purpose 
  //! This class serves as base class for all representations of a local 
  //! coordinate system. In openCFS, all methods for integration and geometry
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
