#ifndef TRIVIAL_CARTESIAN_COORD_SYSTEM_HH
#define TRIVIAL_CARTESIAN_COORD_SYSTEM_HH

#include "CoordSystem.hh"

namespace CoupledField {


  //! Describes a simplistic Cartesian coordinate system.
  class TrivialCartesianCoordSystem : public CoordSystem{

  public:

    //! Default constructor
    //! \param name (in) name of the coordinate system which is used to 
    //!                  find the additional parameters in the .xml-file
    //! \param ptGrid (in) pointer to finite element grid object
    //! \param myParamNode (in) pointer to parameter node of current coosy
    TrivialCartesianCoordSystem(const std::string & name, Grid * ptGrid,
                                PtrParamNode myParamNode );
    
    //! Destructor
    virtual ~TrivialCartesianCoordSystem();

    //! Print information about coordinate system to info node
    void ToInfo( PtrParamNode in );
    
    //! Transform local into global coordinate

    //! This method transforms a point given in local coordinates into a
    //! point w.r.t to the global cartesian x,y,z coordinates.
    //! \param glob (out) point w.r.t. to global cartesian coordinate system
    //! \param loc (in) point w.r.t. to local coordinate system
    void Local2GlobalCoord( Vector<Double> & glob, 
                            const Vector<Double> & loc ) const;
    
    //! Transform global into local coordinate

    //! This method transforms a point given in global cartesian coordinates
    //! into a representation w.r.t to the local coordinate system.
    //! \param loc (out) point w.r.t. to local coordinate system
    //! \param glob (in) point w.r.t. to global cartesian coordinate system
    void Global2LocalCoord( Vector<Double> &loc,  
                            const Vector<Double> & glob ) const;

    //! Returns the global rotation matrix for a given point

    //! This method returns the rotation matrix defining defining a rotation,
    //! by which the global coordinate system has to be rotated, so that
    //! so that it represents the current one in that point.
    //! \param rotMatrix rotation matrix for global point
    //! \param point point w.r.t. to global Cartesian coordinate system
    void GetGlobRotationMatrix(Matrix<Double>& rotMat, const Vector<Double>& /*point*/) const {
       rotMat = invRotationMat_;
     }

    //! Returns the full 3x3 global rotation matrix for a given point

    //! This method returns the full 3x3 rotation matrix defining defining 
    //! a rotation, by which the global coordinate system has to be rotated, 
    //! it represents the current one in that point.
    //! \param rotMatrix rotation matrix for global point
    //! \param point point w.r.t. to global Cartesian coordinate system
    void GetFullGlobRotationMatrix(Matrix<Double>& rotMat, const Vector<Double>& /*point*/) const {
       rotMat = invRotationMatFull_;
     }

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
    void Local2GlobalVector( Vector<Double> & globVec, 
                             const Vector<Double> & locVec, 
                             const Vector<Double> & globModelPoint ) const;

    void Local2GlobalVector( Vector<Complex> & globVec, 
                             const Vector<Complex> & locVec, 
                             const Vector<Double> & globModelPoint ) const;
    //@}
    
    //! Returns for a given coordinate name the according index

    //! This method returns for a given coordinate name (r,phi,z)
    //! the according index in the local vector representation.
    //! \param dof (in) name of a coordinate component
    //! \return index of the coordinate component
    UInt GetVecComponent( const std::string & dof ) const;

    //! Returns for a given coordinate index the according name

    //! This method returns for a given coordinate index (1,2,3)
    //! the according name (x,y,z).
    //! \param dof (in) index of the coordinate component
    //! \return name of the coordinate component
    const std::string GetDofName( const UInt dof ) const;


  protected:

    //! Calculate rotation matrix
    void CalcRotationMatrix();
    
    //! Internal implementation of Local2GlobalVectorInt
    template<class TYPE>
    void Local2GlobalVectorInt( Vector<TYPE> & globVec, 
                                const Vector<TYPE> & locVec, 
                                const Vector<Double> & globModelPoint ) const;
    
    //! global vector pointing in local x-direction
    Vector<Double> axisFactors_;

    //! global vector pointing in local y-direction
    Vector<UInt> axisMap_;

    //! true when origin is zero and axes are identity — enables zero-cost coord transform
    bool isIdentity_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class TrivialCartesianCoordSystem
  //! 
  //! \purpose 
  //! Thic class implements the representation of cartesian coordinate
  //! system. It is defined by an point of origin, and a mapping of how
  //! global axes map to local axes.
  //! 
  //! \collab 
  //! Ref. to base class CoordSystem.
  //!
  //! \implement 
  //! 
  //! \status not tested
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace

#endif
