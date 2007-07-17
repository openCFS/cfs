// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLINCURLCURLNODEINT
#define FILE_NLINCURLCURLNODEINT

#include "Utils/ApproxData.hh"
#include "baseForm.hh"
#include "curlCurlNodeInt.hh"

namespace CoupledField
{

  //! nonlinear curl-curl operator for magnetics 2D
  class nLinCurlCurlNode2DInt : public CurlCurlNode2DInt
  {
  public:


    //! constructor
    nLinCurlCurlNode2DInt(  BaseMaterial* matData, bool axi=false, 
                            bool coordUpdate = false );

    //! destructor 
    virtual ~nLinCurlCurlNode2DInt();

    //! Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Double>& stiffMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
    //! sets type of nonlinear algorithm
    virtual void SetNonLinMethod(std::string atype);

  protected: 

    Double ComputeDiffReluctivity( UInt nrEl, Vector<Double>& val );

  private:
    Double startmatVal_;       //!<  start value for reluctivity
    ApproxData *nlinFnc_;      //!< pointer to BH approximation object
    Vector<Double> magPot_;    //!< magnetic vector potential at nodes
    NonLinMethod nonLinType_;  //!< type of nonlinear algorithm
    bool isHysteresis_ ;       //!< magnetic hystersis is considered
  };


  //! nonlinear curl-curl operator for magnetics 3D
  class nLinCurlCurlNode3DInt : public CurlCurlNode3DInt
  {
  public:


    //! constructor
    nLinCurlCurlNode3DInt( BaseMaterial* matData, bool coordUpdate = false );

    //! destructor
    virtual ~nLinCurlCurlNode3DInt();

    //! Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Double>& stiffMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
    //! sets type of nonlinear algorithm
    virtual void SetNonLinMethod(std::string atype);

  protected: 

  private:

    Double startmatVal_;       //!<  start value for reluctivity
    ApproxData *nlinFnc_;      //!< pointer to BH approximation object
    Vector<Double> magPot_;    //!< magnetic vector potential at nodes
    NonLinMethod nonLinType_;  //!< type of nonlinear algorithm

  };

}

#endif // FILE_NLINCURLCURLNODEINT
