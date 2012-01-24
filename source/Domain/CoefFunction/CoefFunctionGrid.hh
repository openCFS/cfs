// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionGrid.hh 
 *       \brief    Coefficient function which obtains and interpolates values from another Grid
 *
 *       \date     01/23/2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef COEFFUNCTION_HH
#define COEFFUNCTION_HH

#include "CoefFunction.hh"

namespace CoupledField{

template<class DATA_TYPE>
class CoefFunctionGrid : public CoefFunction{
  public:

    CoefFunctionGrid();

    ~CoefFunctionGrid();

    // ========================
    //  ACCESS METHODS
    // ========================
    //@{ \name Access Methods

    virtual void GetTensor(Matrix<DATA_TYPE>& CoefMat,
                           const LocPointMapped& lpm );

    virtual void GetVector(Vector<DATA_TYPE>& CoefMat,
                           const LocPointMapped& lpm );

    virtual void GetScalar(DATA_TYPE& CoefMat,
                           const LocPointMapped& lpm );
    //@}

    //! Dump coefficient function to string
    virtual std::string ToString() const;

  protected:
  private:

};

/*! \class CoefFunctionGrid
 *     \brief Coefficient function which obtains and interpolates values from another Grid
 *     \tparam DATA_TYPE Can be Complex or Double
 *     @author A. Hueppe
 *     @date 01/2012
 *
 */

}
#endif
