// ================================================================================================
/*!
 *       \file     SMSM.hh
 *       \brief    Simplified Multiscale Model (without texture yet)
 *
 *       \date     June 12, 2024
 *       \author   Klaus Roppert, TU Graz
 */
//================================================================================================

#pragma once
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include <list>

#include <boost/utility.hpp>



namespace CoupledField {


  class SMSM {

    public:
      //! Constructor
      SMSM();

      //! Destructor
      virtual ~SMSM();

      void Eval(Double valH, StdVector<Double> dirHloc);

      Vector<Double> GetM(){return MMoy_;};
      Matrix<Double> GetEps(){return epsmumoy_;};
      

    private:

      StdVector<Vector<Double>> TABgamma_;
      Vector<Double> epsmu11_, epsmu12_, epsmu13_, epsmu22_, epsmu23_, epsmu33_, Wan_;
      Double l100_, l111_, AS_, K1_, K2_, Ms_, mu0_;
      Matrix<Double> SIGMAloc_;
      UInt numRows_;
      Vector<Double> MMoy_;
      Matrix<Double> epsmumoy_;
    };
} //end of namespace
