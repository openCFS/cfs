// ================================================================================================
/*!
 *       \file     SMSM.hh
 *       \brief    Simplified Multiscale Model (without texture yet)
 * 
 *                 Integration into the Energy-Based (EB) vector hysteresis framework described in:
 *                   K. Roppert, M. Kaltenbacher, L. Domenig, and L. Daniel,
 *                   "Magneto-Elastic Vector Hysteresis Modelling for Electromagnetic Devices:
 *                   A Combination of a Multiscale Model with the Energy-Based Hysteresis Framework,"
 *                   IEEE Trans. Magn., 2025.
 *                 10.1109/TMAG.2025.3560666
 *                 https://ieeexplore.ieee.org/document/11062564
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
      SMSM(Double Ms, Double AS, Double K1, Double K2, Double lambda100, Double lambda111, UInt dim);

      //! Destructor
      virtual ~SMSM();

      void Eval(Double valH, StdVector<Double> dirHloc);
      void Eval3D(Double valH, StdVector<Double> dirHloc);
      void Eval2D(Double valH, StdVector<Double> dirHloc);

      void Register_stress(Vector<Double> sigma);

      Vector<Double> GetM(){return MMoy_;};
      Matrix<Double> GetEps(){return epsmumoy_;};
      Matrix<Double> GetdMdH(){return dMdH_;};
      Matrix<Double> GetSigma(){return SIGMAloc_;};

    private:
      
      // Discrete set of unit magnetization direction vectors gamma_i (numRows_ x dim_).
      // 3D: read from lookup table TABSPHEREI4S4_2562.txt (2562 points on the unit sphere).
      // 2D: uniformly distributed angles on the unit circle (360*3 points).
      StdVector<Vector<Double>> TABgamma_;

      // omputed magnetostrictive strain tensor components for each direction gamma_i
      Vector<Double> epsmu11_, epsmu12_, epsmu13_, epsmu22_, epsmu23_, epsmu33_;

      // computed magnetization components for each direction gamma_i
      Vector<Double> Wan_;

      // magnetostrictive constants
      Double l100_, l111_;
      
      // statistical parameter of the Boltzmann distribution
      Double AS_; 

      // anisotropy constants
      Double K1_, K2_;

      // saturation magnetization and vacuum permeability 
      Double Ms_, mu0_;

      // local mechanical (Cauchy) stress tensor (set by Register_stress(), used in Eval())
      Matrix<Double> SIGMAloc_;

      // number of discrete magnetization directions (rows of TABgamma_)
      UInt numRows_;

      // volume-averaged magnetization vector
      Vector<Double> MMoy_;

      // volume-averaged magnetostrictive strain tensor
      Matrix<Double> epsmumoy_;

      // differential susceptibility tensor dM/dH
      Matrix<Double> dMdH_;

      // spatial dimension (2 or 3)
      UInt dim_;
    };
} //end of namespace
