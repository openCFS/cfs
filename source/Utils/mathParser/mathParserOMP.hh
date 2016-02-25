// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     mathParserOMP.hh
 *       \brief    <Description>
 *
 *       \date     Mar 3, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef MATHPARSEROMP_HH_
#define MATHPARSEROMP_HH_

#include "mathParser.hh"
#include "Utils/ThreadLocalStorage.hh"
#include <omp.h>

namespace CoupledField{


/* This path of making the Mathparser calls more thread safe is
 * chosen because it influences the rest of the code only insignificantly
 * on the downside, we end up with a somehow fragile implementation.
 * Basically the functionality of the standard parser is used and
 * only functions which might be called from parallel regions are
 * overwritten. Each thread gets its own handle but on the outside of
 * this class we only give one handle from thread 0. On call we map
 * internally to the real thread-local handle.
 * For future developments we might want to rethink the concept
 * and take the effort of altering every parser call in CFS
 */
//! Thread safe wrapper class for MathParser
class MathParserOMP : public MathParser{

public:
  //! constructor
  MathParserOMP();

  virtual ~MathParserOMP();

  virtual HandleType GetNewHandle( bool setDefaults = false);

  virtual void ReleaseHandle( HandleType handle );

  virtual void SetExpr( HandleType handle, const std::string &expr );

  virtual Double Eval( HandleType handle );

  virtual void EvalVector( HandleType handle, Vector<Double>& vec );

  virtual void EvalDivVector( HandleType handle, Double& divergence );

  virtual void EvalMatrix( HandleType handle, Matrix<Double>& matrix,
                       UInt numRows = 0, UInt numCols = 0 );

  virtual void Dump( std::ostream& os );

  virtual void SetValue( HandleType handle,
                     const std::string &varName,
                     Double val );

  virtual void RegisterExternalVar( HandleType handle,
                                const std::string& varName,
                                Double * ptVar );

  virtual void SetCoordinates( HandleType handle,
                                   const CoordSystem &coosy,
                                   const Vector<Double> &globCoord );

  virtual boost::signals2::connection
      AddExpChangeCallBack( const MathParserSignal::slot_function_type
                            &subscriber,
                            HandleType handle );
private:

  std::map<HandleType, CfsTLS<HandleType> > threadHandles_;

  bool isShared_;

};

}
#endif /* MATHPARSEROMP_HH_ */
