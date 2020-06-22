/*
 * registerfunc.cc
 *
 *  Common definition of math functions for MathParser and PyMuParser
 *
 *  Created on: 12.06.2020
 *      Author: Jens Grabinger <jensemann126@gmail.com>
 */

#include "registerfunc.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/Interpolate1D.hh"

namespace CoupledField {

// Static Aliases For Operators <,<=, >, >=
static Double Op_gt( Double a, Double b ) { return a > b; }
static Double Op_ge( Double a, Double b ) { return a >= b; }
static Double Op_le( Double a, Double b ) { return a <= b; }
static Double Op_lt( Double a, Double b ) { return a < b; }

// Legacy if-function
static Double IfFunc( Double condition, Double thenVal, Double elseVal ) {
  if (condition > 0.0 )
    return thenVal;
  else
    return elseVal;
}

// Register our functions and operators with a muParser instance
void RegisterFunctions(mu::Parser &parser) {

  // Register alias functions for >, >=, <=, <
  parser.DefineOprt( "gt", Op_gt, 2);
  parser.DefineOprt( "ge", Op_ge, 2);
  parser.DefineOprt( "le", Op_le, 2);
  parser.DefineOprt( "lt", Op_lt, 2);

  // Register constant variables
  parser.DefineConst("pi", (double) M_PI);

  // Register if-function for legacy support
  parser.DefineFun("if", IfFunc, false );

  // Register functions from within CFS
  parser.DefineFun("sample1D", Interpolate1D::Interpolate, false );

  // Register signal generating functions
  parser.DefineFun("sinBurst", SinBurst, false );
  parser.DefineFun("fadeIn", FadeIn, false );
  parser.DefineFun("spike", Spike, false );
  parser.DefineFun("chirp", Chirp, false );
  parser.DefineFun("cosPulseComb", CosPulseComb, false );
  parser.DefineFun("squareBurst", SquarePulse, false );
  parser.DefineFun("gauss", Gauss, false );
  parser.DefineFun("triangle", Triangle, false );

  // Register general functions
  parser.DefineFun("mod", Mod, false );
  parser.DefineFun("besselCylJ", BesselCylJ, false );
  parser.DefineFun("besselCylY", BesselCylY, false );
  parser.DefineFun("besselSphJ", BesselSphJ, false );
  parser.DefineFun("besselSphY", BesselSphY, false );
}

} // end of namespace
