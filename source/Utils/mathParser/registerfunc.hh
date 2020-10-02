/*
 * registerfunc.hh
 *
 *  Common definition of math functions for MathParser and PyMuParser
 *
 *  Created on: 12.06.2020
 *      Author: Jens Grabinger <jensemann126@gmail.com>
 */

#ifndef SOURCE_UTILS_MATHPARSER_REGISTERFUNC_HH_
#define SOURCE_UTILS_MATHPARSER_REGISTERFUNC_HH_

#include "muParser.h"

namespace CoupledField {

//! Register our functions and operators with a muParser instance
void RegisterFunctions(mu::Parser &parser);

}

#endif /* SOURCE_UTILS_MATHPARSER_REGISTERFUNC_HH_ */
