// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include "multigrid/smoother.hh"

/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_SMOOTHER
#define DEBUG_SMOOTHER
#endif // DEBUG_SMOOTHER
#define  debug  &std::cerr
#endif // DEBUG_TO_CERR
/**********************************************************/

namespace OLAS {
/**********************************************************/

template <typename T>
Smoother<T>::Smoother()
           : prepared_(false)
{
}

/**********************************************************/

template <typename T>
Smoother<T>::~Smoother()
{
}

/**********************************************************/

template <typename T>
void Smoother<T>::Reset()
{

    prepared_ = false;
}

/**********************************************************/

template <typename T>
bool Smoother<T>::Setup( const CRS_Matrix<T>& matrix )
{

    return prepared_ = true;
}

/**********************************************************/
} // namespace OLAS

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
