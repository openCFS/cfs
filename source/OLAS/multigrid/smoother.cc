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
    ENTER_FCN("Smoother::Smoother");
}

/**********************************************************/

template <typename T>
Smoother<T>::~Smoother()
{
    ENTER_FCN("Smoother::~Smoother");
}

/**********************************************************/

template <typename T>
void Smoother<T>::Reset()
{
    ENTER_FCN("Smoother::Reset");

    prepared_ = false;
}

/**********************************************************/

template <typename T>
bool Smoother<T>::Setup( const CRS_Matrix<T>& matrix )
{
    ENTER_FCN("Smoother::Setup");

    return prepared_ = true;
}

/**********************************************************/
} // namespace OLAS

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
