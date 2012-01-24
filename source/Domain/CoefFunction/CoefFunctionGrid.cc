// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGrid.cc 
 *       \brief    Implementation file for the Grid interpolation CoefFunction
 *
 *       \date     Jan. 23, 2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#include "CoefFunctionGrid.hh"

namespace CoupledField{

template<class DATA_TYPE>
CoefFunctionGrid<DATA_TYPE>::CoefFunctionGrid(){

}

template<class DATA_TYPE>
CoefFunctionGrid<DATA_TYPE>::~CoefFunctionGrid(){

}

template<class DATA_TYPE>
void CoefFunctionGrid<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                                            const LocPointMapped& lpm ){

}

template<class DATA_TYPE>
void CoefFunctionGrid<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                                            const LocPointMapped& lpm ){

}

template<class DATA_TYPE>
void CoefFunctionGrid<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                                           const LocPointMapped& lpm ){

}

template<class DATA_TYPE>
std::string CoefFunctionGrid<DATA_TYPE>::ToString(){

}


}
