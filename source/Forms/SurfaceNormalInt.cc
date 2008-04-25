// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "SurfaceNormalInt.hh"
#include "Domain/surfElem.hh"

using namespace CoupledField;

#include "DataInOut/Logging/cfslog.hh"

DECLARE_LOG(forms)

SurfaceNormalInt::SurfaceNormalInt( BaseMaterial* matData, SubTensorType type) : BDBInt(matData,type) 
{
  volElem = NULL;
  name_ = "SurfaceNormalInt";
  subTensorType_ = type;

  nrDofs_ = 3;
  dimD_ = 12;
  
  nrDofsPerNode_ = 3;
}

SurfaceNormalInt::~SurfaceNormalInt() 
{
}

void SurfaceNormalInt::PostInit(Elem* elem)
{
  // this is a first runtime check!
  SurfElem* se = dynamic_cast<SurfElem*>(elem);
  
  // store the corresponding volume type base element
  volElem = se->ptVolElem1->ptElem;
  assert(volElem != NULL);
}


void SurfaceNormalInt::calcBMat( Matrix<Double> &bMat, UInt ip,
    Matrix<Double> &ptCoord ) 
{
  //ptelem->GetShFnc()
  
  // We are here a surface element but use the ansatz functions from
  // the corresponding volume element!
  const unsigned int numFncs  = ptelem->GetNumFncs( ansatzFct1_ );
  const unsigned int nrDofs   = getNrDofs();  

  // check the last element wie insert!
  assert(dimD_ = numFncs * nrDofs);
  
  bMat.Resize(nrDofs, dimD_);  
  bMat.Init();

  Vector<Double> shapeFncAtIp; // (N1 N2 ... Nn)^T where n = numFncs
  ptelem->GetShFncAtIp(shapeFncAtIp, ip, it1_.GetElem() );

  assert(shapeFncAtIp.GetSize() == numFncs);

  LOG_DBG3(forms) << "SurfaceNormalInt::calcBMat(): Ansatz functions are " << shapeFncAtIp.ToString();
  
  // make
  // N1 0  0  N2 0  0  N3 0
  // 0  N1 0  0  N2 0  0  N3
  // 0  0  N1 0  0  N2 0  0  N3 ...
  for(unsigned int b = 0; b < numFncs; b++) // block
    for(unsigned int d = 0; d < nrDofs; d++)
      bMat[d][b*nrDofs + d] = shapeFncAtIp[b];

  LOG_DBG3(forms) << "calcBMat at it " << ip << " -> " << bMat.ToString();
  //std::cout << bMat.ToString() << std::endl;
}



template <class T>
void SurfaceNormalInt::CalcNormalMatrix( Matrix<T> &out )
{
  const unsigned int size   = getNrDofs();

  out.Resize(size, size);
  out.Init();
  // todo! easy to find real normal!
  out[2][2] = 1.0;
}

