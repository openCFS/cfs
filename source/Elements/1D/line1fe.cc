// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "line1fe.hh"
#include "Domain/elem.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField
{

  // declare class specific logging stream
  DECLARE_LOG(line1fe)
  DEFINE_LOG(line1fe,"line1fe")

  Line1FE::Line1FE(IntegMethod method, int order) : LineFE()
  {

    Init(method, order);
  }

  Line1FE::~Line1FE()
  {
    if(sDerivAtIp_) delete[] sDerivAtIp_;
    if(sShFcnAtIp_) delete[] sShFcnAtIp_;
  }

  void Line1FE::Init(IntegMethod method, int order)
  {
    NumNodes_ = 2;

    CommonInit(method, order);
    SetEdgeIndices();
    sDerivAtIp_ = new Vector<Double>[1];
    sShFcnAtIp_ = new Vector<Double>[1];
  }

  void Line1FE::SetCornerCoords()
  {

    LCornerCoords_.Resize(Dim_,NumNodes_);

    LCornerCoords_[0][0] =   1;
    LCornerCoords_[0][1] =  -1;

  }

  void Line1FE::SetEdgeIndices () {

    edgeIndices_ = new StdVector<UInt>[NumEdges_];
    edgeIndices_[0].Resize(2);

    // edge 1
    edgeIndices_[0][0] = 1;
    edgeIndices_[0][1] = 2;

  }


  void Line1FE :: CalcShapeFnc(Vector<Double> & Shape,
                               const Vector<Double> & LCoord,
                               const Elem* elem, UInt dof,
                               AnsatzFct::FctEntityType type)
  {

    if(  actFct_->GetType() == AnsatzFct::LAGRANGE ||
         type == AnsatzFct::NODE ) {
      Shape.Resize(NumNodes_);

      for( UInt i=0; i<NumNodes_; i++)
        Shape[i] = 0.5*(1+LCornerCoords_[0][i] * LCoord[0]);
    }else if(  actFct_->GetType() == AnsatzFct::SPECTRAL) {
      // ===============
      //  Spectral PART
      // ===============
      CalcSpectralShFct(Shape,LCoord,elem,dof,type);

    }else{
      // Get number of ansatz functions
      shared_ptr<LegendreFct> legFct =
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_);

      // ===============
      //  LEGENDRE PART
      // ===============

      // Get number of ansatz functions
      UInt totalFcns = GetNumFncs( actFct_ );

      Shape.Resize( totalFcns );

      // --------------------
      //  a) nodal functions
      // --------------------
         // Offset for different functions
      UInt offset = 0;

      // First of all, calculate all nodal function derivatives
      for( UInt i=0; i<NumNodes_; i++)
        Shape[i] = 0.5*(1+LCornerCoords_[0][i] * LCoord[0]);

      offset = NumNodes_;
      // --------------------
      //  b) edge functions
      // --------------------
      // Obtain order of element
      Integer order = dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_)->GetIsoOrder();
      Double val, factor, deriv;

      factor = elem->edges[0] < 0 ? -1.0 : 1.0;
      for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {
        EvalPolynom( val, deriv, iDof, lCoeff_[iDof],
                     factor*LCoord[0] );
        Shape[offset] = val;
      }
    }


  }

  void Line1FE::CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv,
                                         const Vector<Double> & LCoord,
                                         const Elem* elem, UInt dof,
                                         AnsatzFct::FctEntityType type)
  {
    if(  actFct_->GetType() == AnsatzFct::LAGRANGE ||
         type == AnsatzFct::NODE ) {
      LDeriv.Resize(NumNodes_,Dim_);

      for( UInt i=0; i<NumNodes_; i++) {
          LDeriv[i][0] = 0.5*LCornerCoords_[0][i];
        }
      return;

    }else if(  actFct_->GetType() == AnsatzFct::SPECTRAL) {
      // ===============
      //  Spectral PART
      // ===============
      CalcSpectralDerivFct(LDeriv,LCoord,elem,dof,type);
      return;

    }else{
      // ===============
      //  LEGENDRE PART
      // ===============
      shared_ptr<LegendreFct> legFct =
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_);

      // Get number of ansatz functions
      UInt totalFcns = GetNumFncs( actFct_ );
      LDeriv.Resize( totalFcns, Dim_ );

      // --------------------
      //  a) nodal functions
      // --------------------
      UInt offset = 0;

      // First of all, calculate all nodal function derivatives
      for( UInt i=0; i<NumNodes_; i++) {
          LDeriv[i][0] = 0.5*LCornerCoords_[0][i];
        }

      // -------------------
      //  b) edge functions
      // -------------------

      offset = NumNodes_;
      Double val, deriv, factor;
      // Obtain order of element
      Integer order = dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct_)->GetIsoOrder();

      factor = elem->edges[0] < 0 ? -1.0 : 1.0;
      for( Integer iDof = 2; iDof <= order; iDof++, offset++ ) {
        EvalPolynom( val, deriv, iDof, lCoeff_[iDof],
                     factor*LCoord[0] );
        LDeriv[offset][0] =  deriv * factor;
      }
    }
  }

  void Line1FE::Global2LocalCoords(Matrix<Double> & localCoords,
                                     const Matrix<Double> & globalCoords,
                                     const Matrix<Double> & coordMat)
  {
    Vector<Double> c0, c1; // endpoint-coordinates
    Vector<Double> diff0, diff1;
    Vector<Double> dummy;
    Double s;
    Double dist, dist1, fac;
    UInt globDim = globalCoords.GetNumRows();

    // Get coordinates of the endpoints
    c0.Resize(globDim);
    //    c1.Resize(Dim_);
    //    dummy.Resize(Dim_);
    diff0.Resize(globDim);
    diff1.Resize(globDim);

    //    std::cout << "SIMON> Line1FE->Global2Local(): coordMat " << coordMat << std::endl;
    //    std::cout << "SIMON> Line1FE->Global2Local(): globalCoords " << globalCoords << std::endl;

    for(UInt i = 0; i < globDim; i++)
    {
        c0[i] = coordMat[i][0];
        //        c1[i] = coordMat[i][1];
        diff1[i] = coordMat[i][1] - c0[i];
    }


    //    diff1 = c1 - c0;
    dist = diff1.NormL2();
    fac = 1.0 / dist;
    diff1 *= fac;

    localCoords.Resize(Dim_, globalCoords.GetNumCols());

    for(UInt i=0; i < globalCoords.GetNumCols(); i++)
    {
        for(UInt j = 0; j < globDim; j++)
        {
            diff0[j] = globalCoords[j][i] - c0[j];
        }
        //        diff0 = globalCoords[i] - c0;
        dist1 = diff0.NormL2();
        diff0.Inner(diff1, s);
        s *= fac;

        //        std::cout << "SIMON> Line1FE->Global2Local(): s " << s << std::endl;

        for(UInt j = 0; j < Dim_; j++)
        {
            localCoords[j][i] = LCornerCoords_[j][0] * s + LCornerCoords_[j][1] * (1-s);
            //            std::cout << "SIMON> Line1FE->Global2Local(): localcoord " <<localCoords[j][i]<< std::endl;
        }
    }
  }

  void  Line1FE::GetNumFncs( StdVector<UInt>& numFcns,
                             const shared_ptr<AnsatzFct>& fcnType,
                             AnsatzFct::FctEntityType fctEntityType,
                             UInt dof ) {

    // Check ansatzFctType
    if( fcnType->GetType() == AnsatzFct::LAGRANGE ) {
      numFcns.Resize( NumNodes_ );
      numFcns.Init(1);
    }else if(fcnType->GetType() == AnsatzFct::SPECTRAL){

      // Remember approximation order
      Integer order =  dynamic_pointer_cast<SpectralFct, AnsatzFct> (fcnType)->GetOrder();

      if( fctEntityType == AnsatzFct::NODE ) {
        numFcns.Resize( NumNodes_ );
        numFcns.Init( 1 );

      } else if( fctEntityType == AnsatzFct::EDGE ) {
        numFcns.Resize( NumEdges_ );
        numFcns.Init( order - 1 );
      }
      else if( fctEntityType == AnsatzFct::INTERIOR ) {
        numFcns.Resize(1);
        numFcns.Init(0);
      }
      else if( fctEntityType == AnsatzFct::FACE ) {
        numFcns.Resize(1);
        numFcns.Init(0);
      } else {
        EXCEPTION( "Not yet implemented!" );
      }

    } else if ( fcnType->GetType() == AnsatzFct::LEGENDRE ) {

      // Remember approximation order
      Integer order =  dynamic_pointer_cast<LegendreFct, AnsatzFct>
        (fcnType)->GetIsoOrder();

      if( fcnType->IsIsotropic() == true ) {
        // Check for subentity-type
        if( fctEntityType == AnsatzFct::NODE ) {
          numFcns.Resize( NumNodes_ );
          numFcns.Init( 1 );

        } else if( fctEntityType == AnsatzFct::EDGE ) {
          numFcns.Resize( NumEdges_ );
          numFcns.Init( order - 1 );
        }
        else if( fctEntityType == AnsatzFct::INTERIOR ) {
          numFcns.Resize(1);
          numFcns.Init(0);
        } else if( fctEntityType == AnsatzFct::FACE ) {
          numFcns.Resize(1);
          numFcns.Init(0);
        } else {
          EXCEPTION( "Not yet implemented!" );
        }
      }else {
         // *** anisotropic case ***
        // Remember approximation order
        shared_ptr<LegendreFct> const & legFct =
          dynamic_pointer_cast<LegendreFct, AnsatzFct> (fcnType);
        Matrix<UInt> const & order =  legFct->GetAnisoOrder();

        // Check for subentity-type
        if( fctEntityType == AnsatzFct::NODE ) {
          numFcns.Resize( NumNodes_ );
          numFcns.Init( 1 );

        } else if( fctEntityType == AnsatzFct::EDGE ) {
          numFcns.Resize( NumEdges_ );
          numFcns.Init(0);
          // xi-direction
          numFcns[0]  = order[0][dof]-1;
        } else if( fctEntityType == AnsatzFct::FACE ) {
          numFcns.Resize(1);
          numFcns.Init(0);
        }else if( fctEntityType == AnsatzFct::INTERIOR ) {
          numFcns.Resize(1);
          numFcns.Init(0);
        }

      }
    } else {
      EXCEPTION( "AnsatzFcnType '" << fcnType->GetType()
                 << "' is not known!" );
    }
  }

  UInt Line1FE::GetNumFncs( const shared_ptr<AnsatzFct>& type ) {

    // TODO: FOr anisotropic functions we have
    // to incorporate the dof in the determination
    // Check ansatzFctType
    if( type->GetType() == AnsatzFct::LAGRANGE ) {
      return NumNodes_;

    }else if(type->GetType() == AnsatzFct::SPECTRAL){
      Integer order =  dynamic_pointer_cast<SpectralFct, AnsatzFct> (type)->GetOrder();
      return (order+1);

    } else if ( type->GetType() == AnsatzFct::LEGENDRE) {
      if( type->IsIsotropic() == true ) {

      // Remember approximation order
      Integer order =  dynamic_pointer_cast<LegendreFct, AnsatzFct>
        (type)->GetIsoOrder();

        // edge functions
        UInt numEdgeFncs = NumEdges_* (order-1);

        LOG_DBG3(line1fe) << "total number of unknowns: "
                          <<  NumNodes_ + numEdgeFncs
                          << std::endl;

        actNumFcns_ = (NumNodes_ + numEdgeFncs);
        return actNumFcns_;

      } else {
        // *** anisotropic case ***
        Integer max_1;
        shared_ptr<LegendreFct> legFct =
          dynamic_pointer_cast<LegendreFct, AnsatzFct>(type);

        max_1 = legFct->GetMaxOrderLocDir( 0 );

        // a) nodes
        UInt numNodeModes = 2;

        // b) edges
        UInt numEdgeModes = 0;
        numEdgeModes += max_1 > 0 ? (max_1-1)*2 : 0;

        actNumFcns_ = numNodeModes + numEdgeModes ;
        return actNumFcns_;

      }
    }
    return 0;
  }

  void Line1FE::SetAnsatzFct( shared_ptr<AnsatzFct>& actFct,
                              bool setIntPoints ) {

    // Check if this ansatz fct was already set
    if( actFct_ == actFct ) {
      return;
    }
    actFct_ = actFct;

    if( actFct->GetType() == AnsatzFct::SPECTRAL ) {
      // If not, get order of functions
      shared_ptr<SpectralFct> sFct = dynamic_pointer_cast<SpectralFct, AnsatzFct>(actFct);
      IntegMethod = LOBATTO;
      IntegOrder = 2* sFct->GetOrder() - 1;
      sShFcnAtIp_[0].Resize(sFct->GetOrder()+1);
      sDerivAtIp_[0].Resize(sFct->GetOrder()+1);

      // get the values by IntegMethod and IntegOrder
      SetIntPoints();

      // ... then calc shape function values at integration points
      // for subsequent calls to CalcJacobian() ....
      SetShapeFncAtIp();
      SetShapeFncDerivAtIp();

    }else if( actFct->GetType() == AnsatzFct::LEGENDRE
              && setIntPoints == true ) {
      // If not, get order of functions
      shared_ptr<LegendreFct> legFct =
        dynamic_pointer_cast<LegendreFct, AnsatzFct>(actFct);
      // check if isotropic order is present
      if( legFct->IsIsotropic() ) {
        // Prevent integration of first order, as this may
        // cause non-reasonable results
        //if( legFct->GetIsoOrder() > 1 ) {
        IntegMethod = CARTESIAN;
        //IntegOrder = legFct->GetIsoOrder() *2;
        IntegOrder =  EncodeCartesianOrder( legFct->GetIsoOrder() *2,
                                            legFct->GetIsoOrder() *2,
                                            legFct->GetIsoOrder() *2);
      } else {
        IntegMethod = CARTESIAN;
        IntegOrder =  EncodeCartesianOrder( legFct->GetMaxOrder() * 2,
                                            legFct->GetMaxOrder() * 2,
                                            legFct->GetMaxOrder() * 2);
      }
      // get the values by IntegMethod and IntegOrder
      SetIntPoints();
      // ... then calc shape function values at integration points
      // for subsequent calls to CalcJacobian() ....
      SetShapeFncAtIp();
      SetShapeFncDerivAtIp();
    }
  }

  void Line1FE::CalcSpectralShFct(Vector<Double> & Shape,
                                  const Vector<Double> & LCoord,
                                  const Elem* elem, UInt dof,
                                  AnsatzFct::FctEntityType type ){
      shared_ptr<SpectralFct> myFct = dynamic_pointer_cast<SpectralFct, AnsatzFct>(actFct_);
      UInt order = myFct->GetOrder();

      Shape.Resize( (order+1) );
      Shape.Init();
      //now get the shape functions and the derivatives for the given coordinates
      sShFcnAtIp_[0].Init();
      myFct->EvaluatePolynomial( sShFcnAtIp_[0], LCoord[0] );

      UInt counter = 0;
      Shape[counter++] = sShFcnAtIp_[0][0];
      Shape[counter++] = sShFcnAtIp_[0][order];

      Integer factor = elem->edges[0];
      if(factor < 0){
        for ( UInt i= 1 ; i< order  ;i++ )
          Shape[counter++] = sShFcnAtIp_[0][order - i];
      }else{
        for ( UInt i= 1; i< order ;i++ )
          Shape[counter++] = sShFcnAtIp_[0][i];
      }
      return;
  }

  void Line1FE::CalcSpectralDerivFct( Matrix<Double> & LDeriv,
                                       const Vector<Double> & LCoord,
                                       const Elem* elem, UInt dof,
                                       AnsatzFct::FctEntityType type){
      // ===============
      //  SPECTRAL PART
      // ===============

      //Check if we are only interested in the nodal values
      //now do the calculations for the spectral fem
      shared_ptr<SpectralFct> myFct = dynamic_pointer_cast<SpectralFct, AnsatzFct>(actFct_);
      UInt order = myFct->GetOrder();

      sDerivAtIp_[0].Init();

      LDeriv.Resize( (order+1), Dim_ );
      //now get the shape functions and the derivatives for the given coordinates
      myFct->EvaluateDerivPolynomial( sDerivAtIp_[0],LCoord[0] );
      LDeriv[0][0] = sDerivAtIp_[0][0];
      LDeriv[1][0] = sDerivAtIp_[0][order];

      LDeriv.Resize(NumNodes_,Dim_);
      UInt counter = 2;
      if(elem->edges[0] < 0){
        for ( UInt i=1;i<order;i++ )
          LDeriv[counter++][0] = sDerivAtIp_[0][order - i];
      }else{
        for ( UInt i=1;i<order ;i++ )
          LDeriv[counter++][0] = sDerivAtIp_[0][i];
      }
      return;

    }

} // end of namespace
