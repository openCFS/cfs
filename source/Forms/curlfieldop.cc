#include "Forms/curlfieldop.hh"

#include <string>

#include "Elements/basefe.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "General/environment.hh"
#include "Utils/vector.hh"
#include "Matrix/matrix.hh"

#include <PDE/StdPDE.hh>

namespace CoupledField
{

  CurlEdgeOp::CurlEdgeOp(Grid * ptGrid, 
                         StdPDE * ptPDE,
                         shared_ptr<EqnMap> eqnMap,
                         NodeStoreSol<Double> & aSol,
                         BaseSystem * algsys, bool coordUpdate) 
    : BaseOperator(ptGrid, ptPDE, eqnMap, coordUpdate ), algsys_(algsys)
  {
    ENTER_FCN( "CurlEdgeOp::CurlEdgeOp" );

    sol_ = &aSol;
  }


  CurlEdgeOp::~CurlEdgeOp()
  {
    ENTER_FCN( "CurlEdgeOp::~CurlEdgeOp" );
  }


  void CurlEdgeOp::CalcElemCurlEdge(Vector<Double> & curlField, 
                                    const EntityIterator& ent,
                                    const Vector<Double> & LCoord)
  {
    ENTER_FCN( "CurlEdgeOp::CalcElemCurlEdgeOp" );
  
    Error("CurlEdgeOp::CalcElemCurlEdgeOp: Not working due to EQN-class!",
          __FILE__, __LINE__);
    //  ShortInt dim = ptElement->ptElem->GetDim();
  
    //   curlField.Resize(dim);
    //   for (UInt i=0; i<dim; i++)
    //     curlField[i] = 0;

    //   UInt nrEdges = ptElement->ptElem->GetNumEdges();
    //   UInt nrNodes = ptElement->ptElem->GetNumNodes();

    //   Matrix<Double> CornerCoords; 
    //   ptGrid_->GetCoordNodesElemMat(ptElement->connect, CornerCoords, level_);


    //   Matrix<Double> curlOnEdges;


    //   StdVector<Matrix<Double>* > deriv;
    //   deriv.Resize(nrEdges);
    //   for (UInt actEdge=0; actEdge < nrEdges; actEdge++)
    //     deriv[actEdge] = new Matrix<Double>;
  

    //   ptElement->ptElem->GetEdgeGlobalDerivShapeFnc (deriv, LCoord, CornerCoords);

    //   curlOnEdges.Resize(dim, nrEdges);
  
    //   for (UInt actEdge=0; actEdge < nrEdges; actEdge++)
    //     for (UInt actDim=0; actDim < dim; actDim++)
    //       curlOnEdges[actDim][actEdge] = 
    //      (*deriv[actEdge])[(actDim+2)%dim][(actDim+1)%dim] -
    //      (*deriv[actEdge])[(actDim+1)%dim][(actDim+2)%dim];


    //   StdVector<Double> sol(nrEdges);
    //   // global edge index
    //   StdVector<UInt> epos(nrEdges);
    //   StdVector<UInt> esign(nrEdges);


    //   StdVector<UInt> pos(nrNodes);
  
    //   for (UInt i=0; i < nrNodes; i++)
    //      pos[i] = (*ptMesh2PDENode_)[ptElement->connect[i]-1];
  
    //   ptElement->ptElem->GetGlobalEdgeIndices(epos, &pos[0], algsys_);

    // #ifdef DEBUG
    //   (*debug) << "CurlOP pos \n" << pos << std::endl
    //         << "epos \n " << epos << std::endl;
  
    // #endif

    //   for (UInt j=0; j<nrEdges; j++)
    //     {
    //       esign[j] = epos[j]/abs(epos[j]);
    //       epos[j]  = abs(epos[j]);
    //       sol[j] = (*sol_)[epos[j]-1] * esign[j];
    //     }
  
  
    //   // loop over edge curls
    //   for( UInt i=0; i<dim; i++ )
    //     {
    //       curlField[i]=0;
      
    //       for( UInt j=0; j < nrEdges; j++ )
    //      curlField[i] += curlOnEdges[i][j] * sol[j];
    //     }
  
  }




  void CurlEdgeOp::CalcElemMagVec(Vector<Double> & magVecPot, 
                                  const EntityIterator& ent,
                                  const Vector<Double> & lCoord)
  {
    ENTER_FCN( "CurlEdgeOp::CalcElemMagVec" );
  
    Error( "CurlEdgeOp::CalcElemMagVec: Not working due to EQN-class",
           __FILE__, __LINE__);
    // UInt nrEdges = ptElement->ptElem->GetNumEdges();
    //   UInt nrNodes = ptElement->ptElem->GetNumNodes();
    //   BaseFE * ptElem = ptElement->ptElem;
    //   ShortInt dim = ptElem->GetDim();

    //   Matrix<Double> cornerCoords; 
    //   Matrix<Double> shape;

    //   StdVector<Double> sol(nrEdges);
    //   // global edge index
    //   StdVector<UInt> epos(nrEdges);
    //   StdVector<UInt> esign(nrEdges);
    //   Vector<UInt> pos(nrNodes);



    //   magVecPot.Resize(dim);

    //   for (UInt i=0; i<dim; i++)
    //     magVecPot[i] = 0;

    //   ptGrid_->GetCoordNodesElemMat(ptElement->connect, cornerCoords, level_);
  
    //   ptElem->CalcEdgeShapeFnc(shape, lCoord, cornerCoords);

  
    //   for (UInt i=0; i < nrNodes; i++)
    //      pos[i] = (*ptMesh2PDENode_)[ptElement->connect[i]-1];
  
    //   ptElem->GetGlobalEdgeIndices(epos, &pos[0], algsys_);


    //   for (UInt j=0; j<nrEdges; j++)
    //     {
    //       esign[j] = epos[j]/abs(epos[j]);
    //       epos[j]  = abs(epos[j]);
    //       sol[j] = (*sol_)[epos[j]-1] * esign[j];
    //     }
  

  
  
    //   // loop over edge curls
    //   // magVecPot = sol * shape;
    //   for( UInt j=0; j<dim; j++ )  
    //     {
    //       magVecPot[j]=0;
    //       for( UInt i=0; i < nrEdges; i++ )
    //      magVecPot[j] += shape[i][j] * sol[i];
    //     }
  }


  //=========================== CurlNodeOperator-Class=======================================

  template<class TYPE>
  CurlNodeOp<TYPE>::CurlNodeOp(Grid * ptGrid, 
                               StdPDE * ptPDE,
                               shared_ptr<EqnMap> eqnMap,
                               NodeStoreSol<TYPE> & aSol,
                               bool coordUpdate) 
    : BaseOperator(ptGrid, ptPDE , eqnMap, false, coordUpdate )
  {
    ENTER_FCN( "CurlNodeOp::CurlNodeOp" );

    sol_ = &aSol;
  }

  template<class TYPE>
  CurlNodeOp<TYPE>::~CurlNodeOp()
  {
    ENTER_FCN( "CurlNodeOp::~CurlNodeOp" );

  }

  template<class TYPE>
  void CurlNodeOp<TYPE>::CalcElemCurlNode(Vector<TYPE> & B, 
                                          const EntityIterator& ent,
                                          const Vector<Double> & lCoord ) 
  {
    ENTER_FCN( "CurlNodeOp::CalcElemCurlNode" );
  
    UInt dim;
    TYPE solEntry;
    dim = ent.GetElem()->ptElem->GetDim();
    Matrix<Double> cornerCoords;
    ptGrid_->GetElemNodesCoord( cornerCoords, ent.GetElem()->connect,
                                coordUpdate_ );
    BaseFE * elem = ent.GetElem()->ptElem;
    if (dim ==2) {
        B.Resize(dim);
        B.Init();

        UInt nShFnc = 0;
        nShFnc = elem->GetNumNodes();

       
        //ptElement->ptElem->SetAnsatzFct( rsult_->fctType );
        
        Matrix<Double> CornerCoords; 
        ptGrid_->GetElemNodesCoord( CornerCoords, ent.GetElem()->connect, 
                                    coordUpdate_ );
      
        Matrix<Double> GlobalGradient;
      
        elem->GetGlobDerivShFnc(GlobalGradient, lCoord, 
                                cornerCoords, ent.GetElem() );
      
        if (isaxi_)
          {
            Vector<Double> ShpFncAtIp;
            Vector<Double> CoordAtIP;
            elem->GetShFnc(ShpFncAtIp,lCoord,ent.GetElem());
            CoordAtIP = CornerCoords * ShpFncAtIp;
            for (UInt i=0; i<nShFnc; i++)
              GlobalGradient[i][0] += ShpFncAtIp[i] / CoordAtIP[0];
          }
      
        // loop over shape functions
        for( UInt i=0; i<dim; i++ )
          for( UInt j=0; j<nShFnc; j++ )
            {
              sol_->Get(ent.GetElem()->connect[j]-1,0,solEntry);
              //const Double solEntry = (*sol_)(ptElement->connect[j],1);
              //B[i] += GlobalGradient[j][i] * (*sol_)[(*ptMesh2PDENode_) [ptElement->connect[j]-1]-1];
              B[i] += GlobalGradient[j][i] * solEntry;
            }
        //account, that we compute the curl!
        TYPE temp = B[0];
        if (isaxi_)
          {
            B[0] = -B[1];
            B[1] = temp;
          }
        else
          {
            B[0] = B[1];
            B[1] = -temp;
          }
      }

    else if (dim ==3)
      Error("CalcElemCurlNode for 3D not implemented",__FILE__,__LINE__);
 
  }

  // explicit teplate instantiation for gcc compiler
#ifdef __GNUC__
  template class CurlNodeOp<Double>;
  template class CurlNodeOp<Complex>;
#endif


} // namespace CoupledField
