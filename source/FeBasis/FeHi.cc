#include "FeHi.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/ElemMapping/EdgeFace.hh"

namespace CoupledField {
  
// declare class specific logging stream
DEFINE_LOG(feHi, "feHi")


  FeHi::FeHi(Elem::FEType feType  ) {
    myFeType_ = feType;
    elemShape_ = Elem::shapes[feType];

    updateUnknowns_ = true;
    isoOrder_ = 0;
    maxOrder_ = 0;
    isIsotropic_ = false;  // Must initialize to false so SetIsoOrder() doesn't early-exit before orderEdge_ is set
  }

  FeHi::FeHi(const FeHi& other){
    //tremendous amount of copies...
    this->myFeType_ = other.myFeType_;
    this->elemShape_ = Elem::shapes[this->myFeType_];
    this->updateUnknowns_ = other.updateUnknowns_;
    this->isoOrder_ = other.isoOrder_;
    this->maxOrder_ = other.maxOrder_;
    this->isIsotropic_ = other.isIsotropic_;
    this->entityFncs_ = other.entityFncs_;
    this->isoOrder_ = other.isoOrder_;
    this->anisoOrder_ = other.anisoOrder_;
    this->orderEdge_ = other.orderEdge_;
    this->orderFace_ = other.orderFace_;
    this->orderInner_ = other.orderInner_;
  }
  
  FeHi::~FeHi() {
    
  }
  
  
  void FeHi::SetIsoOrder( UInt order ) {
     
     LOG_DBG3(feHi) << "SetIsoOrder " << order
         << " for H1Hi elem of type "
         << Elem::feType.ToString(myFeType_);
     // just change, if order is different from previously set one
     if( order == this->isoOrder_ && isIsotropic_ ) return;

     orderEdge_.Resize(elemShape_.numEdges);
     orderFace_.Resize(elemShape_.numFaces);
     maxOrder_ = 0;
     
     // set order for edges
     orderEdge_.Init(order);
     
     // set order for faces
     boost::array<UInt,2> faceOrder = {{order, order}};
     orderFace_.Init(faceOrder);

     // set order for inner
     boost::array<UInt, 3> innerOrder = {{order, order, order}}; 
     orderInner_ = innerOrder;
     
     maxOrder_ = std::max( maxOrder_, order);
     updateUnknowns_ = true;
     isIsotropic_ = true;
     isoOrder_ = order;
   }
   
   
   void FeHi::SetAnisoOrder( const StdVector<UInt>& order,
                               const Elem* ptElem ) {
     LOG_DBG3(feHi) << "SetAnisoOrder for H1Hi elem #" << ptElem->elemNum
             << " of type " 
             << Elem::feType.ToString(myFeType_)
             << " to order " << order.Serialize();
     
     // check initially, if element support anisotropic 
     // polynomial order
     // -> this is performed later on
     
     anisoOrder_ = order;
     
     // resize edge and face arrays
     orderEdge_.Resize(elemShape_.numEdges);
     orderFace_.Resize(elemShape_.numFaces);
     orderEdge_.Init();
     orderFace_.Init();
     // -------
     //  EDGES
     // -------
     Integer locDir = -1;
     for( UInt iEdge = 0; iEdge < elemShape_.numEdges; ++iEdge ) {
       LOG_DBG3(feHi) << "\tTreating edge #" << std::abs(ptElem->extended->edges[iEdge]) << std::endl;
       locDir = elemShape_.edgeLocDirs[iEdge];
       SetEdgeOrder( iEdge, order[locDir]);
       LOG_DBG3(feHi) << "\t\tdir: #" << ptElem->connect[elemShape_.edgeNodes[iEdge][0]-1]
                        << " -> #" << ptElem->connect[elemShape_.edgeNodes[iEdge][1]-1]
                        << ", order: " << order[locDir];
       
     }
     
     // -------
     //  FACES
     // -------
     if( elemShape_.numFaces > 0 ) {
       
       // variables for face-local directions
       UInt locDir1 = 0, locDir2 = 0;
       boost::array<UInt,2> faceOrder;
       for( UInt iFace = 0; iFace < elemShape_.numFaces; ++iFace ) {
         LOG_DBG3(feHi) << "\tTreating face #" << ptElem->extended->faces[iFace] << std::endl;
         
         const StdVector<UInt>& unsorted = elemShape_.faceNodes[iFace];
         
         // check: currently, we can only set an anisotropic order,
         // if the face has 4 nodes
         if( unsorted.GetSize() != 4 ) {
           EXCEPTION( "Anisotropic order only supported for quad-faces");
         }
         
         // use the face flags of the current face to match the element-
         // local directions (xi,eta,zeta) to the face local ones.
         // Note: This only works for quad-shaped faces
         locDir1 = elemShape_.faceLocDirs[iFace][0];
         locDir2 = elemShape_.faceLocDirs[iFace][1];
          // check, if directions have to get interchanged
         if( !ptElem->extended->faceFlags[iFace][2]) {
           std::swap( locDir1, locDir2 );
         }
         faceOrder[0] = order[locDir1];
         faceOrder[1] = order[locDir2];
         
         // Logging stuff
         StdVector<UInt> ind;
         Face::GetSortedIndices( ind, unsorted,  4, ptElem->extended->faceFlags[iFace]);
         LOG_DBG3(feHi) << "\t\tdir1: node #"
                          << ptElem->connect[ind[1]] 
                          << "-> #" << ptElem->connect[ind[0]] 
                          << ", order: " << faceOrder[0]; 
         LOG_DBG3(feHi) << "\t\tdir2: node #" 
                          << ptElem->connect[ind[3]] 
                          << "-> #" << ptElem->connect[ind[0]] 
                          << ", order: " << faceOrder[1];
         
         SetFaceOrder( iFace, faceOrder );
       }
     }
     // -------
     //  INNER
     // -------
     if( elemShape_.dim == 3 ) {
       LOG_DBG3(feHi) << "\tTreating interior" << std::endl;
       LOG_DBG3(feHi) << "\t\tSetting Order to "  << order.Serialize() 
               << std::endl;
       boost::array<UInt,3> orderInner;
       orderInner[0] = order[0];
       orderInner[1] = order[1];
       orderInner[2] = order[2];
       SetInteriorOrder( orderInner );
     }
     
     // set status to anisotropic and determine maximal order
     isIsotropic_ = false;
     maxOrder_ = *(std::max_element(order.Begin(), order.End()));
     updateUnknowns_ = true;
   }
   
   void FeHi::GetNodesExceedingOrder( std::set<UInt>& nodes, 
                                        const StdVector<UInt>& order,
                                        const Elem* ptElem ) {

     Matrix<UInt> fullOrder;
     GetPolyOrderOfNodes( fullOrder, ptElem );

     UInt numFncs = fullOrder.GetNumRows();
     UInt dim = fullOrder.GetNumCols();
     nodes.clear();

     // loop over all function
     for( UInt iFnc = 0; iFnc < numFncs; ++iFnc ) {
       // loop over all directions
       for( UInt iDim = 0; iDim < dim; ++iDim ) {
         // check, if component exceeds the maximum allowed order
         if( fullOrder[iFnc][iDim] > order[iDim] ) {
           nodes.insert(iFnc);
           break;
         }
       }
     }
   }
     
   void FeHi::SetEdgeOrder( UInt edgeNum, UInt order ) {
     assert( edgeNum <= elemShape_.numEdges);

     if( orderEdge_[edgeNum] == order) return;
     
     LOG_DBG3(feHi) << "\tSetting order of edge #" << edgeNum << " to " << order;
     orderEdge_[edgeNum] = order;
     maxOrder_ = std::max( maxOrder_, order);
     isIsotropic_ = false;
     updateUnknowns_ = true;
   }
   
   void FeHi::SetFaceOrder( UInt faceNum,
                              const boost::array<UInt,2>& order ) {
     assert( faceNum <= elemShape_.numFaces);

     if( orderFace_[faceNum] == order ) return;
     
     orderFace_[faceNum] = order;
     maxOrder_ = std::max( maxOrder_, 
                           std::max( order[0], order[1]) );
     isIsotropic_ = false;
     updateUnknowns_ = true;
   }

   void FeHi::SetInteriorOrder( const boost::array<UInt,3>& order ) {

     if( orderInner_ == order ) return;
     orderInner_ = order;
     maxOrder_ = std::max( maxOrder_, 
                           *std::max_element( order.begin(), order.end() ) );
     isIsotropic_ = false;
     updateUnknowns_ = true;
   }
} // end of namespace
