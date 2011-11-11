#ifndef FILE_CFS_EDGE_HH
#define FILE_CFS_EDGE_HH

#include <bitset>

#include "General/Environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  struct Edge {

    //! Node numbers defining an edge.
    //! Per definition the first node number
    //! has to be smaller than the second one
    UInt nodes[2];

  };

  struct Face {

    //! Normalize face orientation to match global one and return flagset
    
    //! This method re-orientates a given face to match the global orientation
    //! of the face. It returns the three-entry bitset which contains the
    //! orientation flags as defined in 
    //! Solin,Segeth: Higher-Order Finite Elements, p.165
    //! See also the faceFlags bitset of the class Elem.
    void Normalize( std::bitset<3>& flags); 

    
    //! Return permutation vector for face nodes according to bit-field
    static void GetSortedIndices( StdVector<UInt>& sorted, 
                                  const StdVector<UInt>& unsorted,
                                  UInt numVertices,
                                  const std::bitset<3>& flags );
    
    //! Global node numbers defining a face
    StdVector<UInt> nodes;

  private:

    //! static map for bits of quadrilateral faces
    static char quadBits[4][4];
    
    //! static map for bits of triangular face
    static char triaBits[3][3];
    
  };

  //! external operator for comparing two edges
  bool operator<( const Edge& a, const Edge& b );

  //! external operator for comparing two faces
  bool operator<( const Face& a, const Face& b );

}



#endif
