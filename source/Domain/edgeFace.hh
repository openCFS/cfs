#ifndef FILE_CFS_EDGE_HH
#define FILE_CFS_EDGE_HH

#include <bitset>

#include "General/environment.hh"
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
    //! of the face. It eturns the three-entry bitset which contains the
    //! orientation flags as defined in 
    //! Solin,Segeth: Higher-Order Finite Elements, p.165
    void Normalize( std::bitset<3>& flags); 

    //! Global node numbers definig a face
    StdVector<UInt> nodes;

  private:

    //! static map for bits of quadrilateral faces
    static Char quadBits[4][4];
    
    //! static map for bits of triangular face
    static Char triaBits[3][3];
    
  };

  //! external operator for comparing two edges
  bool operator<( const Edge& a, const Edge& b );

  //! external operator for comparing two faces
  bool operator<( const Face& a, const Face& b );

}



#endif
